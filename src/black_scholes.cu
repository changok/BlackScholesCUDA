#include <iostream>
#include <limits>
#include <cstddef>
#include <cassert>
#include <cmath>

#include <stdio.h>

#include <cuda.h>
#include <curand_kernel.h>

#include "black_scholes.cuh"

using namespace std;

const int WINDOW_WIDTH = 128;

__global__ void setup_rnd_kernel ( curandState * state, unsigned long seed )
{
    long id = (blockIdx.x * WINDOW_WIDTH) + threadIdx.x;
    curand_init ( seed, id, 0, &state[id] );
} 

__device__ double black_scholes_value (const double S,
             const double E, const double r, const double sigma,
             const double T, const double random_number) {
    const double current_value = S * exp ( (r - (sigma*sigma) / 2.0) * T + 
                       sigma * sqrt (T) * random_number );
    return exp (-r * T) * 
      ((current_value - E < 0.0) ? 0.0 : current_value - E);
}

__global__ void black_scholes_kernel(const double S, const double E, 
            const double r, const double sigma, const double T,
            const long M, double* blockMeans, double* cudaTrials,
            curandState* rnds, double* prns) {
    
    __shared__ double sum_of_trials[WINDOW_WIDTH];
    
    unsigned int gId = (blockIdx.x * WINDOW_WIDTH) + threadIdx.x; 
    unsigned int tId = threadIdx.x;

    // Do the Black-Scholes iterations
    //const double random_number = 1.0; 
    curandState localState = rnds[gId];
    double RANDOM = curand_uniform( &localState );

//test
prns[gId] = RANDOM;

    rnds[gId] = localState;

    double value = black_scholes_value (S, E, r, sigma, T, RANDOM);
    cudaTrials[gId] = value;
    
    // we need to keep origianl trial values for calculatng standard deviation
    sum_of_trials[tId] = value;

    for(unsigned int stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
        __syncthreads();
        sum_of_trials[tId] += sum_of_trials[tId + stride];
    }

    // Pack the OUT values into the args struct 
    if(tId == 0) {
        blockMeans[blockIdx.x] = sum_of_trials[0]/(double)M;
    }
}

__global__ void black_scholes_variance_kernel(const double mean,
            const long M, double* cudaTrials, double* cudaVariances) {
    
    __shared__ double variances[WINDOW_WIDTH];
    
    unsigned int gId = (blockIdx.x * WINDOW_WIDTH) + threadIdx.x; 
    unsigned int tId = threadIdx.x;
    
    variances[tId] = cudaTrials[gId];
    variances[tId] = variances[tId] - mean;
    variances[tId] = variances[tId] *  variances[tId] / (double)M;

    for(unsigned int stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
        __syncthreads();
        variances[tId] += variances[tId + stride];
    }

    if(tId == 0) {
        cudaVariances[blockIdx.x] = variances[0];
    }
}

cit black_scholes(const double S, const double E, const double r,
                   const double sigma, const double T, const long M) {
    cit interval;
    int num_of_blocks = M/WINDOW_WIDTH;
    double* means = new double[num_of_blocks];
    double stddev = 0.0;
    double conf_width = 0.0;

    assert (M > 0);
    double* trials = new double[M]; //Array containing the results of each of the M trials.
    assert (trials != NULL);
    long size = M * sizeof(double);

    dim3 dimGrid(num_of_blocks);
    dim3 dimBlock(WINDOW_WIDTH);

    curandState* devStates;
    cudaMalloc((void **) &devStates, M * sizeof(curandState));
    setup_rnd_kernel<<<dimGrid, dimBlock>>>(devStates, time(NULL));

    double* blockMeans;
    cudaMalloc((void**) &blockMeans, num_of_blocks * sizeof(double));
    
    double* cudaTrials;
    cudaMalloc((void**) &cudaTrials, size);

//test
double* prns;
cudaMalloc((void**) &prns, size);
double* prns_table = new double[size];

    black_scholes_kernel<<<dimGrid, dimBlock>>>(S, E, r, sigma, T, M, blockMeans, cudaTrials, devStates, prns);
    
    cudaMemcpy(means, blockMeans, num_of_blocks * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(blockMeans);

//test
cudaMemcpy(prns_table, prns, size, cudaMemcpyDeviceToHost);
for (long i = 0; i < M; i++) {
	printf("%lf, ", prns_table[i]);
	//cout << prns_table[i] << ", ";
}
cout << endl;
    
    double mean = 0.0;
    // combine results from each threads
    for (long i = 0; i < num_of_blocks; i++) {
        mean += means[i];
    }
    
    stddev = black_scholes_stddev (mean, M, cudaTrials);
    cudaMemcpy(trials, cudaTrials, size, cudaMemcpyDeviceToHost);
    cudaFree(cudaTrials);

    conf_width = 1.96 * stddev / sqrt ((double) M);
    interval.min = mean - conf_width;
    interval.max = mean + conf_width;

    cudaFree(devStates);
    delete [] trials;
    delete [] means;
    
    return interval;
}

double black_scholes_stddev (const double mean, const long M, double* cudaTrials) {
    const long numOfBlocks = M/WINDOW_WIDTH;
    double* variances = new double[numOfBlocks]; 
    double* cudaVariances;
    cudaMalloc((void**) &cudaVariances, numOfBlocks * sizeof(double));
    
    dim3 dimGrid(numOfBlocks);
    dim3 dimBlock(WINDOW_WIDTH);
    black_scholes_variance_kernel<<<dimGrid, dimBlock>>>(mean, M, cudaTrials, cudaVariances);
    cudaMemcpy(variances, cudaVariances, numOfBlocks * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(cudaVariances);
    
    double variance = 0.0;
    for(long idx=0; idx<numOfBlocks; idx++) {
        variance += variances[idx];
    }
    
    return sqrt(variance);
}

/**
 * Compute the standard deviation of trials[0 .. M-1].
 */
/*
double black_scholes_stddev (const double mean, const long M, double* trials) {
    double variance = 0.0;
    long k;
    
    for (k = 0; k < M; k++) {
        const double diff = trials[k] - mean;
        
        // Just like when computing the mean, we scale each term of this
        // sum in order to avoid overflow.
        //
        variance += diff * diff / (double) M;
    }
    
    return sqrt (variance);
}
*/
