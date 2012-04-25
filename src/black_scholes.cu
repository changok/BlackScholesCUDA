#include <iostream>
#include <cstddef>
#include <cassert>
#include <cmath>

#include <stdio.h>

#include <cuda.h>
#include <curand_kernel.h>

#include "black_scholes.cuh"

using namespace std;

const int WINDOW_WIDTH = 128;

//__global__ void setup_rnd_kernel ( curandState * state, unsigned long seed )
__global__ void setup_rnd_kernel ( curandState * state, time_t seed )
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

__device__ gaussrand_result_t gaussrand (curandState* localState) {
  gaussrand_result_t result;

  //if (u1==u2) return -1.0;

  double v1, v2, s;
  do {
    v1 = 2.0 * curand_uniform(localState) - 1.0;
    v2 = 2.0 * curand_uniform(localState) - 1.0;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1 || s== 0);

  double w = sqrt ( (-2.0 * log (s)) / s);
  result.grand1 = v1 * w;
  result.grand2 = v2 * w;

  return result;
}

__global__ void black_scholes_kernel(const double S, const double E, 
            const double r, const double sigma, const double T,
            const long M, double* blockMeans, double* cudaTrials,
            curandState* randStates, double* debug) {
    
    __shared__ double sum_of_trials[WINDOW_WIDTH];
   
    /* Since we decreased ACTUAL size of block to the half
	   , gId which is thread position in KERNEL will be jumped
	   by WINDOW_WIDTH/2 at each end of tId in a block */
    unsigned int gId = (blockIdx.x * WINDOW_WIDTH) + (threadIdx.x); 
    unsigned int tId = threadIdx.x;

    // Do the Black-Scholes iterations
    //const double random_number = 1.0; 
    curandState localState = randStates[gId];
    double RANDOM = curand_uniform( &localState );
    //double unirand1 = curand_uniform( &localState );
    //double unirand2 = curand_uniform( &localState );
    gaussrand_result_t gresult = gaussrand (&localState);
// prng test
//debug[gId] = grand1;
debug[gId] = gresult.grand1;
//debug[gId+(WINDOW_WIDTH/2)] = gresult.grand2;
// prng test end
    randStates[gId] = localState;

    //double value = black_scholes_value (S, E, r, sigma, T, RANDOM);
    //cudaTrials[gId] = value;
    gresult.grand1 = black_scholes_value (S, E, r, sigma, T, gresult.grand1);
    cudaTrials[gId] = gresult.grand1;
//    gresult.grand2 = black_scholes_value (S, E, r, sigma, T, gresult.grand2);
    //cudaTrials[gId+(WINDOW_WIDTH/2)] = gresult.grand2;
    
    // we need to keep origianl trial values for calculatng standard deviation
    //sum_of_trials[tId] = value;
    sum_of_trials[tId] = gresult.grand1;
//    sum_of_trials[tId+(WINDOW_WIDTH/2)] = gresult.grand2;

// shouldn't it be WINDOW_WIDTH instead of blockDim.x?

	// don't have to do "blockDim.x >> 1" because of faked size of block
 //   for(unsigned int stride = blockDim.x; stride > 0; stride >>= 1) {
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
    variances[tId] = (variances[tId] *  variances[tId] )/ (double)M;

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
    int num_of_blocks = max((long)1,M/WINDOW_WIDTH);
    double* means = new double[num_of_blocks];
    double stddev = 0.0;
    double conf_width = 0.0;

    assert (M > 0);
    double* trials = new double[M]; //Array containing the results of each of the M trials.
    long size = M * sizeof(double);
    assert (trials != NULL);

    dim3 dimGrid(num_of_blocks);
    
    /* The number of threads in each block was supposed to be WINDOW_WIDTH(128). 
       However, the feature of gaussian randum number generator forced us to make 
       work per each thread to be twice of black scholes value operation instead of one. 
       So, we decide to decrease the ACTUAL number of thread(WINDOW_WIDTH) by the half of
       , and pretend to operate with WINDOW_WIDTH number of threads in a block. */
    //dim3 dimBlock(WINDOW_WIDTH/2);
    dim3 dimBlock(WINDOW_WIDTH);
	////////////////////////////////////////
	// below pretend working with MOCKED thread number

    curandState* randStates;
    cudaMalloc((void **) &randStates, M * sizeof(curandState));
    setup_rnd_kernel<<<dimGrid, dimBlock>>>(randStates, time(NULL));

    double* blockMeans;
    cudaMalloc((void**) &blockMeans, num_of_blocks * sizeof(double));
    
    double* cudaTrials;
    cudaMalloc((void**) &cudaTrials, size);

// prng test
    double* hostRands = new double[M];
    double* cudaRands;
    cudaMalloc((void**) &cudaRands, size);
    cudaMemset(cudaRands, 0, size);
// prng test end

    black_scholes_kernel<<<dimGrid, dimBlock>>>(S, E, r, sigma, T, M, blockMeans, cudaTrials, randStates, cudaRands);
    
    cudaMemcpy(means, blockMeans, num_of_blocks * sizeof(double), cudaMemcpyDeviceToHost);

// prng test
    cudaMemcpy(hostRands, cudaRands, size, cudaMemcpyDeviceToHost);
    for (int i = 0; i < M; i++) 
      printf("%lf, ", hostRands[i]);
    cout << endl;
// prng test end
    
    double mean = 0.0;
    // combine results from each threads
    for (long i = 0; i < num_of_blocks; i++) {
        mean += means[i];
    }
    
    stddev = black_scholes_stddev (mean, M, cudaTrials);
    cudaMemcpy(trials, cudaTrials, size, cudaMemcpyDeviceToHost);

	/* confidence interval */
    conf_width = 1.96 * stddev / sqrt ((double) M);
    interval.min = mean - conf_width;
    interval.max = mean + conf_width;

	/* clean up */
	cudaFree(cudaRands);
    cudaFree(cudaTrials);
    cudaFree(blockMeans);
    cudaFree(randStates);
    delete [] trials;
    delete [] means;
    
    return interval;
}

double black_scholes_stddev (const double mean, const long M, double* cudaTrials) {
    const long num_of_blocks = max((long)1,M/WINDOW_WIDTH);
    double* variances = new double[num_of_blocks]; 
    double* cudaVariances;
    cudaMalloc((void**) &cudaVariances, num_of_blocks * sizeof(double));
    
    dim3 dimGrid(num_of_blocks);
    dim3 dimBlock(WINDOW_WIDTH);
    black_scholes_variance_kernel<<<dimGrid, dimBlock>>>(mean, M, cudaTrials, cudaVariances);
    cudaMemcpy(variances, cudaVariances, num_of_blocks * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(cudaVariances);
    
    double variance = 0.0;
    for(long idx=0; idx<num_of_blocks; idx++) {
        variance += variances[idx];
    }

	/* clean up */
	cudaFree(cudaVariances);
    
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
