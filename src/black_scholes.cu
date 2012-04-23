#include <cstddef>
#include <cassert>
#include <cmath>

#include "cuda.h"
#include "black_scholes.h"

const int TILE_WIDTH = 128;

__device__ static inline double black_scholes_value (const double S,
             const double E, const double r, const double sigma,
             const double T, const double random_number) {
    const double current_value = S * exp ( (r - (sigma*sigma) / 2.0) * T + 
                       sigma * sqrt (T) * random_number );
    return exp (-r * T) * 
      ((current_value - E < 0.0) ? 0.0 : current_value - E);
}

__global__ void black_scholes_kernel(const double S, const double E, 
            const double r, const double sigma, const double T,
            const long M, double* blockMeans, double* cudaTrials ) {
    
    __shared__ double trials[TILE_WIDTH];
    __shared__ double sum_of_trials[TILE_WIDTH];
    
    const unsigned int tId = threadIdx.x;

    /* Do the Black-Scholes iterations */
    // M/nthreads
    const double random_number = 1.0; 
    trials[tId] = black_scholes_value (S, E, r, sigma, T, random_number);
    
    // we need origianl trial values to calculate standard deviation
    sum_of_trials[tId] = trials[tId];

    /*
     * We scale each term of the sum in order to avoid overflow. 
     * This ensures that mean is never larger than the max
     * element of trials[0 .. M-1].
     */
    for(unsigned int stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
        __syncthreads();
        sum_of_trials[tId] += sum_of_trials[tId + stride];
    }

    /* Pack the OUT values into the args struct */
    if(tId == 0) {
        blockMeans[blockIdx.x] = sum_of_trials[0]/(double)M;
    }
}

cit black_scholes(const double S, const double E, const double r,
                   const double sigma, const double T, const long M) {
    cit interval;
    long num_of_blocks = M/TILE_WIDTH;
    double* means = new double[num_of_blocks];
    double stddev = 0.0;
    double conf_width = 0.0;

    assert (M > 0);
    double* trials = new double[M]; //Array containing the results of each of the M trials.
    long size = M * sizeof(double);
    assert (trials != NULL);

    double* blockMeans;
    cudaMalloc((void**) &blockMeans, size);
    
    double* cudaTrials;
    cudaMalloc((void**) &cudaTrials, size);

    dim3 dimGrid(num_of_blocks);
    dim3 dimBlock(TILE_WIDTH);

    black_scholes_kernel<<<dimGrid, dimBlock>>>(S, E, r, sigma, T, M, blockMeans, cudaTrials);
    //black_scholes_kernel (S, E, r, sigma, T, M, cudaMean, cudaTrials);
    
    cudaMemcpy(&means, &blockMeans, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(blockMeans);
    
    cudaMemcpy(&trials, &cudaTrials, size, cudaMemcpyDeviceToHost);
    cudaFree(cudaTrials);
    
    double mean = 0.0;
    // combine results from each threads
    for (long i = 0; i < num_of_blocks; i++) {
        mean += means[i];
    }
    
    stddev = black_scholes_stddev (mean, M, trials);

    conf_width = 1.96 * stddev / sqrt ((double) M);
    interval.min = mean - conf_width;
    interval.max = mean + conf_width;

    delete [] trials;
    delete [] means;
    
    return interval;
}

/**
 * Compute the standard deviation of trials[0 .. M-1].
 */
static double black_scholes_stddev (const double mean, const long M, const double* trials) {
    double variance = 0.0;
    long k;
    
    for (k = 0; k < M; k++) {
        const double diff = trials[k] - mean;
        /*
        * Just like when computing the mean, we scale each term of this
        * sum in order to avoid overflow.
        */
        variance += diff * diff / (double) M;
    }
    
    return sqrt (variance);
}