#include <iostream>
#include <cstddef>
#include <cassert>
#include <cmath>

#include <stdio.h>

#include <cuda.h>
#include <curand_kernel.h>
#include "cutil.h"
#include "cutil_inline_runtime.h"

#include "common.h"
#include "black_scholes.cuh"
#include "timer.h"
#include "bsconfig.h"

using namespace std;

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

// standard normal distributed random number [0~1]
__device__ gaussrand_result_t gaussrand (curandState* localState) {
  gaussrand_result_t result;

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

__global__ void black_scholes_kernel(double* blockMeans, double* cudaTrials,
            curandState* randStates, const double* fixedRands, double* debug, BSConfig config) {
    
    __shared__ double means[WINDOW_WIDTH];
   
    const long NUM_OF_TOT_THREAD = gridDim.x * blockDim.x;
    const long LOOP_SIZE = config.M / NUM_OF_TOT_THREAD;
    const unsigned int GID = (blockIdx.x * blockDim.x) * LOOP_SIZE + threadIdx.x * LOOP_SIZE;
    const unsigned int TID = threadIdx.x;

	
    curandState localState = randStates[(blockIdx.x * blockDim.x) + threadIdx.x];
    gaussrand_result_t gresult;

    means[TID] = 0.0;

    // Do the Black-Scholes iterations
    for(long trial = 0; trial < LOOP_SIZE; trial++) {
        double value = 0.0;
        if (trial%2 == 0) {
            if (config.RND_MODE == 1) {
                gresult.grand1 = 1.0;
                gresult.grand2 = 1.0;
            }
            // use pre-generated random number
            else if (config.RND_MODE == 2) {
                gresult.grand1 = fixedRands[GID + trial];
                gresult.grand2 = fixedRands[GID + trial+1];
            }
            // use gaussian random number (standard normal distributed)
            else {
                gresult = gaussrand (&localState);
            }

            if(config.DEBUG_LEVEL == 1) {
                debug[GID + trial] = gresult.grand1;
            }
            value = black_scholes_value (config.S, config.E, config.r, config.sigma, config.T, gresult.grand1);
        } else {
            if(config.DEBUG_LEVEL == 1) {
                debug[GID + trial] = gresult.grand2;
            }
            value = black_scholes_value (config.S, config.E, config.r, config.sigma, config.T, gresult.grand2);
        }

        // we need to keep origianl trial values for calculatng standard deviation
        // for current calculation, we use trials
        // Also, to prevent overflow caused by adding, divide the value by M in advance
        means[TID] += value/config.M;
        cudaTrials[GID + trial] = value;
    }

    for(unsigned int stride = blockDim.x>>1; stride > 0; stride >>= 1) {
        __syncthreads();
		if (TID < stride)
		    means[TID] += means[TID + stride];
    }

    if(TID == 0) {
        blockMeans[blockIdx.x] = means[0];
    }
}

__device__ void trunc(double* target) {
    if (*target < 0.0000000005 && *target > 0)
        *target = 0.0;
    if (*target > -0.0000000005 && *target < 0)
        *target = 0.0;
}

__global__ void black_scholes_variance_kernel(const double mean,
            const long M, double* cudaTrials, double* cudaVariances) {
    
    __shared__ double variances[WINDOW_WIDTH];
    
    unsigned int gId = (blockIdx.x * WINDOW_WIDTH) + threadIdx.x; 
    unsigned int tId = threadIdx.x;
    
    variances[tId] = cudaTrials[gId];
    variances[tId] = variances[tId] - mean;

    // Meaningless value such as 1.1E-15 could lead invalid result
    // when number of trial is so high. Thus, truncate all after the 10th
    // decimal place. Even though we truncate them, the result still in
    // acceptable valid range
    trunc(&variances[tId]);

    variances[tId] = (variances[tId] *  variances[tId])/ (double)(M-1);

    for(unsigned int stride = blockDim.x>>1; stride > 0; stride >>= 1) {
        __syncthreads();
		if (stride > tId)
        	variances[tId] += variances[tId + stride];
    }

    if(tId == 0) {
        cudaVariances[blockIdx.x] = variances[0];
    }
}

Result black_scholes(double* cudafixedRands, BSConfig config) {
    Result result;

    double* means = new double[config.totalNumOfBlocks()];
    double conf_width = 0.0;
	double t1, t2;

    assert (config.M > 0);
    long size = config.M * sizeof(double);

    dim3 dimGrid(config.totalNumOfBlocks());
    dim3 dimBlock(WINDOW_WIDTH);

    // part5_start
	t1 = 0; t1 = get_seconds();
    curandState* randStates;
    cutilSafeCall(cudaMalloc((void **) &randStates, config.totalNumOfThread() * sizeof(curandState)));

    setup_rnd_kernel<<<dimGrid, dimBlock>>>(randStates, time(NULL));
	t2 = 0; t2 = get_seconds();
	result.init_seeds_setup_time = t2 - t1;
	// part5_end

	// part3_begin
	t1 = 0; t1 = get_seconds();

    double* blockMeans;
    cutilSafeCall(cudaMalloc((void**) &blockMeans, config.totalNumOfBlocks() * sizeof(double)));
    
    double* cudaTrials;
    cutilSafeCall(cudaMalloc((void**) &cudaTrials, size));

    double* hostDebug = NULL;
    double* cudaDebug;
//    if (config.DEBUG_LEVEL == 1) {
        hostDebug = new double[config.M];
        cutilSafeCall(cudaMalloc((void**) &cudaDebug, size));
//    }

    black_scholes_kernel<<<dimGrid, dimBlock>>>(blockMeans, cudaTrials, randStates, cudafixedRands, cudaDebug, config);

    cutilSafeCall(cudaMemcpy(means, blockMeans, config.totalNumOfBlocks() * sizeof(double), cudaMemcpyDeviceToHost));

    if (config.DEBUG_LEVEL == 2) {
        cudaMemcpy(hostDebug, cudaDebug, size, cudaMemcpyDeviceToHost);
        for (int i = 0; i < config.M; i++) {
            printf("r%d: %lf, ", i, hostDebug[i]);
        }
        puts("\n");

        for (int i = 0; i < config.totalNumOfBlocks(); i++) {
            printf("m%d: %lf, ", i, means[i]);
        }
        puts("");

        double* t = new double[config.M];
        cutilSafeCall(cudaMemcpy(t, cudaTrials, size, cudaMemcpyDeviceToHost));
        for (int i = 0; i < config.M; i++) {
            printf("t%d: %lf, ", i, t[i]);
        }
        puts("");

        delete [] t;
    }

	t2 =0; t2 = get_seconds();
	result.black_sholes_kernel_time = t2 - t1;
	// part3_end

	// part4_begin
    t1 = 0;
	t1 = get_seconds();
    result.mean = 0.0;

    // combine results from each blocks
    for (long i = 0; i < config.totalNumOfBlocks(); i++) {
        result.mean += means[i];
    }

	result.stddev = black_scholes_stddev (result.mean, config.M, cudaTrials);

	t2 = 0;
	t2 = get_seconds();
	result.calc_stddev_time = t2 - t1;
	// part4_end

    // confidence interval
    conf_width = 1.96 * result.stddev / sqrt ((double) config.M);
    result.min = result.mean - conf_width;
    result.max = result.mean + conf_width;

    /* clean up */
    cudaFree(cudaDebug);
    cudaFree(cudaTrials);
    cudaFree(blockMeans);
    cudaFree(randStates);

    if(hostDebug != NULL) delete [] hostDebug;
    if(means != NULL) delete [] means;

    return result;
}

double black_scholes_stddev (const double mean, const long M, double* cudaTrials) {
    double* variances = new double[M/WINDOW_WIDTH];
    double* cudaVariances;
    cutilSafeCall(cudaMalloc((void**) &cudaVariances, M/WINDOW_WIDTH * sizeof(double)));

    dim3 dimGrid(M/WINDOW_WIDTH);
    dim3 dimBlock(WINDOW_WIDTH);
    black_scholes_variance_kernel<<<dimGrid, dimBlock>>>(mean, M, cudaTrials, cudaVariances);

    cutilSafeCall(cudaMemcpy(variances, cudaVariances, M/WINDOW_WIDTH * sizeof(double), cudaMemcpyDeviceToHost));
    double variance = 0.0;
    for(long idx=0; idx<M/WINDOW_WIDTH; idx++) {
        variance += variances[idx];
    }
    cout << endl;

    cudaFree(cudaVariances);
    delete [] variances;
    
    return sqrt(variance);
}
