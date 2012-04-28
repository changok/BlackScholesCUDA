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

__device__ double black_scholes_value (BSConfig cf, const double random_number) {
    const double current_value = cf.S * exp ( (cf.r - (cf.sigma*cf.sigma) / 2.0) * cf.T +
                                 cf.sigma * sqrt (cf.T) * random_number );
    return exp (-cf.r * cf.T) *
      ((current_value - cf.E < 0.0) ? 0.0 : current_value - cf.E);
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
   
    const long LOOP_SIZE = config.M / (BLOCK_SIZE * WINDOW_WIDTH);
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

                if(config.DEBUG_LEVEL == 2) {
                    debug[GID + trial] = gresult.grand1;
                    debug[GID + trial + 1] = gresult.grand2;
                }
            }
            // use gaussian random number (standard normal distributed)
            else {
                gresult = gaussrand (&localState);
            }

            value = black_scholes_value (config, gresult.grand1);
        } else {
            value = black_scholes_value (config, gresult.grand2);
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
    else if (*target > -0.0000000005 && *target < 0)
        *target = 0.0;
}

__global__ void black_scholes_variance_kernel(const long M, const double mean,
            double* cudaTrials, double* cudaVariances, double* debug) {
    
    __shared__ double variances[WINDOW_WIDTH];
    
    long gId = (blockIdx.x * WINDOW_WIDTH) + threadIdx.x;
    unsigned int tId = threadIdx.x;
    
    variances[tId] = cudaTrials[gId];
    variances[tId] = variances[tId] - mean;

    // Meaningless value such as 1.1E-15 could lead invalid result
    // when number of trial is so high. Thus, truncate all after the 10th
    // decimal place. Even though we truncate them, the result still in
    // acceptable valid range
    trunc(&variances[tId]);

    variances[tId] = (variances[tId] *  variances[tId]) / (double)(M-1);
    debug[gId] = cudaTrials[gId];

    for(unsigned int stride = WINDOW_WIDTH>>1; stride > 0; stride >>= 1) {
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
	t1 = get_seconds();
    curandState* randStates;
    cutilSafeCall(cudaMalloc((void **) &randStates, config.totalNumOfThread() * sizeof(curandState)));

    setup_rnd_kernel<<<dimGrid, dimBlock>>>(randStates, time(NULL));
	t2 = get_seconds();
	result.init_seeds_setup_time = t2 - t1;
	// part5_end

	// part3_begin
	t1 = get_seconds();

    double* blockMeans;
    cutilSafeCall(cudaMalloc((void**) &blockMeans, config.totalNumOfBlocks() * sizeof(double)));
    
    double* cudaTrials;
    cutilSafeCall(cudaMalloc((void**) &cudaTrials, size));

    double* hostDebug = new double[config.M];
    double* cudaDebug;
    cutilSafeCall(cudaMalloc((void**) &cudaDebug, size));

    black_scholes_kernel<<<dimGrid, dimBlock>>>(blockMeans, cudaTrials, randStates, cudafixedRands, cudaDebug, config);

    cutilSafeCall(cudaMemcpy(means, blockMeans, config.totalNumOfBlocks() * sizeof(double), cudaMemcpyDeviceToHost));

    if (config.DEBUG_LEVEL == 2) {
        cudaMemcpy(hostDebug, cudaDebug, size, cudaMemcpyDeviceToHost);
        for (int i = 0; i < config.M; i++) {
            if(i < 10 || i > (config.M - 10))
                printf("RND[%d]: %lf\n", i, hostDebug[i]);
        }
        puts("\n");

        for (int i = 0; i < config.totalNumOfBlocks(); i++) {
            if(i < 10 || i > (config.M - 10))
                printf("MEAN[%d]: %lf\n", i, means[i]);
        }
        puts("");

        double* t = new double[config.M];
        cutilSafeCall(cudaMemcpy(t, cudaTrials, size, cudaMemcpyDeviceToHost));
        for (int i = 0; i < config.M; i++) {
            if(i < 10 || i > (config.M - 10))
                printf("TRIAL[%d]: %lf\n", i, t[i]);
        }
        puts("");

        delete [] t;
    }

	t2 = get_seconds();
	result.black_sholes_kernel_time = t2 - t1;
	// part3_end

	// part4_begin
	t1 = get_seconds();
    result.mean = 0.0;

    // combine results from each blocks
    for (long i = 0; i < config.totalNumOfBlocks(); i++) {
        result.mean += means[i];
    }

	result.stddev = black_scholes_stddev(result.mean, config, cudaTrials);

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

double black_scholes_stddev (const double mean, BSConfig config, double* cudaTrials) {
    double* variances = new double[config.M/WINDOW_WIDTH];
    double* cudaVariances;
    cutilSafeCall(cudaMalloc((void**) &cudaVariances, (config.M/WINDOW_WIDTH) * sizeof(double)));

    dim3 dimGrid(config.M/WINDOW_WIDTH);
    dim3 dimBlock(WINDOW_WIDTH);

    double* debug = new double[config.M];

    double* cudaDebug;
    cutilSafeCall(cudaMalloc((void**) &cudaDebug, config.M * sizeof(double)));
    black_scholes_variance_kernel<<<dimGrid, dimBlock>>>(config.M, mean, cudaTrials, cudaVariances, cudaDebug);

    cutilSafeCall(cudaMemcpy(variances, cudaVariances, (config.M/WINDOW_WIDTH) * sizeof(double), cudaMemcpyDeviceToHost));
    cutilSafeCall(cudaMemcpy(debug, cudaDebug, config.M * sizeof(double), cudaMemcpyDeviceToHost));

    double variance = 0.0;
    for(long idx=0; idx<config.M/WINDOW_WIDTH; idx++) {
        if(config.DEBUG_LEVEL == 2)
            cout << "VARI[" << idx << "]: " << variances[idx] << endl;
        variance += variances[idx];
    }
    cout << endl;

    for(long idx=0; idx<config.M; idx++) {
        if(config.DEBUG_LEVEL == 2) {
            if(idx < 10 || idx > (config.M - 10))
                cout << "TRI[" << idx << "]: " << debug[idx] << endl;
        }
    }

    cudaFree(cudaDebug);
    cudaFree(cudaVariances);

    delete [] variances;
    delete [] debug;
    
    return sqrt(variance);
}
