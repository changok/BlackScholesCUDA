#include <iostream>
#include <cstddef>
#include <cassert>
#include <cmath>

#include <stdio.h>

#include <cuda.h>
#include <curand_kernel.h>

#include "common.h"
#include "black_scholes.cuh"
#include "timer.h"

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

__global__ void black_scholes_kernel(const double S, const double E,
            const double r, const double sigma, const double T,
            const long M, double* blockMeans, double* cudaTrials,
            curandState* randStates, const int mode, const double* fixedRands,
            double* debug, int debug_mode) {
    
    __shared__ double means[WINDOW_WIDTH];
   
    const long NUM_OF_TOT_THREAD = gridDim.x * blockDim.x;
    const long LOOP_SIZE = M / NUM_OF_TOT_THREAD;
    const unsigned int GID = (blockIdx.x * blockDim.x) * LOOP_SIZE + threadIdx.x * LOOP_SIZE;
    const unsigned int TID = threadIdx.x;

	
    curandState localState = randStates[(blockIdx.x * blockDim.x) + threadIdx.x];
    gaussrand_result_t gresult;

    means[TID] = 0.0;

    // Do the Black-Scholes iterations
    for(long trial = 0; trial < LOOP_SIZE; trial++) {
        double value = 0.0;
        if (trial%2 == 0) {
            if (mode == 1) {
                gresult.grand1 = 1.0;
                gresult.grand2 = 1.0;
            }
            // use pre-generated random number
            else if (mode == 2) {
                gresult.grand1 = fixedRands[GID + trial];
                gresult.grand2 = fixedRands[GID + trial+1];
            }
            // use gaussian random number (standard normal distributed)
            else {
                gresult = gaussrand (&localState);
            }

            if(debug_mode == 1) {
                debug[GID + trial] = gresult.grand1;
            }
            value = black_scholes_value (S, E, r, sigma, T, gresult.grand1);
        } else {
            if(debug_mode == 1) {
                debug[GID + trial] = gresult.grand2;
            }
            value = black_scholes_value (S, E, r, sigma, T, gresult.grand2);
        }

        // we need to keep origianl trial values for calculatng standard deviation
        // for current calculation, we use trials
        // Also, to prevent overflow caused by adding, divide the value by M in advance
        means[TID] += value/M;
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

    //*target = *target > 0.0 ? floor(*target + 0.00000000005) : ceil(*target - 0.0000000005);
}

__global__ void black_scholes_variance_kernel(const double mean,
            const long M, double* cudaTrials, double* cudaVariances) {
    
    __shared__ double variances[WINDOW_WIDTH];
    
    unsigned int gId = (blockIdx.x * WINDOW_WIDTH) + threadIdx.x; 
    unsigned int tId = threadIdx.x;
    
    variances[tId] = cudaTrials[gId];
    variances[tId] = variances[tId] - mean;
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

cit black_scholes(const double S, const double E, const double r,
                   const double sigma, const double T, const long M,
                   const int mode, double* cudafixedRands, int debug_mode) {
    cit interval;
    const long num_of_blocks = min((M/WINDOW_WIDTH), BLOCK_SIZE);
    long num_of_tot_threads = num_of_blocks * WINDOW_WIDTH;

    double* means = new double[num_of_blocks];
    double stddev = 0.0;
    double conf_width = 0.0;
	double t1, t2;

// part1_begin
	t1 = get_seconds();

    assert (M > 0);
    long size = M * sizeof(double);

    dim3 dimGrid(num_of_blocks);
    dim3 dimBlock(WINDOW_WIDTH);

	t2 = get_seconds();
	interval.t1 = t2-t1;	// init time
// part1_end

    // below pretend working with MOCKED thread number
// part2_begin
	t1 = 0; t1 = get_seconds();
    curandState* randStates;
    cudaMalloc((void **) &randStates, num_of_tot_threads * sizeof(curandState));
	t2 = 0; t2 = get_seconds();
    interval.t2 = t2- t1;	// setup_rnd_kernel time
// part2_end

// part5_begin
	t1 = 0; t1 = get_seconds();
    setup_rnd_kernel<<<dimGrid, dimBlock>>>(randStates, time(NULL));
	t2 = 0; t2 = get_seconds();
	interval.t5 = t2-t1;
// part5_end

// part3_begin
	t1 = 0; t1 = get_seconds();

    double* blockMeans;
    cudaMalloc((void**) &blockMeans, num_of_blocks * sizeof(double));
    
    double* cudaTrials;
    cudaMalloc((void**) &cudaTrials, size);

    double* hostDebug = new double[M];
    double* cudaDebug;
    if (debug_mode == 1) {
        cudaMalloc((void**) &cudaDebug, size);
    }

    black_scholes_kernel<<<dimGrid, dimBlock>>>(S, E, r, sigma, T, M, blockMeans, cudaTrials, randStates, mode, cudafixedRands, cudaDebug, debug_mode);

    cudaMemcpy(means, blockMeans, num_of_blocks * sizeof(double), cudaMemcpyDeviceToHost);

    if (debug_mode == 1) {
        cudaMemcpy(hostDebug, cudaDebug, size, cudaMemcpyDeviceToHost);
        for (int i = 0; i < M; i++) {
            printf("r%d: %lf, ", i, hostDebug[i]);
        }
        puts("\n");

        for (int i = 0; i < num_of_blocks; i++) {
            printf("m%d: %lf, ", i, means[i]);
        }
        puts("");

        double* t = new double[M];
        cudaMemcpy(t, cudaTrials, size, cudaMemcpyDeviceToHost);
        for (int i = 0; i < M; i++) {
            printf("t%d: %lf, ", i, t[i]);
        }
        puts("");

        delete [] t;
    }

	t2 =0; t2 = get_seconds();
	interval.t3 = t2-t1;	// black_scholes_kernel time
// part3_end

// part4_begin
    t1 = 0;
	t1 = get_seconds();
    double mean = 0.0;

    // combine results from each blocks
    for (long i = 0; i < num_of_blocks; i++) {
        mean += means[i];
    }
   
	printf("mean: %.20lf\n", mean);

    stddev = black_scholes_stddev (mean, M, cudaTrials);
    cout << "StdDev: " << stddev << endl;


	t2 = 0;
	t2 = get_seconds();
	interval.t4 = t2-t1;
// part4_end

    /* confidence interval */
    conf_width = 1.96 * stddev / sqrt ((double) M);
    interval.min = mean - conf_width;
    interval.max = mean + conf_width;

    /* clean up */
	cudaFree(cudaDebug);
    cudaFree(cudaTrials);
    cudaFree(blockMeans);
    cudaFree(randStates);

    delete [] hostDebug;
    delete [] means;
    
    return interval;
}

double black_scholes_stddev (const double mean, const long M, double* cudaTrials) {
    double* variances = new double[M/WINDOW_WIDTH];
    double* cudaVariances;
    cudaMalloc((void**) &cudaVariances, M/WINDOW_WIDTH * sizeof(double));
    
//    dim3 dimGrid(M/WINDOW_WIDTH);
    dim3 dimGrid(M/WINDOW_WIDTH);
    dim3 dimBlock(WINDOW_WIDTH);
    black_scholes_variance_kernel<<<dimGrid, dimBlock>>>(mean, M, cudaTrials, cudaVariances);

    cudaMemcpy(variances, cudaVariances, M/WINDOW_WIDTH * sizeof(double), cudaMemcpyDeviceToHost);
    
    double variance = 0.0;
    for(long idx=0; idx<M/WINDOW_WIDTH; idx++) {
        cout << "std: " << variances[idx] << ", ";
        variance += variances[idx];
    }
    cout << endl;

    cudaFree(cudaVariances);
    
    return sqrt(variance);
}
