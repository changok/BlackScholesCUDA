#include <iostream>
#include <limits.h>
#include <cstdlib>
#include <stdio.h>

#include <cuda.h>
#include <curand_kernel.h>
#include "cutil.h"
#include "cutil_inline_runtime.h"

#include "common.h"
#include "parser.h"
#include "timer.h"
#include "black_scholes.cuh"
#include "bsconfig.h"

using namespace std;

/**
 * Usage: ./blackScholes <filename> [Trials(M)] [Random Mode] [Debug Flag]
 *
 * <filename> (don't include the angle brackets) is the name of 
 * a data file in the current directory containing the parameters
 * for the Black-Scholes simulation.  It has exactly six lines 
 * with no white space.  Put each parameter one to a line, with
 * an endline after it.  Here are the parameters:
 *
 * S
 * E
 * r
 * sigma
 * T
 *
 * [Random Type] (don't include the brackets) is used for specify random
 * number generator type. It can be omitted.
 * 0, or nothing: return Gaussian Number (Standard Normal Distributed Random Number)
 * 1: Test purpose generator. Always returns 1
 * 2: Test purpose generator. It returns one element from pre-generated sequence
 *    To run program under this mode, un-comment below line in Makefile
 *
 *    #FLAG = -D__GOGO_DEBUG__
 *
 *    Then, compile again.
 *
 * [Debug Flag] print out with debug mode
 * 0: Default. No show any additional information
 * 1: Show consumed time to process each major parts
 * 2: Verbose mode. To run program under this mode, un-comment below line in Makefile
 *
 *    #FLAG = -D__GOGO_DEBUG__
 *
 *    Then, compile again.
 */
int main(int argc, char* argv[]) {
    BSConfig config;
    char* filename = NULL;
    double t1, t2;

    if (argc < 2) {
        cerr << "Usage: ./blackScholes <filename> [M:Number of Trials] [Mode]" << endl << endl;
        exit(EXIT_FAILURE);
    }
    filename = argv[1];
    parse_parameters(&config, filename);

    if (argv[2] != NULL) {
	  cout << "Number of Trials[M] : " << argv[2] << endl;
	  config.M = to_long(argv[2]);
    }

    if (argv[3] != NULL) {
        config.RND_MODE = to_int(argv[3]);
#ifdef __GOGO_DEBUG__
		if (config.RND_MODE > 2) {
			cerr << "Available Random Mode: [0], [1], [2]" << endl << endl;
			exit(EXIT_FAILURE);
		}
#else
		if (config.RND_MODE > 1) {
            cerr << "Available Random Mode: [0], [1]" << endl << endl;
            exit(EXIT_FAILURE);
        }
#endif
    }

	if (argv[4] != NULL) {
		config.DEBUG_LEVEL = to_int(argv[4]);
#ifdef __GOGO_DEBUG__
		if (config.DEBUG_LEVEL > 2) {
			cout << "Only three debug mode are possible[0][1][2], default is 0. Thus set as 0" << endl;
			config.DEBUG_LEVEL = 0;
		} else if (config.DEBUG_LEVEL == 2) {
			cout << "Verbose Debug Mode ON" << endl;
		} else if (config.DEBUG_LEVEL == 1) {
			cout << "Debug Mode ON" << endl;
		}
#else
		if (config.DEBUG_LEVEL > 1) {
            cout << "Only two debug mode are possible[0][1], default is 0. Thus set as 0" << endl;
            config.DEBUG_LEVEL = 0;
		} else if (config.DEBUG_LEVEL == 1) {
            cout << "Debug Mode ON" << endl;
        }
#endif
	}

    /*
     * Make sure init_timer() is only called by one thread,
     * before all the other threads run!
     */
    init_timer();

    /*
     * Run the benchmark and time it.
     */

	if (config.M < WINDOW_WIDTH) {
		cout << "M(trials) is smaller than minimum requirement(" << WINDOW_WIDTH << "). So, automatically set as the minimum." << endl;
		config.M = WINDOW_WIDTH;
	}

#ifdef __GOGO_DEBUG__
	// pre-generated fixed numbers as random for correctness test : mode[2]
	double* fixedRands = new double[config.M];
	double* cudafixedRands;
    for (int i = 0; i < config.M; i++) {
      fixedRands[i] = i/(double)config.M;
    }

    cutilSafeCall(cudaMalloc((void**) &cudafixedRands, config.M*sizeof(double)));
    cutilSafeCall(cudaMemcpy(cudafixedRands, fixedRands, config.M*sizeof(double), cudaMemcpyHostToDevice));
#endif

    t1 = get_seconds();
    /*
     * In the parallel case, you may want to set prng_stream_spawn_time to
     * the max of all the prng_stream_spawn_times, or just take a representative
     * sample...
     */

#ifdef __GOGO_DEBUG__
    Result result = black_scholes(cudafixedRands, config);
#else
    Result result = black_scholes(config);
#endif
    //Result result = black_scholes(S, E, r, sigma, T, M, mode, cudafixedRands, debug_mode, config);

    t2 = get_seconds();

    /*
     * A fun fact about C string literals (i.e., strings enclosed in
     * double quotes) is that the C preprocessor automatically
     * concatenates them if they are separated only by whitespace.
     */

    cout << "Black-Scholes in GPU benchmark :" << endl
            << "------------------------" << endl
            << "S        " << config.S << endl
            << "E        " << config.E << endl
            << "r        " << config.r << endl
            << "sigma    " << config.sigma << endl
			<< "T        " << config.T << endl
			<< "M        " << config.M << endl;

    cout.precision(7);
    cout << endl;
    cout << "Mean: " << result.mean << endl;
    cout << "Standard Deviation: " << result.stddev << endl;
    cout << "Confidence interval: [" << result.min << ", " << result.max << "]" << endl;

    cout.precision(10);
    cout << "Overall Simulation Time     : " << t2 - t1 << " seconds" << endl;

    if(config.DEBUG_LEVEL > 0) {
        cout << "Random streams Initializing : " << result.init_seeds_setup_time << " seconds" << endl;
        cout << "Black Scholes Kernel        : " << result.black_sholes_kernel_time << " seconds" << endl;
        cout << "Standard Deviation Kernel   : " << result.calc_stddev_time << " seconds" << endl;
    }

#ifdef __GOGO_DEBUG__
	cudaFree (cudafixedRands);
	if(fixedRands != NULL) delete [] fixedRands;
#endif

    return 0;
}

