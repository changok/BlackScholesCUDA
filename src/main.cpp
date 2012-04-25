#include <iostream>
#include <limits.h>
#include <cstdlib>
#include <stdio.h>

#include "common.h"
#include "parser.h"
#include "timer.h"
#include "black_scholes.cuh"

using namespace std;

/**
 * Usage: ./blackScholes <filename> [Random Mode]
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
 * M
 *
 * [Random Type] (don't include the brackets) is used for specify random
 * number generator type. It can be omitted.
 * 0: Normal Generator
 * 1: Test purpose generator. Always returns 1
 * 2: Test purpose generator. It returns one element from pre-generated sequence
 */
int main(int argc, char* argv[]) {
    double S, E, r, sigma, T;
    long M = 0;
    char* filename = NULL;
    double t1, t2;
//  double prng_stream_spawn_time;

    if (argc < 2) {
        //cerr << "Usage: ./blackScholes <filename> [Random Mode]" << endl << endl;
        cerr << "Usage: ./blackScholes <filename> [Number of Trials]" << endl << endl;
        exit(EXIT_FAILURE);
    }
    filename = argv[1];
    parse_parameters(&S, &E, &r, &sigma, &T, &M, filename);

    if (argv[2] != NULL) {
cout << "Number of Trials[M] : " << argv[2] << endl;
      M = to_long(argv[2]);
    }

    if (argc == 4) {
        rnd_mode = to_int(argv[3]);
    }

    /*
     * Make sure init_timer() is only called by one thread,
     * before all the other threads run!
     */
    init_timer();

    /*
     * Run the benchmark and time it.
     */
    t1 = get_seconds();
    /*
     * In the parallel case, you may want to set prng_stream_spawn_time to
     * the max of all the prng_stream_spawn_times, or just take a representative
     * sample...
     */
    cit interval = black_scholes(S, E, r, sigma, T, M);
    t2 = get_seconds();

    /*
     * A fun fact about C string literals (i.e., strings enclosed in
     * double quotes) is that the C preprocessor automatically
     * concatenates them if they are separated only by whitespace.
     */

    cout << "Black-Scholes in GPU benchmark :" << endl
            << "------------------------" << endl
            << "S        " << S << endl
            << "E        " << E << endl
            << "r        " << r << endl
            << "sigma    " << sigma << endl 
			<< "T        " << T << endl 
			<< "M        " << M << endl;

    cout << "Confidence interval: " << interval.min << ", " << interval.max << endl;
    cout << "Total simulation time: " << t2 - t1 << " seconds" << endl;
    cout << "part1 time: "; printf("%.20lf seconds\n", interval.t1);
    cout << "part2 time: "; printf("%.20lf seconds\n", interval.t2);
    cout << "part3 time: "; printf("%.20lf seconds\n", interval.t3);
    cout << "part4 time: "; printf("%.20lf seconds\n", interval.t4);
    cout << "part5 time: "; printf("%.20lf seconds\n", interval.t5);
    //cout << "setup_rnd_kernel time: " << interval.t2 << " seconds" << endl;
    //cout << "black_scholes_kernel time: " << interval.t2 << " seconds" << endl;
    //cout << "black_scholes_kernel time: " << interval.t2 << " seconds" << endl;
    //cout << "black_scholes_stddev time: " << interval.t3 << " seconds" << endl;
//  cout << "PRNG stream spawn time: " << prng_stream_spawn_time << " seconds" << endl;
    return 0;
}

