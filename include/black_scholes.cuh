#ifndef _black_scholes_h
#define _black_scholes_h

#include "common.h"

typedef struct __gaussrand_result_t {
  double grand1, grand2;
} gaussrand_result_t;

typedef struct __black_scholes_args_t {
  double S, E, r, sigma, T; //[IN] Various parameters of the Black-Scholes MC method.
  long M; //[IN] Number of Black-Scholes MC iterations.
  double mean; //[OUT] Arithmetic mean of trials[0 .. M-1].
  double variance; //[OUT] variance of trials[0 .. M-1].
} black_scholes_args_t;

typedef black_scholes_args_t bca_t;

/**
 * Run the Black-Scholes MC simulation using the parameters S, E, r,
 * sigma, and T, with M total trials.  
 *
 */
cit black_scholes(const double S, const double E, const double r,
                  const double sigma, const double T, const long M);

void deinit_black_scholes_args (bca_t* args);
double black_scholes_stddev (const double mean, const long M, double* trials);
//static double black_scholes_stddev (const double mean, const long M, const double* trials);

#endif /* _black_scholes_h */
