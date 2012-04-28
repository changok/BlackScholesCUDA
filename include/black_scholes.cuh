#ifndef _black_scholes_h
#define _black_scholes_h

#include "common.h"
#include "bsconfig.h"

typedef struct __gaussrand_result_t {
  double grand1, grand2;
} gaussrand_result_t;

Result black_scholes(double* cudafixedRands, BSConfig config);

double black_scholes_stddev (const double mean, const long M, double* trials);

#endif /* _black_scholes_h */
