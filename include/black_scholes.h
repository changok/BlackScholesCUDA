#ifndef _black_scholes_h
#define _black_scholes_h

#include "common.h"

/**
 * Run the Black-Scholes MC simulation using the parameters S, E, r,
 * sigma, and T, with M total trials.  
 *
 */
cit* black_scholes(const double S, const double E, const double r,
                  const double sigma, const double T, const long M);

#endif /* _black_scholes_h */
