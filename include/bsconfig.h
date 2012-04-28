/*
 * bsconfig.h
 *
 *  Created on: Apr 27, 2012
 *      Author: jinyoung
 */

#ifndef BSCONFIG_H_
#define BSCONFIG_H_

#include <iostream>

class BSConfig {
public:
    double S, E, r, sigma, T; //[IN] Various parameters of the Black-Scholes MC method.
    long M; //[IN] Number of Black-Scholes MC iterations.
    int DEBUG_LEVEL;
    int RND_MODE;

    long totalNumOfBlocks();
    long totalNumOfThread();

    BSConfig();
};


#endif /* BSCONFIG_H_ */
