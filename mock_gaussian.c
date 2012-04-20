/*
 * mock_gaussian.c
 *
 *  Created on: Apr 20, 2012
 *      Author: jinyoung
 */

#include "mock_gaussian.h"

double gaussrand_only1 (const double_generator_one_input_t f,
        void* f_state,
        gaussrand_state_t* gaussrand_state) {
    return (double)1;
}

double gaussrand_pre_generated (const double_generator_one_input_t f,
        void* f_state,
        gaussrand_state_t* gaussrand_state) {
    return 0;
}
