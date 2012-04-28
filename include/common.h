/*
 * common.h
 *
 *  Created on: Apr 22, 2012
 *      Author: jinyoung
 */

#ifndef COMMON_H_
#define COMMON_H_

static int rnd_mode = 0;

/* type for debug */
typedef struct __debug_t {
  int nthreads;
  int nblocks;
  double* rands;
} debug_t;

/**
 * Confidence interval [min,max].
 *
 * The typedef lets you refer to this type without saying "struct"
 * in front of it, but it means that you have to include this header
 * file if you ever use this datatype.
 */
typedef struct __confidence_interval_t {
  double min, max;
  double t1, t2, t3, t4, t5;
} confidence_interval_t;

typedef confidence_interval_t cit;

typedef struct __run_result__ {
    double min, max; //[OUT] Arithmetic mean of trials[0 .. M-1].
    double mean;
    double stddev;
    double mem_for_rnd_init_time;
    double init_seeds_setup_time;
    double calc_stddev_time;
    double black_sholes_kernel_time;
} Result;


const long BLOCK_SIZE = 128;
const int WINDOW_WIDTH = 256;

#endif /* COMMON_H_ */
