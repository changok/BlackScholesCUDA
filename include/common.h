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

const long BLOCK_SIZE = 128;
const int WINDOW_WIDTH = 256;

#endif /* COMMON_H_ */
