/*
 * common.h
 *
 *  Created on: Apr 22, 2012
 *      Author: jinyoung
 */

#ifndef COMMON_H_
#define COMMON_H_

static int rnd_mode = 0;

/**
 * Confidence interval [min,max].
 *
 * The typedef lets you refer to this type without saying "struct"
 * in front of it, but it means that you have to include this header
 * file if you ever use this datatype.
 */
typedef struct __confidence_interval_t {
  double min, max;
} confidence_interval_t;

const typedef confidence_interval_t cit;

#endif /* COMMON_H_ */
