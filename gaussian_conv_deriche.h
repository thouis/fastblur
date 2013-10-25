/**
 * \file gaussian_conv_deriche.h
 * \brief Deriche's approximation of Gaussian convolution
 * \author Pascal Getreuer <getreuer@gmail.com>
 * 
 * Copyright (c) 2012-2013, Pascal Getreuer
 * All rights reserved.
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as 
 * published by the Free Software Foundation, either version 3 of the 
 * License, or (at your option) any later version, or the terms of the 
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along with this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

/** 
 * \defgroup deriche_gaussian Deriche Gaussian convolution
 * 
 * Deriche uses a sum of causal and an anticausal recursive filters to 
 * approximate the Gaussian.  The filters were designed with orders K = 2, 3,
 * and 4, allowing to trade between better accuracy vs. faster computation
 * time.
 * 
 * The process to use these functions is the following:
 *    -# deriche_precomp() to precompute coefficients for the convolution
 *    -# deriche_gaussian_conv() or deriche_gaussian_conv_image() to perform
 *       the convolution itself (may be called multiple times if desired)
 *
 * \par Example
\code
    deriche_coeffs c;
    num *buffer;
    
    deriche_precomp(&c, sigma, K, tol);
    buffer = (num *)malloc(sizeof(num) * 2 * N);
    deriche_gaussian_conv(c, dest, buffer, src, N, stride);
    free(buffer);
\endcode
 * 
 * \{
 */

#ifndef _GAUSSIAN_CONV_DERICHE_H_
#define _GAUSSIAN_CONV_DERICHE_H_

#include "num.h"

/** \brief Minimum Deriche filter order */
#define DERICHE_MIN_K       2
/** \brief Maximum Deriche filter order */
#define DERICHE_MAX_K       4
/** \brief Test whether a given K value is a valid Deriche filter order */
#define DERICHE_VALID_K(K)  (DERICHE_MIN_K <= (K) && (K) <= DERICHE_MAX_K)

/** \brief Coefficients for Deriche Gaussian approximation */
typedef struct deriche_coeffs
{
    num a[DERICHE_MAX_K + 1];             /**< Denominator coeffs          */
    num b_causal[DERICHE_MAX_K];          /**< Causal numerator            */
    num b_anticausal[DERICHE_MAX_K + 1];  /**< Anticausal numerator        */
    num sum_causal;                       /**< Causal filter sum           */
    num sum_anticausal;                   /**< Anticausal filter sum       */
    num sigma;                            /**< Gaussian standard deviation */
    int K;                                /**< Filter order = 2, 3, or 4   */
    num tol;                              /**< Boundary accuracy           */
    long max_iter;
} deriche_coeffs;

void deriche_precomp(deriche_coeffs *c, double sigma, int K, num tol);
void deriche_gaussian_conv(deriche_coeffs c,
    num *dest, num *buffer, const num *src, long N, long stride);
void deriche_gaussian_conv_image(deriche_coeffs c,
    num *dest, num *buffer, const num *src, 
    int width, int height, int num_channels);

/** \} */
#endif /* _GAUSSIAN_CONV_DERICHE_H_ */
