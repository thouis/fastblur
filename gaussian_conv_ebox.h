/**
 * \file gaussian_conv_ebox.h
 * \brief Extended box filter Gaussian approximation
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
 * \defgroup ebox_gaussian Extended box filter Gaussian convolution
 * 
 * This file implments the extended box filter approximation of Gaussian
 * convolution proposed by Gwosdek, Grewenig, Bruhn, and Weickert.
 * 
 * The process to use these functions is the following:
 *    -# ebox_precomp() to precompute coefficients for the convolution
 *    -# ebox_gaussian_conv() or ebox_gaussian_conv_image() to perform
 *       the convolution itself (may be called multiple times if desired)
 *
 * \par Example
\code
    ebox_coeffs c;
    num *buffer;
    
    ebox_precomp(&c, sigma, K);
    buffer = (num *)malloc(sizeof(num) * N);
    ebox_gaussian_conv(c, dest, buffer, src, N, stride);
    free(buffer);
\endcode
 * 
 * \{
 */

#ifndef _GAUSSIAN_CONV_EBOX_H_
#define _GAUSSIAN_CONV_EBOX_H_

#include "num.h"

/** \brief Coefficients for extended box filter Gaussian approximation */
typedef struct ebox_coeffs_
{
    num c_1;        /**< Outer box weight           */
    num c_2;        /**< Inner box weight           */
    long r;         /**< Inner box radius           */
    int K;          /**< Number of filtering passes */
} ebox_coeffs;

void ebox_precomp(ebox_coeffs *c, double sigma, int K);
void ebox_gaussian_conv(ebox_coeffs c, num *dest, num *buffer,
    const num *src, long N, long stride);
void ebox_gaussian_conv_image(ebox_coeffs c, num *dest, num *buffer,
    const num *src, int width, int height, int num_channels);

/** \} */
#endif /* _GAUSSIAN_CONV_EBOX_H_ */
