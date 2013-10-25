/**
 * \file gaussian_conv_box.h
 * \brief Box filtering approximation of Gaussian convolution
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
 * \defgroup box_gaussian Box filter Gaussian convolution
 * 
 * This file implments iterated box filter approximation of Gaussian
 * convolution.
 *
 * \par Example
\code
    num *buffer;
    
    buffer = (num *)malloc(sizeof(num) * N);
    box_gaussian_conv(dest, buffer, src, N, stride, sigma, K);
    free(buffer);
\endcode
 * 
 * \{
 */
#ifndef _GAUSSIAN_CONV_BOX_H_
#define _GAUSSIAN_CONV_BOX_H_

#include "num.h"

void box_gaussian_conv(num *dest, num *buffer, const num *src,
    long N, long stride, num sigma, int K);
void box_gaussian_conv_image(num *dest, num *buffer, const num *src,
    int width, int height, int num_channels, num sigma, int K);

/** \} */
#endif /* _GAUSSIAN_CONV_BOX_H_ */
