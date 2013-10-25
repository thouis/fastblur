/**
 * \file gaussian_conv_box.c
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

#include "gaussian_conv_box.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "filter_util.h"

/**
 * \brief Perform one pass of box filtering
 * \param dest          destination array
 * \param dest_stride   stride between successive samples of dest
 * \param src_stride    src array (must be distinct from dest)
 * \param src_stride    stride between successive samples of src
 * \param N             number of samples
 * \param r             radius of the box filter
 * 
 * Performs one pass of box filtering with radius r (diameter 2r+1).
 */
static void box_filter(num *dest, long dest_stride, 
    const num *src, long src_stride, long N, long r)
{
    long n;
    num accum;
    
    assert(dest && src && dest != src && N > 0 && r >= 0);
    
    for(n = -r, accum = 0; n <= r; n++)
        accum += src[src_stride * extension(N, n)];
    
    dest[0] = accum;
    
    for(n = 1; n < N; n++)
    {
        accum += src[src_stride * extension(N, n + r)]
            - src[src_stride * extension(N, n - r - 1)];
        dest[dest_stride * n] = accum;
    }
    
    return;
}

/** 
 * \brief Box filtering approximation of Gaussian convolution
 * \param dest_data     destination array
 * \param buffer_data   array with space for at least N samples
 * \param src           src array
 * \param N             number of samples
 * \param stride        stride between successive samples
 * \param sigma         Gaussian standard deviation in pixels
 * \param K             number of box filter passes
 * 
 * This routine performs iterated box filtering approximation of Gaussian
 * convolution.  Well's approximation 
 * \f$ \sigma^2 = \tfrac{1}{12} K \bigl((2r+1)^2 - 1\bigr) \f$
 * is used to select the box filter radius.
 */ 
void box_gaussian_conv(num *dest_data, num *buffer_data, const num *src,
    long N, long stride, num sigma, int K)
{  
    struct
    {
        num *data;
        long stride;
    } dest, buffer, cur, next;
    num scale;
    long r;
    int step;
    
    assert(dest_data && buffer_data && src && dest_data != buffer_data
        && N > 0 && sigma > 0 && K > 0);
    
    r = (long)(0.5 * sqrt((12.0 * sigma * sigma) / K + 1.0));
    scale = (num)(1.0 / pow(2*r + 1, K));
    dest.data = dest_data;
    dest.stride = stride;
    buffer.data = buffer_data;
    buffer.stride = (buffer_data == src) ? stride : 1;
    
    next = (buffer_data == src || (dest_data != src && K % 2 == 1)) 
        ? dest : buffer;
    box_filter(next.data, next.stride, src, stride, N, r);
    
    for(step = 1; step < K; step++)
    {
        cur = next;
        next = (cur.data == buffer_data) ? dest : buffer;
        box_filter(next.data, next.stride, cur.data, cur.stride, N, r);
    }
    
    if(next.data != dest_data)
    {
        long n, i;
        
        for(n = i = 0; n < N; n++, i += stride)
            dest_data[i] = buffer_data[n] * scale;
    }
    else
    {
        long i, i_end = stride * N;
        
        for(i = 0; i < i_end; i += stride)
            dest_data[i] *= scale;
    }
    
    return;
}

/** 
 * \brief Box filtering approximation of 2D Gaussian convolution
 * \param dest          destination image
 * \param buffer        array with at least max(width,height) samples
 * \param src           src image
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 * \param sigma         Gaussian standard deviation in pixels
 * \param K             number of box filter passes
 */
void box_gaussian_conv_image(num *dest, num *buffer, const num *src,
    int width, int height, int num_channels, num sigma, int K)
{
    const long num_pixels = ((long)width) * ((long)height);
    int x, y, channel;
    
    assert(dest && buffer && src && dest != buffer
        && num_pixels > 0 && sigma > 0 && K > 0);
    
    for(channel = 0; channel < num_channels; channel++,
        dest += num_pixels, src += num_pixels)
    {
        num *dest_y = dest;
        const num *src_y = src;
            
        for(y = 0; y < height; y++, dest_y += width, src_y += width)
            box_gaussian_conv(dest_y, buffer, src_y, 
                width, 1, sigma, K);
        
        for(x = 0; x < width; x++)
            box_gaussian_conv(dest + x, buffer, dest + x, 
                height, width, sigma, K);
    }
    
    return;
}
