/**
 * \file gaussian_conv_ebox.c
 * \brief Gaussian convolution with extended box filters
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
 * 
 * box filter approximation references:
 * 
 * Wells, W.M.: Efficient synthesis of Gaussian filters by cascaded uniform filters.
 * IEEE Transactions on Pattern Analysis and Machine Intelligence 8(2), 234--239
 * (Mar 1986)
 * 
 * R. Rau,  J.H. McClellan, "Efficient approximation of Gaussian filters,"
 * IEEE Transactions on Signal Processing, Volume 45, Issue 2, February 1997,
 * Page 468-471.
 * 
 * Peter Kovesi: Fast Almost-Gaussian Filtering. DICTA 2010: 121-125
 * 
 * Pascal Gwosdek, Sven Grewenig, Andr\'es Bruhn, and Joachim Weickert, 
 * "Theoretical Foundations of Gaussian Convolution by Extended Box Filtering,"
 * SSVM 2011, 447-458.
 * 
 * "Efficient and Accurate Gaussian Image Filtering Using Running Sums"
 * Elhanan Elboher, Michael Werman, arXiv:1107.4958
 * 
 */

#include "gaussian_conv_ebox.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include "filter_util.h"

/**
 * \brief Precompute coefficients for extended box filtering
 * \param c         ebox_coeffs pointer to hold precomputed coefficients
 * \param sigma     Gaussian standard deviation
 * \param K         number of filtering passes
 * 
 * This routine precomputes the coefficients for extended box filtering
 * Gaussian convolution approximation.
 */
void ebox_precomp(ebox_coeffs *c, double sigma, int K)
{
    double alpha;
    
    assert(c && sigma > 0 && K > 0);
    
    c->r = (long)(0.5 * sqrt((12.0 * sigma * sigma) / K + 1.0) - 0.5);
    alpha = (2 * c->r + 1) * (c->r * (c->r + 1) - 3.0 * sigma * sigma / K)
        / (6.0 * (sigma * sigma / K - (c->r + 1) * (c->r + 1)));
    c->c_1 = alpha / (2.0 * (alpha + c->r) + 1);
    c->c_2 = (1.0 - alpha) / (2.0 * (alpha + c->r) + 1);
    c->K = K;
    return;
}

/**
 * \brief Perform one pass of extended box filtering
 * \param dest          destination array
 * \param dest_stride   stride between successive samples of dest
 * \param src_stride    src array (must be distinct from dest)
 * \param src_stride    stride between successive samples of src
 * \param N             number of samples
 * \param r             radius of the inner box
 * \param c_1           weight of the outer box
 * \param c_2           weight of the inner box
 */
static void ebox_filter(num *dest, long dest_stride,
    const num *src, long src_stride, long N, long r, num c_1, num c_2)
{
    long n;
    num accum = 0;
    
    assert(dest && src && dest != src && N > 0 && r >= 0);
    
    for(n = -r; n <= r; n++)
        accum += src[src_stride * extension(N, n)];
    
    dest[0] = accum = c_1 * (src[src_stride * extension(N, r + 1)]
        + src[src_stride * extension(N, -r - 1)])
        + (c_1 + c_2) * accum;
    
    for(n = 1; n < N; n++)
    {
        accum += c_1 * (src[src_stride * extension(N, n + r + 1)]
            - src[src_stride * extension(N, n - r - 2)])
            + c_2 * (src[src_stride * extension(N, n + r)]
            - src[src_stride * extension(N, n - r - 1)]);
        dest[dest_stride * n] = accum;
    }
    
    return;
}

/** 
 * \brief Extended box filtering approximation of Gaussian convolution
 * \param c             ebox_coeffs created by ebox_precomp()
 * \param dest_data     destination array
 * \param buffer_data   array with space for at least N samples
 * \param src           src array
 * \param N             number of samples
 * \param stride        stride between successive samples
 */ 
void ebox_gaussian_conv(ebox_coeffs c, num *dest_data, num *buffer_data,
    const num *src, long N, long stride)
{  
    struct
    {
        num *data;
        long stride;
    } dest, buffer, cur, next;
    int step;
    
    assert(dest_data && buffer_data && src
        && dest_data != buffer_data && N > 0);
    
    dest.data = dest_data;
    dest.stride = stride;
    buffer.data = buffer_data;
    buffer.stride = (buffer_data == src) ? stride : 1;
    
    next = (buffer_data == src || (dest_data != src && c.K % 2 == 1)) 
        ? dest : buffer;
    ebox_filter(next.data, next.stride, src, stride,
        N, c.r, c.c_1, c.c_2);
    
    for(step = 1; step < c.K; step++)
    {
        cur = next;
        next = (cur.data == buffer_data) ? dest : buffer;
        ebox_filter(next.data, next.stride, cur.data, cur.stride,
            N, c.r, c.c_1, c.c_2);
    }
    
    if(next.data != dest_data)
    {
        if(stride == 1)
            memcpy(dest_data, buffer_data, sizeof(num) * N);
        else
        {
            long n, i;
            
            for(n = i = 0; n < N; n++, i += stride)
                dest_data[i] = buffer_data[n];
        }
    }
    
    return;
}

/** 
 * \brief Extended box filtering approximation of 2D Gaussian convolution
 * \param c             ebox_coeffs created by ebox_precomp()
 * \param dest          destination image
 * \param buffer        array with at least max(width,height) samples
 * \param src           src image
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 */
void ebox_gaussian_conv_image(ebox_coeffs c, num *dest, num *buffer,
    const num *src, int width, int height, int num_channels)
{
    const long num_pixels = ((long)width) * ((long)height);
    int x, y, channel;
    
    assert(dest && buffer && src && dest != buffer && num_pixels > 0);
    
    for(channel = 0; channel < num_channels; channel++,
        dest += num_pixels, src += num_pixels)
    {
        num *dest_y = dest;
        const num *src_y = src;
            
        for(y = 0; y < height; y++, dest_y += width, src_y += width)
            ebox_gaussian_conv(c, dest_y, buffer, src_y, width, 1);
        
        for(x = 0; x < width; x++)
            ebox_gaussian_conv(c, dest + x, buffer, dest + x, height, width);
    }
    
    return;
}
