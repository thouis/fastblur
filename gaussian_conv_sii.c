/**
 * \file gaussian_conv_sii.c
 * \brief Gaussian convolution using stacked integral images
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

#include "gaussian_conv_sii.h"
#include <assert.h>
#include <stdio.h>
#include "filter_util.h"

#ifndef M_PI
/** \brief The constant pi */
#define M_PI        3.14159265358979323846264338327950288
#endif

/**
 * \brief Precompute filter coefficients for SII Gaussian convolution
 * \param c         sii_coeffs pointer to hold precomputed coefficients
 * \param sigma     Gaussian standard deviation
 * \param K         number of boxes = 3, 4, or 5
 * \return 1 on success, 0 on failure
 */
void sii_precomp(sii_coeffs *c, double sigma, int K)
{
    const double sigma0 = 100.0 / M_PI;
    static const short radii0[SII_MAX_K - SII_MIN_K + 1][SII_MAX_K] = 
        {{76, 46, 23, 0, 0},
         {82, 56, 37, 19, 0}, 
         {85, 61, 44, 30, 16}};
    static const float weights0[SII_MAX_K - SII_MIN_K + 1][SII_MAX_K] =
        {{0.1618f, 0.5502f, 0.9495f, 0, 0},
         {0.0976f, 0.3376f, 0.6700f, 0.9649f, 0},
         {0.0739f, 0.2534f, 0.5031f, 0.7596f, 0.9738f}};
    const int i = K - SII_MIN_K;
    double sum;
    int k;
    
    assert(c && sigma > 0 && SII_VALID_K(K));
    c->K = K;
    
    for(k = 0, sum = 0; k < K; k++)
    {
        c->radii[k] = (long)(radii0[i][k] * (sigma / sigma0) + 0.5);
        sum += weights0[i][k] * (2 * c->radii[k] + 1);
    }
    
    for(k = 0; k < K; k++)
        c->weights[k] = (num)(weights0[i][k] / sum);
    
    return;
}

/**
 * \brief Determines the buffer size needed for SII Gaussian convolution
 * \param c     sii_coeffs created by sii_precomp()
 * \param N     number of samples
 * \return required buffer size in units of num samples
 */
long sii_buffer_size(sii_coeffs c, long N)
{
    long pad = c.radii[0] + 1;
    return N + 2 * pad;
}

/**
 * \brief Gaussian convolution SII approximation
 * \param c         sii_coeffs created by sii_precomp()
 * \param dest      output convolved data
 * \param buffer    array with space for sii_buffer_size() samples
 * \param src       data to be convolved
 * \param N         number of samples
 * \param stride    stride between successive samples
 * 
 * The buffer should have space for at least sii_buffer_size(c,N) samples.
 */
void sii_gaussian_conv(sii_coeffs c, num *dest, num *buffer,
    const num *src, long N, long stride)
{
    num accum;
    long pad, n;
    int k;
    
    assert(dest && buffer && src && dest != buffer && src != buffer
        && N > 0 && stride != 0);
    
    pad = c.radii[0] + 1;
    buffer += pad;
    
    /* Compute cumulative sum */
    for(n = -pad, accum = 0; n < N + pad; n++)
        buffer[n] = accum += src[stride * extension(N, n)];
    
    /*for(n = -pad, accum = 0; n < 0; n++)
        sum[n] = accum += src[stride * extension(N, n)];
    
    for(; n < N; n++)
        sum[n] = accum += src[stride * n];
    
    for(; n < N + pad; n++)
        sum[n] = accum += src[stride * extension(N, n)];*/
    
    /* Compute summed box filters */
    for(n = 0; n < N; n++, dest += stride)
    {
        accum = c.weights[0] * (buffer[n + c.radii[0]] 
            - buffer[n - c.radii[0] - 1]);
        
        for(k = 1; k < c.K; k++)
            accum += c.weights[k] * (buffer[n + c.radii[k]] 
                - buffer[n - c.radii[k] - 1]);
        
        *dest = accum;
    }
    
    return;
}

/**
 * \brief 2D Gaussian convolution SII approximation
 * \param c             sii_coeffs created by sii_precomp()
 * \param dest          output convolved data
 * \param buffer        array with space for sii_buffer_size() samples
 * \param src           image to be convolved
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 * 
 * The buffer should have space for at least 
 * sii_buffer_size(c,max(width,height)) samples.
 */
void sii_gaussian_conv_image(sii_coeffs c, num *dest, num *buffer, 
    const num *src, int width, int height, int num_channels)
{
    long num_pixels = ((long)width) * ((long)height);
    int x, y, channel;
    
    assert(dest && buffer && src && num_pixels > 0);
    
    for(channel = 0; channel < num_channels; channel++,
        dest += num_pixels, src += num_pixels)
    {
        num *dest_y = dest;
        const num *src_y = src;
            
        for(y = 0; y < height; y++, dest_y += width, src_y += width)
            sii_gaussian_conv(c,
                dest_y, buffer, src_y, width, 1);
        
        for(x = 0; x < width; x++)
            sii_gaussian_conv(c, 
                dest + x, buffer, dest + x, height, width);
    }
    
    return;
}
