/**
 * \file gaussian_conv_am.c
 * \brief Alvarez-Mazorra approximate Gaussian convolution
 * \author Pascal Getreuer <getreuer@gmail.com>
 * 
 * Copyright (c) 2011-2013, Pascal Getreuer
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

#include <assert.h>
#include <math.h>
#include "filter_util.h"
#include "gaussian_conv_am.h"

/**
 * \brief Handling of the left boundary for Alvarez-Mazorra
 * \param data the signal data
 * \param N number of elements
 * \param stride the stride between successive samples
 * \param nu filter parameter nu
 * \param num_terms number of terms to use to approximate infinite sum
 * \return the sum approximating the first filtered sample value
 * 
 * This routine approximates the infinite sum
 * \f$ u_0 = \sum_{j=0}^\infty \nu^j \Tilde{x}_j \f$
 * by adding the first \c num_terms terms.
 */
static num am_left_boundary(const num *data, long N, long stride,
    num nu, long num_terms)
{
    num h = 1, accum = data[0];
    long m;
    
    for(m = 1; m < num_terms; m++)
    {
        h *= nu;
        accum += h * data[stride * extension(N, -m)];
    }
    
    return accum;
}

/**
 * \brief Gaussian convolution with Alvarez-Mazorra
 * \param dest the output convolved data
 * \param src the data to be convolved, modified in-place if src = dest
 * \param N number of elements
 * \param stride the stride between successive samples
 * \param sigma the standard deviation of the Gaussian in pixels
 * \param K number of timesteps, more steps implies better accuracy
 * \param tol for symmetric boundaries, accuracy in evaluating left sum
 * \param use_adjusted_q if nonzero, use proposed regression for q
 *
 * Implements the fast approximate Gaussian convolution algorithm of Alvarez
 * and Mazorra, where the Gaussian is approximated by a cascade of first-order
 * recursive filters.  Boundaries are handled with half-sample symmetric 
 * extension, and \c tol specifies the accuracy in approximating an infinite 
 * sum on the left boundary.
 * 
 * Gaussian convolution is approached as approximating the heat equation and 
 * each timestep is performed with an efficient recursive computation.  Using
 * more steps yields a more accurate approximation of the Gaussian.  
 * 
 * Reasonable values for the parameters are \c K = 4, \c tol = 1e-3.
 *
 * Reference:
 * Alvarez, Mazorra, "Signal and Image Restoration using Shock Filters and
 * Anisotropic Diffusion," SIAM J. on Numerical Analysis, vol. 31, no. 2, 
 * pp. 590-605, 1994.
 */
void am_gaussian_conv(num *dest, const num *src, long N, long stride,
    num sigma, int K, num tol, int use_adjusted_q)
{
    const long stride_N = stride * N;
    double q, lambda, dnu;
    num nu, scale;
    long i, M;
    int pass;
    
    assert(dest && src && N > 0 && stride != 0 && sigma > 0
        && K > 0 && tol > 0);
    
    if(use_adjusted_q)
        q = sigma * (1.0 + (0.3148 * K + 0.6121) 
                / ((K + 0.8268) * (K + 0.8268)));
    else
        q = sigma;
    
    lambda = (q * q) / (2.0 * K);
    dnu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);    
    nu = (num)dnu;
    M = (long)ceil(log((1.0 - dnu)*tol) / log(dnu));
    scale = (num)(pow(dnu / lambda, K));
    
    for(i = 0; i < stride_N; i += stride)
        dest[i] = src[i] * scale;
    
    for(pass = 0; pass < K; pass++)
    {
        dest[0] = am_left_boundary(dest, N, stride, nu, M);
        
        for(i = stride; i < stride_N; i += stride)  /* Causal filter */
            dest[i] += nu * dest[i - stride];
        
        i -= stride;        
        dest[i] /= (1 - nu);
        
        for(; i > 0; i -= stride)                   /* Anticausal filter */
            dest[i - stride] += nu * dest[i];
    }
    
    return;
}

/**
 * \brief 2D Gaussian convolution with Alvarez-Mazorra
 * \param dest the output convolved data
 * \param src the data to be convolved, modified in-place if src = dest
 * \param width, height, num_channels the image dimensions
 * \param sigma the standard deviation of the Gaussian in pixels
 * \param K number of timesteps, more steps implies better accuracy
 * \param tol for symmetric boundaries, accuracy in evaluating left sum
 * \param use_adjusted_q if nonzero, use proposed regression for q
 */
void am_gaussian_conv_image(num *dest, const num *src, 
    int width, int height, int num_channels,
    num sigma, int K, num tol, int use_adjusted_q)
{
    const long num_pixels = ((long)width) * ((long)height);
    int x, y, channel;
    
    assert(dest && src && num_pixels > 0 && sigma > 0
        && K > 0 && tol > 0);
    
    for(channel = 0; channel < num_channels; channel++,
        dest += num_pixels, src += num_pixels)
    {
        num *dest_y = dest;
        const num *src_y = src;
            
        for(y = 0; y < height; y++, dest_y += width, src_y += width)
            am_gaussian_conv(dest_y, src_y, width, 1, 
                sigma, K, tol, use_adjusted_q);
        
        for(x = 0; x < width; x++)
            am_gaussian_conv(dest + x, dest + x, height, width,
                sigma, K, tol, use_adjusted_q);
    }
    
    return;
}
