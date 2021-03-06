/**
 * \file gaussian_conv_fir.c
 * \brief Gaussian convolution using FIR filters
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

#include "gaussian_conv_fir.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "filter_util.h"
#include "inverfc_acklam.h"

/**
 * \brief Construct truncated Gaussian filter for FIR convolution
 * \param sigma the Gaussian standard deviation
 * \param r the radius of the filter
 * \return pointer to filter, or NULL on failure
 * \ingroup fir_gaussian
 * 
 * This routine constructs a truncated Gaussian \f$ g^\text{trunc}\f$
 * for a specified radius \c r.  The filter should be released by the caller
 * with \c free().
 */
static num *make_g_trunc(num sigma, long r)
{
    num *g_trunc = NULL;
    
    if((g_trunc = (num *)malloc(sizeof(num) * (r + 1))))
    {
        num accum = g_trunc[0] = 1.0;
        long m;
        
        for(m = 1; m <= r; m++)
        {
            num temp = m / sigma;
            g_trunc[m] = exp(-0.5*temp*temp);
            accum += 2.0 * g_trunc[m];
        }
        
        for(m = 0; m <= r; m++)
            g_trunc[m] /= accum;
    }
    
    return g_trunc;
}

/**
 * \brief Convolution with a symmetric filter
 * \param dest destination (must be distinct from src)
 * \param src source signal
 * \param N signal length
 * \param stride the stride between successive src and dest samples
 * \param h symmetric filter, an array of length r + 1
 * \param r radius of filter h
 * \ingroup fir_gaussian
 * 
 * This routine computes the convolution of \c src and \c h according to
 * \f[ \mathrm{dest[stride*n]} = \sum_{|m| \le r} h_{|m|} \,
 \mathrm{src[stride*(n-m)]}, \f]
 * where \c src is extrapolated with half-sample symmetry.
 */
static void conv_sym(num *dest, const num *src, long N, long stride, 
    const num *h, long r)
{
    long n;
    
    for(n = 0; n < N; n++)
    {        
        num accum = h[0] * src[stride * n];
        long m;
        
        for(m = 1; m <= r; m++)
            accum += h[m] * (src[stride * extension(N, n - m)]
                + src[stride * extension(N, n + m)]);
        
        dest[stride * n] = accum;
    }
    
    return;
}

/**
 * \brief Precompute filter coefficients for FIR filtering
 * \param c         fir_coeffs pointer to hold precomputed coefficients
 * \param sigma     Gaussian standard deviation
 * \param tol       filter accuracy (smaller tol implies larger filter)
 */
int fir_precomp(fir_coeffs *c, double sigma, num tol)
{
    assert(c && sigma > 0.0 && 0.0 < tol && tol < 1.0);    
    c->radius = (long)ceil(M_SQRT2 * sigma * inverfc(0.5 * tol));
    return (c->g_trunc = make_g_trunc(sigma, c->radius)) ? 1 : 0;
}

/**
 * \brief FIR Gaussian convolution
 * \param c         fir_coeffs created by fir_precomp()
 * \param dest      output convolved data (must be distinct from src)
 * \param src       data to be convolved
 * \param N         number of samples
 * \param stride    stride between successive samples
 */
void fir_gaussian_conv(fir_coeffs c, num *dest, const num *src,
    long N, long stride)
{
    assert(c.g_trunc && dest && src && dest != src && N > 0 && stride != 0);
    conv_sym(dest, src, N, stride, c.g_trunc, c.radius);
    return;
}

/** 
 * \brief FIR filtering approximation of 2D Gaussian convolution
 * \param c             fir_coeffs created by fir_precomp()
 * \param dest          destination image
 * \param buffer        array with at least width samples
 * \param src           src image
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 */
void fir_gaussian_conv_image(fir_coeffs c, num *dest, num *buffer,
    const num *src, int width, int height, int num_channels)
{
    const long num_pixels = ((long)width) * ((long)height);
    int x, y, channel;
    
    assert(c.g_trunc && dest && buffer 
        && src && dest != src && num_pixels > 0);
    
    for(channel = 0; channel < num_channels; channel++,
        dest += num_pixels, src += num_pixels)
    {
        num *dest_y = dest;
        const num *src_y = src;
            
        for(x = 0; x < width; x++)
            conv_sym(dest + x, src + x, height, width, c.g_trunc, c.radius);
        
        for(y = 0; y < height; y++, dest_y += width, src_y += width)
        {
            conv_sym(buffer, dest_y, width, 1, c.g_trunc, c.radius);
            memcpy(dest_y, buffer, sizeof(num) * width);
        }
    }
    
    return;
}

/** 
 * \brief Release memory associated with fir_coeffs struct 
 * \param c    fir_coeffs created by fir_precomp()
 */
void fir_free(fir_coeffs *c)
{
    if(c && c->g_trunc)
    {
        free(c->g_trunc);
        c->g_trunc = NULL;
    }
    
    return;
}
