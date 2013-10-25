/**
 * \file gaussian_conv_deriche.c
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

#include "gaussian_conv_deriche.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "filter_util.h"
#include "complex_arith.h"
#include "gaussian_short_conv.h"

#ifndef M_SQRT2PI
/** \brief The constant sqrt(2 pi) */
#define M_SQRT2PI   2.50662827463100050241576528481104525
#endif

static void make_filter(num *result_b, num *result_a, 
    const complex *alpha, const complex *beta, int K, double sigma);

/**
 * \brief Precompute coefficients for Deriche's Gaussian approximation
 * \param c         deriche_coeffs pointer to hold precomputed coefficients
 * \param sigma     Gaussian standard deviation
 * \param K         filter order = 2, 3, or 4
 * \param tol       accuracy for filter initailization at the boundaries
 * 
 * This routine precomputes the recursive filter coefficients for
 * Deriche's Gaussian convolution approximation.
 */
void deriche_precomp(deriche_coeffs *c, double sigma, int K, num tol)
{    
    static const complex alpha[DERICHE_MAX_K - DERICHE_MIN_K + 1][4] = {
        {{0.48145, 0.971}, {0.48145, -0.971}},
        {{-0.44645, 0.5105}, {-0.44645, -0.5105}, {1.898, 0}},
        {{0.84, 1.8675}, {0.84, -1.8675}, 
            {-0.34015, -0.1299}, {-0.34015, 0.1299}}
        };
    static const complex lambda[DERICHE_MAX_K - DERICHE_MIN_K + 1][4] = {
        {{1.26, 0.8448}, {1.26, -0.8448}},
        {{1.512, 1.475}, {1.512, -1.475}, {1.556, 0}},
        {{1.783, 0.6318}, {1.783, -0.6318}, 
            {1.723, 1.997}, {1.723, -1.997}}
        };
    complex beta[DERICHE_MAX_K];
    
    int k;
    double accum, accum_denom;
    
    assert(c && sigma > 0 && DERICHE_VALID_K(K) && 0 < tol && tol < 1);
    
    for(k = 0; k < K; k++)
    {
        double temp = exp(-lambda[K - DERICHE_MIN_K][k].real / sigma);
        beta[k] = make_complex(
            -temp * cos(lambda[K - DERICHE_MIN_K][k].imag / sigma),
            temp * sin(lambda[K - DERICHE_MIN_K][k].imag / sigma));
    }
    
    /* Compute the causal filter coefficients */
    make_filter(c->b_causal, c->a, alpha[K - DERICHE_MIN_K], beta, K, sigma);
    
    /* Numerator coefficients of the anticausal filter */
    c->b_anticausal[0] = (num)(0.0);
    
    for(k = 1; k < K; k++)
        c->b_anticausal[k] = c->b_causal[k] - c->a[k] * c->b_causal[0];
    
    c->b_anticausal[K] = -c->a[K] * c->b_causal[0];
    
    /* Impulse response sums */
    for(k = 1, accum_denom = 1.0; k <= K; k++)
        accum_denom += c->a[k];
    
    for(k = 0, accum = 0.0; k < K; k++)
        accum += c->b_causal[k];
    
    c->sum_causal = (num)(accum / accum_denom);
    
    for(k = 1, accum = 0.0; k <= K; k++)
        accum += c->b_anticausal[k];
    
    c->sum_anticausal = (num)(accum / accum_denom);
    
    c->sigma = (num)sigma;
    c->K = K;    
    c->tol = tol;
    c->max_iter = (num)ceil(10 * sigma);
    return;
}

/**
 * \brief Make Deriche filter from alpha and beta coefficients 
 * \param result_b      resulting numerator filter coefficients
 * \param result_a      resulting denominator filter coefficients
 * \param alpha, beta   input coefficients
 * \param K             number of terms
 * \param sigma         Gaussian sigma parameter
 * \ingroup deriche_gaussian
 * 
 * This routine performs the algebraic rearrangement 
 * \f[ \sum_{k=0}^{K-1}\frac{\alpha_k}{1+\beta_k z^{-1}}=\frac{1}{\sqrt{2\pi
\sigma^2}}\frac{\sum_{k=0}^{K-1}b_k z^{-k}}{1+\sum_{k=1}^{K}a_k z^{-k}} \f]
 * to obtain the numerator and denominator coefficients for the causal filter
 * in Deriche's Gaussian approximation.
 * 
 * The routine initializes b/a as the 0th term,
 * \f[ \frac{b(z)}{a(z)} = \frac{\alpha_0}{1 + \beta_0 z^{-1}}, \f]
 * then the kth term is added according to
 * \f[ \frac{b(z)}{a(z)}\leftarrow\frac{b(z)}{a(z)}+\frac{\alpha_k}{1+\beta_k
z^{-1}}=\frac{b(z)(1+\beta_kz^{-1})+a(z)\alpha_k}{a(z)(1+\beta_kz^{-1})}. \f]
 */
static void make_filter(num *result_b, num *result_a, 
    const complex *alpha, const complex *beta, int K, double sigma)
{
    const double denom = sigma * M_SQRT2PI;
    complex b[DERICHE_MAX_K], a[DERICHE_MAX_K + 1];
    int k, j;
        
    b[0] = alpha[0];    /* Initialize b/a = alpha[0] / (1 + beta[0] z^-1) */
    a[0] = make_complex(1, 0);
    a[1] = beta[0];
    
    for(k = 1; k < K; k++)
    {   /* Add kth term, b/a += alpha[k] / (1 + beta[k] z^-1) */
        b[k] = c_mul(beta[k], b[k-1]);
        
        for(j = k - 1; j > 0; j--)
            b[j] = c_add(b[j], c_mul(beta[k], b[j - 1]));
        
        for(j = 0; j <= k; j++)
            b[j] = c_add(b[j], c_mul(alpha[k], a[j]));
        
        a[k + 1] = c_mul(beta[k], a[k]);
        
        for(j = k; j > 0; j--)
            a[j] = c_add(a[j], c_mul(beta[k], a[j - 1]));
    }
    
    for(k = 0; k < K; k++)
    {
        result_b[k] = (num)(b[k].real / denom);
        result_a[k + 1] = (num)a[k + 1].real;
    }
    
    return;
}

/**
 * \brief Deriche Gaussian convolution
 * \param c         coefficients precomputed by deriche_precomp()
 * \param dest      output convolved data
 * \param buffer    workspace array with space for at least 2*N elements
 * \param src       data to be convolved
 * \param N         number of samples
 * \param stride    stride between successive samples
 */
void deriche_gaussian_conv(deriche_coeffs c,
    num *dest, num *buffer, const num *src, long N, long stride)
{
    const long stride_2 = stride * 2;
    const long stride_3 = stride * 3;
    const long stride_4 = stride * 4;
    const long stride_N = stride * N;
    num *buffer_l, *buffer_r;
    long i, n;
    
    assert(dest && buffer && src && buffer != src && N > 0 && stride != 0);
    
    if(N <= 4)
    {
        gaussian_short_conv(dest, src, N, stride, c.sigma);
        return;
    }
    
    /* Divide buffer into two buffers each of length N */
    buffer_l = buffer;
    buffer_r = buffer + N;
    
    /* Causal filter */
    init_recursive_filter(buffer_l, src, N, stride,
        c.b_causal, c.K - 1, c.a, c.K, c.sum_causal, c.tol, c.max_iter);
    
    switch(c.K)
    {
    case 2:
        for(n = 2, i = stride_2; n < N; n++, i += stride)
            buffer_l[n] = c.b_causal[0] * src[i]
                + c.b_causal[1] * src[i - stride]
                - c.a[1] * buffer_l[n - 1] 
                - c.a[2] * buffer_l[n - 2];
        break;
    case 3:
        for(n = 3, i = stride_3; n < N; n++, i += stride)
            buffer_l[n] = c.b_causal[0] * src[i]
                + c.b_causal[1] * src[i - stride]
                + c.b_causal[2] * src[i - stride_2]
                - c.a[1] * buffer_l[n - 1] 
                - c.a[2] * buffer_l[n - 2]
                - c.a[3] * buffer_l[n - 3];
        break;
    case 4:
        for(n = 4, i = stride_4; n < N; n++, i += stride)
            buffer_l[n] = c.b_causal[0] * src[i]
                + c.b_causal[1] * src[i - stride]
                + c.b_causal[2] * src[i - stride_2]
                + c.b_causal[3] * src[i - stride_3]
                - c.a[1] * buffer_l[n - 1] 
                - c.a[2] * buffer_l[n - 2]
                - c.a[3] * buffer_l[n - 3]
                - c.a[4] * buffer_l[n - 4];
        break;
    }
    
    /* Anticausal filter */
    init_recursive_filter(buffer_r, src + stride_N - stride, N, -stride,
        c.b_anticausal, c.K, c.a, c.K, c.sum_anticausal, c.tol, c.max_iter);
    
    switch(c.K)
    {
    case 2:
        for(n = 2, i = stride_N - stride_3; n < N; n++, i -= stride)
            buffer_r[n] = c.b_anticausal[1] * src[i + stride]
                + c.b_anticausal[2] * src[i + stride_2]
                - c.a[1] * buffer_r[n - 1] 
                - c.a[2] * buffer_r[n - 2];
        break;
    case 3:
        for(n = 3, i = stride_N - stride_4; n < N; n++, i -= stride)
            buffer_r[n] = c.b_anticausal[1] * src[i + stride]
                + c.b_anticausal[2] * src[i + stride_2]
                + c.b_anticausal[3] * src[i + stride_3]
                - c.a[1] * buffer_r[n - 1] 
                - c.a[2] * buffer_r[n - 2] 
                - c.a[3] * buffer_r[n - 3];
        break;
    case 4:
        for(n = 4, i = stride_N - stride * 5; n < N; n++, i -= stride)
            buffer_r[n] = c.b_anticausal[1] * src[i + stride]
                + c.b_anticausal[2] * src[i + stride_2]
                + c.b_anticausal[3] * src[i + stride_3]
                + c.b_anticausal[4] * src[i + stride_4]
                - c.a[1] * buffer_r[n - 1] 
                - c.a[2] * buffer_r[n - 2] 
                - c.a[3] * buffer_r[n - 3]
                - c.a[4] * buffer_r[n - 4];
        break;
    }
    
    for(n = 0, i = 0; n < N; n++, i += stride)
        dest[i] = buffer_l[n] + buffer_r[N - n - 1];
    
    return;
}

/**
 * \brief Deriche Gaussian 2D convolution
 * \param c             coefficients precomputed by deriche_precomp()
 * \param dest          output convolved data
 * \param buffer        array with at least 2*max(width,height) elements
 * \param src           data to be convolved
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 */
void deriche_gaussian_conv_image(deriche_coeffs c,
    num *dest, num *buffer, const num *src, 
    int width, int height, int num_channels)
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
            deriche_gaussian_conv(c,
                dest_y, buffer, src_y, width, 1);
        
        for(x = 0; x < width; x++)
            deriche_gaussian_conv(c, 
                dest + x, buffer, dest + x, height, width);
    }
    
    return;
}
