/**
 * \file gaussian_conv_vyv.c
 * \brief Vliet-Young-Verbeek approximate Gaussian convolution
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
 * You should have received a copy of these licenses along this program. 
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "filter_util.h"
#include "complex_arith.h"
#include "invert_matrix.h"
#include "gaussian_short_conv.h"
#include "gaussian_conv_vyv.h"
#include <cilk/cilk.h>
#include <stdio.h>

/** \brief Number of newton iterations used to determine q */
#define YVY_NUM_NEWTON_ITERATIONS       6

/** 
 * \brief Compute the variance of the impulse response
 * \param poles0    unscaled pole locations
 * \param q         rescaling parameter
 * \param K         number of poles
 * \return variance acheived by poles = poles0^(1/q)
 * \ingroup vyv_gaussian
 */
static double variance(const complex *poles0, int K, double q)
{
    complex sum = {0, 0};
    int k;
    
    for(k = 0; k < K; k++)
    {
        complex z = c_real_pow(poles0[k], 1/q), denom = z;
        denom.real -= 1;
        /* Compute sum += z / (z - 1)^2 */
        sum = c_add(sum, c_div(z, c_mul(denom, denom)));
    }
    
    return 2 * sum.real;
}

/**
 * \brief Derivative of variance with respect to q
 * \param poles0    unscaled pole locations
 * \param q         rescaling parameter
 * \param K         number of poles
 * \return derivative of variance with respect to q
 * \ingroup vyv_gaussian
 */
static double dq_variance(const complex *poles0, int K, double q)
{
    complex sum = {0, 0};
    int k;
    
    for(k = 0; k < K; k++)
    {   
        complex z = c_real_pow(poles0[k], 1/q), w = z, denom = z;
        w.real += 1;
        denom.real -= 1;
        /* Compute sum += z log(z) (z + 1) / (z - 1)^3 */
        sum = c_add(sum, c_div(c_mul(c_mul(z, c_log(z)), w), 
            c_real_pow(denom, 3)));
    }
    
    return (2 / q) * sum.real;
}

/** 
 * \brief Compute q for a desired sigma using Newton's method 
 * \param poles0    unscaled pole locations
 * \param K         number of poles
 * \param sigma     the desired sigma
 * \param q0        initial estimate of q
 * \return refined value of q
 * \ingroup vyv_gaussian
 */
static double compute_q(const complex *poles0, int K, 
    double sigma, double q0)
{
    double sigma2 = sigma * sigma;
    double q = q0;
    int i;
    
    for(i = 0; i < YVY_NUM_NEWTON_ITERATIONS; i++)
        q -= (variance(poles0, K, q) - sigma2) 
            / dq_variance(poles0, K, q);
    
    return q;
}

/**
 * \brief Expand pole product
 * \param c         resulting filter coefficients
 * \param poles     pole locations
 * \param K         number of poles
 * \ingroup vyv_gaussian
 * 
 * This routine expands the product to obtian the filter coefficients:
 * \f[ \prod_{k=0}^{K-1}\frac{\mathrm{poles}[k]-1}{\mathrm{poles}[k]-z^{-1}}
 = \frac{c[0]}{1+\sum_{k=1}^K c[k] z^{-k}} \f]
 */
static void expand_pole_product(double *c, const complex *poles, int K)
{
    complex denom[VYV_MAX_K + 1];
    int k, j;
    
    assert(K <= VYV_MAX_K);
    denom[0] = poles[0];
    denom[1] = make_complex(-1, 0);
    
    for(k = 1; k < K; k++)
    {
        denom[k + 1] = c_neg(denom[k]);
        
        for(j = k; j > 0; j--)
            denom[j] = c_sub(c_mul(denom[j], poles[k]), denom[j - 1]);
        
        denom[0] = c_mul(denom[0], poles[k]);
    }
    
    for(k = 1; k <= K; k++)
        c[k] = c_div(denom[k], denom[0]).real;
    
    for(c[0] = 1, k = 1; k <= K; k++)
        c[0] += c[k];
    
    return;
}

/**
 * \brief Precomputations for Vliet-Young-Verbeek Gaussian approximation
 * \param c         coefficients structure to populate
 * \param sigma     the standard deviation of the Gaussian in pixels
 * \param K         filter order = 3, 4, or 5
 * \param tol       accuracy for initialization on the left boundary
 */
void vyv_precomp(vyv_coeffs *c, double sigma, int K, num tol)
{
    static const complex poles0[VYV_MAX_K - VYV_MIN_K + 1][5] = {
        {{1.4165, 1.00829}, {1.4165, -1.00829}, {1.86543, 0}}, 
        {{1.13228, 1.28114}, {1.13228, -1.28114},
            {1.78534, 0.46763}, {1.78534, -0.46763}},
        {{0.8643, 1.45389}, {0.8643, -1.45389},
            {1.61433, 0.83134}, {1.61433, -0.83134}, {1.87504, 0}}
        };
    complex poles[VYV_MAX_K];
    double q, filter[VYV_MAX_K + 1];
    double A[VYV_MAX_K * VYV_MAX_K], inv_A[VYV_MAX_K * VYV_MAX_K];
    int i, j, matrix_size;
    
    assert(c && sigma > 0 && VYV_VALID_K(K) && tol >= 0);
    /* Initial estimate of q */
    q = sigma / 2;
    /* Compute accurate value of q using Newton's method */
    q = compute_q(poles0[K - VYV_MIN_K], K, sigma, q);
    
    /* Compute the filter coefficients */
    for(i = 0; i < K; i++)
        poles[i] = c_real_pow(poles0[K - VYV_MIN_K][i], 1/q);
    
    expand_pole_product(filter, poles, K);
    
    /* Compute matrix for handling the right boundary */
    for(i = 0, matrix_size = K * K; i < matrix_size; i++)
        A[i] = 0.0;
        
    for(i = 0; i < K; i++)
        for(A[i + K * i] = 1.0, j = 1; j <= K; j++)
            A[i + K * extension(K, i + j)] += filter[j];
    
    invert_matrix(inv_A, A, K);
    
    for(i = 0; i < matrix_size; i++)
        inv_A[i] *= filter[0];
    
    /* Store precomputations in coeffs struct */
    for(i = 0; i <= K; i++)
        c->filter[i] = (num)filter[i];
    
    for(i = 0; i < matrix_size; i++)
        c->M[i] = (num)inv_A[i];
    
    c->K = K;
    c->sigma = (num)sigma;
    c->tol = tol;
    c->max_iter = (num)(10 * sigma);
    return;
}

/**
 * \brief Gaussian convolution Vliet-Young-Verbeek approximation
 * \param c         vyv_coeffs created by vyv_precomp()
 * \param dest      output convolved data
 * \param src       data to be convolved, modified in-place if src = dest 
 * \param N         number of samples
 * \param stride    stride between successive samples
 */
void vyv_gaussian_conv(vyv_coeffs c,
    num *dest, const num *src, long N, long stride)
{
    const long stride_2 = stride * 2;
    const long stride_3 = stride * 3;
    const long stride_4 = stride * 4;
    const long stride_5 = stride * 5;
    const long stride_N = stride * N;
    num q[VYV_MAX_K];
    long i;
    int m, n;
    
    assert(dest && src && N > 0 && stride != 0);
    
    if(N <= 4)
    {
        gaussian_short_conv(dest, src, N, stride, c.sigma);
        return;
    }
    
    init_recursive_filter(q, src, N, stride,
        c.filter, 0, c.filter, c.K, 1.0f, c.tol, c.max_iter);

    for(m = 0; m < c.K; m++)
        dest[stride * m] = q[m];
    
    /* Apply causal recursive filter on samples n = K,...,N-1 */
    switch(c.K)
    {
    case 3:
        for(i = stride_3; i < stride_N-2*stride; i += 2*stride) {
              dest[i]   =        c.filter[0]*src[i] - (c.filter[1] * dest[i - stride] +
                               c.filter[2] * dest[i - stride_2] +
                               c.filter[3] * dest[i - stride_3]);
              dest[i + stride]   = c.filter[0]*src[i+stride] - (c.filter[1] * dest[i] +
                               c.filter[2] * dest[i - stride] +
                               c.filter[3] * dest[i - stride_2]);

        }
        for(; i < stride_N; i += stride) {
              dest[i]   =        c.filter[0]*src[i] - (c.filter[1] * dest[i - stride] +
                               c.filter[2] * dest[i - stride_2] +
                               c.filter[3] * dest[i - stride_3]);
        }

        break;
    case 4:
        printf("case 4\n");
        for(i = stride_4; i < stride_N; i += stride)
            dest[i] = c.filter[0] * src[i] 
                    - c.filter[1] * dest[i - stride] 
                    - c.filter[2] * dest[i - stride_2] 
                    - c.filter[3] * dest[i - stride_3]
                    - c.filter[4] * dest[i - stride_4];
        break;
    case 5:
        printf("case 5\n");
        for(i = stride_5; i < stride_N; i += stride)
            dest[i] = c.filter[0] * src[i] 
                    - c.filter[1] * dest[i - stride] 
                    - c.filter[2] * dest[i - stride_2] 
                    - c.filter[3] * dest[i - stride_3]
                    - c.filter[4] * dest[i - stride_4]
                    - c.filter[5] * dest[i - stride_5];
        break;
    }
    
    /* Handle right boundary, multiply M with last K samples */
    
    for(m = 0; m < c.K; m++)
        q[m] = dest[stride_N - stride * (c.K - m)];
    
    for(m = 0; m < c.K; m++)
    {
        num accum = (num)0;
        
        for(n = 0; n < c.K; n++)
            accum += c.M[m + c.K*n] * q[n];
        
        dest[stride_N - stride * (c.K - m)] = accum;
    }
        
    /* Apply anticausal recursive filter on samples n = N-1-K,...,0 */
    switch(c.K)
    {
    case 3:

        for(i = stride_N - stride_4; i >= 2*stride; i -= 2*stride){
             dest[i] = c.filter[0] * dest[i]
                    - c.filter[1] * dest[i + stride]
                    - c.filter[2] * dest[i + stride_2]
                    - c.filter[3] * dest[i + stride_3];

             dest[i - stride] = c.filter[0] * dest[i - stride]
                    - c.filter[1] * dest[i]
                    - c.filter[2] * dest[i + stride]
                    - c.filter[3] * dest[i + stride_2];

        }
        for(; i >= 0; i -= stride){
             dest[i] = c.filter[0] * dest[i]
                    - c.filter[1] * dest[i + stride]
                    - c.filter[2] * dest[i + stride_2]
                    - c.filter[3] * dest[i + stride_3];
        }

        break;
    case 4:
        for(i = stride_N - stride_5; i >= 0; i -= stride)
            dest[i] = c.filter[0] * dest[i]
                    - c.filter[1] * dest[i + stride]
                    - c.filter[2] * dest[i + stride_2]
                    - c.filter[3] * dest[i + stride_3]
                    - c.filter[4] * dest[i + stride_4];
        break;
    case 5:
        for(i = stride_N - stride * 6; i >= 0; i -= stride)
            dest[i] = c.filter[0] * dest[i]
                    - c.filter[1] * dest[i + stride]
                    - c.filter[2] * dest[i + stride_2]
                    - c.filter[3] * dest[i + stride_3]
                    - c.filter[4] * dest[i + stride_4]
                    - c.filter[5] * dest[i + stride_5];
        break;
    }

    return;
}

/**
 * \brief 2D Gaussian convolution Vliet-Young-Verbeek approximation
 * \param c             vyv_coeffs created by vyv_precomp()
 * \param dest          destination image data
 * \param src           source image data, modified in-place if src = dest
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 */
void vyv_gaussian_conv_image(vyv_coeffs c, num *dest, const num *src,
    int width, int height, int num_channels)
{
    const long num_pixels = ((long)width) * ((long)height);
    int x, y, channel;
    
    assert(dest && src && num_pixels > 0);
    
    for(channel = 0; channel < num_channels; channel++,
        dest += num_pixels, src += num_pixels)
    {
        num *dest_y = dest;
        const num *src_y = src;
            
        for(y = 0; y < height; y++, dest_y += width, src_y += width)
            vyv_gaussian_conv(c, dest_y, src_y, width, 1);
        
        for(x = 0; x < width; x++)
            vyv_gaussian_conv(c, dest + x, dest + x, height, width);
    }
    
    return;
}
