/**
 * \file strategy_gaussian_conv.c
 * \brief Suite of different 1D Gaussian convolution methods 
 * \author Pascal Getreuer <getreuer@gmail.com>
 * 
 * Copyright (c) 2013, Pascal Getreuer
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

#include "strategy_gaussian_conv.h"
#include <stdlib.h>
#include <string.h>
#include "gaussian_conv_fir.h"
#include "gaussian_conv_dct.h"
#include "gaussian_conv_box.h"
#include "gaussian_conv_ebox.h"
#include "gaussian_conv_sii.h"
#include "gaussian_conv_am.h"
#include "gaussian_conv_deriche.h"
#include "gaussian_conv_vyv.h"

struct gconv_
{
    num *dest;                  /**< destination array                 */
    const num *src;             /**< sourcec array                     */
    num *buffer;                /**< workspace memory (if needed)      */
    long N;                     /**< number of samples                 */
    long stride;                /**< stride between successive samples */
    const char *algo;           /**< algorithm name                    */
    double sigma;               /**< Gaussian standard deviation       */
    int K;                      /**< steps/filter order parameter      */
    num tol;                    /**< accuracy parameter                */
    void *coeffs;               /**< algorithm-specific data           */
    void (*execute)(gconv*);    /**< algorithm execution function      */
    void (*free)(gconv*);       /**< algorithm clean up function       */
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static int _fir_plan(gconv *g);
static int _dct_plan(gconv *g);
static int _box_plan(gconv *g);
static int _ebox_plan(gconv *g);
static int _sii_plan(gconv *g);
static int _am_plan(gconv *g);
static int _am_original_plan(gconv *g);
static int _deriche_plan(gconv *g);
static int _vyv_plan(gconv *g);
#endif

/**
 * \brief Plan a 1D Gaussian convolution
 * \param dest      destination array
 * \param src       source array
 * \param N         number of samples
 * \param stride    stride between successive samples
 * \param algo      Gaussian convolution algorithm
 * \param sigma     Gaussian standard deviation
 * \param K         algorithm steps or filter order parameter
 * \param tol       algorithm accuracy parameter
 * \return gconv pointer or NULL on failure
 */
gconv* gconv_plan(num *dest, const num *src, long N, long stride,
    const char *algo, double sigma, int K, num tol)
{
    struct gconv_algo_entry
    {
        const char *name;
        int (*plan)(gconv*);
    } algos[] = {
        {"fir",         _fir_plan},
        {"dct",         _dct_plan},
        {"box",         _box_plan},
        {"ebox",        _ebox_plan},
        {"sii",         _sii_plan},
        {"am",          _am_plan},
        {"am_original", _am_original_plan},
        {"deriche",     _deriche_plan},
        {"vyv",         _vyv_plan}
        };
    gconv *g;
    size_t i;
    
    if(!dest || !src || N <= 0 || sigma <= 0 || tol < 0 || tol > 1
        || !(g = (gconv *)malloc(sizeof(gconv))))
        return NULL;
    
    g->dest = dest;
    g->src = src;
    g->buffer = NULL;
    g->N = N;
    g->stride = stride;
    g->algo = algo;
    g->K = K;
    g->sigma = sigma;
    g->tol = tol;
    g->coeffs = NULL;
    
    for(i = 0; i < sizeof(algos)/sizeof(*algos); i++)
        if(!strcmp(algo, algos[i].name))
        {
            if(algos[i].plan(g))
                return g;
            else
                break;
        }
    
    gconv_free(g);    
    return NULL;
}

/**
 * \brief Execute a 1D Gaussian convolution
 * \param g     gconv* created by gconv_plan()
 */
void gconv_execute(gconv *g)
{
    g->execute(g);
}

/**
 * \brief Free memory associated with a gconv
 * \param g     gconv* created by gconv_plan()
 */
void gconv_free(gconv *g)
{
    if(g)
    {
        if(g->free)
            g->free(g);
        
        free(g);
    }
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* FIR filtering */
static void _fir_execute(gconv *g);
static void _fir_free(gconv *g);

static int _fir_plan(gconv *g)
{
    g->execute = _fir_execute;
    g->free = _fir_free;
    return (g->coeffs = malloc(sizeof(fir_coeffs)))
        && fir_precomp((fir_coeffs *)g->coeffs, g->sigma, g->tol);
}

static void _fir_execute(gconv *g)
{
    fir_gaussian_conv(*((fir_coeffs *)g->coeffs),
        g->dest, g->src, g->N, g->stride);
}

static void _fir_free(gconv *g)
{
    if(g->coeffs)
    {
        fir_free((fir_coeffs *)g->coeffs);
        free(g->coeffs);
    }
}

/* DCT (discrete cosine transform) based convolution */
static void _dct_execute(gconv *g);
static void _dct_free(gconv *g);

static int _dct_plan(gconv *g)
{
    g->execute = _dct_execute;
    g->free = _dct_free;
    return (g->coeffs = malloc(sizeof(dct_coeffs)))
        && dct_precomp((dct_coeffs *)g->coeffs, g->dest, g->src,
            g->N, g->stride, g->sigma);
}

static void _dct_execute(gconv *g)
{
    dct_gaussian_conv(*((dct_coeffs *)g->coeffs));
}

static void _dct_free(gconv *g)
{
    if(g->coeffs)
    {
        dct_free((dct_coeffs *)g->coeffs);
        free(g->coeffs);
    }
}

/* Box filtering */
static void _box_execute(gconv *g);
static void _box_free(gconv *g);

static int _box_plan(gconv *g)
{
    g->execute = _box_execute;
    g->free = _box_free;
    return ((g->K > 0) && (g->buffer = malloc(sizeof(num)*g->N))) ? 1 : 0;
}

static void _box_execute(gconv *g)
{
    box_gaussian_conv(g->dest, g->buffer, g->src, g->N, g->stride, 
        g->sigma, g->K);
}

static void _box_free(gconv *g)
{
    if(g->buffer)
        free(g->buffer);
} 

/* Extended box filter */
static void _ebox_execute(gconv *g);
static void _ebox_free(gconv *g);

static int _ebox_plan(gconv *g)
{
    g->execute = _ebox_execute;
    g->free = _ebox_free;
    
    if((g->K <= 0) || !(g->buffer = malloc(sizeof(num) * g->N))
        || !(g->coeffs = malloc(sizeof(ebox_coeffs))))
        return 0;
    
    ebox_precomp((ebox_coeffs *)g->coeffs, g->sigma, g->K);
    return 1;
}

static void _ebox_execute(gconv *g)
{    
    ebox_gaussian_conv(*((ebox_coeffs *)g->coeffs), 
        g->dest, g->buffer, g->src, g->N, g->stride);
}

static void _ebox_free(gconv *g)
{
    if(g->buffer)
        free(g->buffer);
    if(g->coeffs)
        free(g->coeffs);
}

/* SII Stacked integral images */
static void _sii_execute(gconv *g);
static void _sii_free(gconv *g);

static int _sii_plan(gconv *g)
{
    g->execute = _sii_execute;
    g->free = _sii_free;
    
    if(!SII_VALID_K(g->K) || !(g->coeffs = malloc(sizeof(sii_coeffs))))
        return 0;
    
    sii_precomp((sii_coeffs *)g->coeffs, g->sigma, g->K);
    
    return (g->buffer = malloc(sizeof(num) * sii_buffer_size(
        *((sii_coeffs *)g->coeffs), g->N))) ? 1 : 0;
}

static void _sii_execute(gconv *g)
{
    sii_gaussian_conv(*((sii_coeffs *)g->coeffs),
        g->dest, g->buffer, g->src, g->N, g->stride);
}

static void _sii_free(gconv *g)
{
    if(g->buffer)
        free(g->buffer);
    if(g->coeffs)
        free(g->coeffs);
}

/* Alvarez-Mazorra with proposed regression on parameter q */
static void _am_execute(gconv *g);

static int _am_plan(gconv *g)
{
    g->execute = _am_execute;
    g->free = NULL;
    return (g->K > 0);
}

static void _am_execute(gconv *g)
{
    am_gaussian_conv(g->dest, g->src, g->N, 1, 
        g->sigma, g->K, g->tol, 1);
}

/* Alvarez-Mazorra, original algorithm */
static void _am_original_execute(gconv *g);

static int _am_original_plan(gconv *g)
{
    g->execute = _am_original_execute;
    g->free = NULL;
    return 1;
}

static void _am_original_execute(gconv *g)
{
    am_gaussian_conv(g->dest, g->src, g->N, 1, 
        g->sigma, g->K, g->tol, 0);
}

/* Deriche recursive filtering */
static void _deriche_execute(gconv *g);
static void _deriche_free(gconv *g);

static int _deriche_plan(gconv *g)
{
    g->execute = _deriche_execute;
    g->free = _deriche_free;
    
    if(!DERICHE_VALID_K(g->K)
        || !(g->buffer = malloc(sizeof(num) * 2 * g->N))
        || !(g->coeffs = malloc(sizeof(deriche_coeffs))))
        return 0;
    
    deriche_precomp((deriche_coeffs *)g->coeffs, g->sigma, g->K, g->tol);    
    return 1;
}

static void _deriche_execute(gconv *g)
{
    deriche_gaussian_conv(*((deriche_coeffs *)g->coeffs),
        g->dest, g->buffer, g->src, g->N, g->stride);
}

static void _deriche_free(gconv *g)
{
    if(g->buffer)
        free(g->buffer);
    if(g->coeffs)
        free(g->coeffs);
}

/* Vliet-Young-Verbeek recursive filtering */
static void _vyv_execute(gconv *g);
static void _vyv_free(gconv *g);

static int _vyv_plan(gconv *g)
{
    g->execute = _vyv_execute;
    g->free = _vyv_free;
    
    if(!VYV_VALID_K(g->K) || !(g->coeffs = malloc(sizeof(vyv_coeffs))))
        return 0;
    
    vyv_precomp((vyv_coeffs *)g->coeffs, g->sigma, g->K, g->tol);    
    return 1;
}

static void _vyv_execute(gconv *g)
{
    vyv_gaussian_conv(*((vyv_coeffs *)g->coeffs),
        g->dest, g->src, g->N, g->stride);
}

static void _vyv_free(gconv *g)
{
    if(g->coeffs)
        free(g->coeffs);
}    
#endif
