/**
 * \file filter_util.h
 * \brief Filtering utility functions
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

#ifndef _FILTER_UTIL_H_
#define _FILTER_UTIL_H_

#include "num.h"

#ifdef __GNUC__
__attribute__((pure,unused))
#endif
/**
 * \brief Half-sample symmetric boundary extension
 * \param N signal length
 * \param n requested sample, possibly outside {0,...,N-1}
 * \return reflected sample in {0,...,N-1}
 */
static long extension(long N, long n)
{
    while(1)        
        if(n < 0)
            n = -1 - n;         /* Reflect over n = -1/2    */
        else if(n >= N)
            n = 2*N - 1 - n;    /* Reflect over n = N - 1/2 */
        else
            break;
        
    return n;
}

void recursive_filter_impulse(num *h, long N, 
    const num *b, int p, const num *a, int q);

void init_recursive_filter(num *dest, const num *src, long N, long stride,
    const num *b, int p, const num *a, int q, 
    num sum, num tol, long max_iter);

void copy_array(num *dest, const num *src, long N, long stride);

void scale_array(num *dest, const num *src,
    long N, long stride, num scale);

#endif /* _FILTER_UTIL_H_ */
