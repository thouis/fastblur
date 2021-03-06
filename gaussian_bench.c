/**
 * \file gaussian_bench.c
 * \brief Benchmark for testing Gaussian convolution algorithms
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
 * You should have received a copy of these licenses along this program. 
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

/**
 * \mainpage
 * \section Overview
 * This work surveys algorithms for computing and efficiently approximating
 * Gaussian convolution,
 * \f[ u(x)=(G_\sigma*f)(x):=\int_{\mathbb{R}^d}G_\sigma(x-y)f(y)\,dy,\f]
 * where f is the input signal, u is the filtered signal, and 
 * \f$ G_\sigma \f$ is the Gaussian with standard deviation \f$ \sigma \f$,
 * \f[ G_\sigma(x)=(2\pi\sigma^2)^{-d/2}\exp\left(-\frac{\lVert x
\rVert_2^2}{2\sigma^2}\right). \f]
 * 
 * \section Algorithms
 *  - \ref fir_gaussian
 *  - \ref dct_gaussian
 *  - \ref box_gaussian
 *  - \ref ebox_gaussian
 *  - \ref sii_gaussian
 *  - \ref am_gaussian
 *  - \ref deriche_gaussian
 *  - \ref vyv_gaussian
 * 
 * \section Usage
 * The gaussian_demo.c program applies 2D convolution using any algorithm 
 * on an image.
\verbatim
Usage: gaussian_conv_demo [options] <input> <output>

Only BMP/JPEG/PNG/TIFF images are supported.

Options:
   -a <algo>     algorithm to use, choices are
                 fir     FIR approximation, tol = kernel accuracy
                 dct     DCT-based convolution
                 box     box filtering, K = # passes
                 ebox    extended box filtering, K = # passes
                 sii     stacked integral images, K = # boxes
                 am      Alvarez-Mazorra recursive filtering,
                         K = # passes, tol = boundary accuracy
                 deriche Deriche recursive filtering,
                         tol = boundary accuracy
                 vyv     Vliet-Young-Verbeek recursive filtering,
                         K = order, tol = boundary accuracy
   -s <number>   sigma, standard deviation of the Gaussian    
   -K <number>   specifies number of steps (box, ebox, sii, am)
   -t <number>   accuracy tolerance (fir, am, deriche, vyv)
\endverbatim
 * 
 * The gaussian_bench.c program tests the speed, accuracy, and impulse 
 * response of any algorithm.
\verbatim
Usage: gaussian_bench [bench type] [options]

Bench type:
   speed         measure computation time
   accuracy      measure L^infty operator norm
   impulse       compute impulse response, written to bench.out

Options:
   -a <algo>     algorithm to use, choices are
                 fir     FIR approximation, tol = kernel accuracy
                 dct     DCT-based convolution
                 box     box filtering, K = # passes
                 ebox    extended box filtering, K = # passes
                 sii     stacked integral images, K = # boxes
                 am      Alvarez-Mazorra recursive filtering,
                         K = # passes, tol = boundary accuracy
                 deriche Deriche recursive filtering,
                         tol = boundary accuracy
                 vyv     Vliet-Young-Verbeek recursive filtering,
                         K = order, tol = boundary accuracy
   -s <number>   sigma, standard deviation of the Gaussian
   -K <number>   specifies number of steps (box, ebox, sii, am)
   -t <number>   accuracy tolerance (fir, am, deriche, vyv)
   -N <number>   signal length
   
   -r <number>   (speed bench) number of runs
   -n <number>   (impulse bench) position of the impulse
\endverbatim
 * 
 * \section Implementation
 * The numeric datatype for computation is #num.  By default, this is double
 * precision.  It can be changed to single precision by defining NUM_SINGLE,
 * see num.h.
 * 
 * The FFTW library is needed for computing DCT transforms in the DCT-based
 * convolution algorithm.  Optionally, libjpeg, libpng, and libtiff can be
 * used to read and write JPEG, PNG, and TIFF image formats (the programs
 * have native support for BMP images).
 * 
 * Please see the readme for further details and compilation instructions.
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include "strategy_gaussian_conv.h"
#include "filter_util.h"

/** \brief Output file for impulse test */
#define OUTPUT_FILE     "impulse.txt"

#ifndef M_SQRT2PI
/** \brief The constant sqrt(2 pi) */
#define M_SQRT2PI   2.50662827463100050241576528481104525
#endif


/** \brief Print program usage help message */
void print_usage()
{
    puts("Gaussian benchmark, P. Getreuer 2012-2013");
#ifdef NUM_SINGLE
    puts("Configuration: single-precision computation");
#else
    puts("Configuration: double-precision computation");
#endif    
    puts("\nUsage: gaussian_bench [bench type] [options] [output]\n");
    puts("Bench type:");
    puts("   speed         measure computation time");
    puts("   accuracy      measure L^infty operator norm");
    puts("   impulse       impulse response, written to " OUTPUT_FILE "\n");
    puts("Options:");
    puts("   -a <algo>     algorithm to use, choices are");
    puts("                 fir     FIR approximation, tol = kernel accuracy");
    puts("                 dct     DCT-based convolution");
    puts("                 box     box filtering, K = # passes");
    puts("                 sii     stacked integral images, K = # boxes");
    puts("                 am      Alvarez-Mazorra recursive filtering,");
    puts("                         K = # passes, tol = boundary accuracy");
    puts("                 deriche Deriche recursive filtering,");
    puts("                         K = order, tol = boundary accuracy");
    puts("                 vyv     Vliet-Young-Verbeek recursive filtering,");
    puts("                         K = order, tol = boundary accuracy");
    puts("                 pyramid Burt-Adelson Gaussian pyramid");    
    puts("   -s <number>   sigma, standard deviation of the Gaussian");
    puts("   -K <number>   specifies number of steps (box, sii, am)");
    puts("   -t <number>   accuracy tolerance (fir, am, deriche, vyv)");
    puts("   -N <number>   signal length\n");
    puts("   -r <number>   (speed bench) number of runs");
    puts("   -n <number>   (impulse bench) position of the impulse\n");
}

/** \brief struct of program parameters */
typedef struct
{    
    /** \brief Benchmark type */
    const char *bench_type;
    
    /** \brief Length of the input signal */
    long N;
    /** \brief Number of runs for speed measurement */
    long num_runs;
    /** \brief Impulse location (default = N/2) */
    long n0;
    
    /** \brief Name of the convolution algorithm */
    const char *algo;
    /** \brief sigma parameter of the Gaussian */
    double sigma;
    /** \brief Parameter K */
    int K;
    /** \brief Tolerance */
    double tol;
} program_params;

int parse_params(program_params *param, int argc, char **argv);

int speed_test(program_params p, num *output, num *input)
{
    gconv *g = NULL;
    unsigned long time_start, time_stop;
    long run;
    
    if(!(g = gconv_plan(output, input, p.N, 1,
        p.algo, p.sigma, p.K, p.tol)))
        return 0;
    
    time_start = millisecond_timer();
    
    for(run = 0; run < p.num_runs; run++)
        gconv_execute(g);
        
    time_stop = millisecond_timer();
    printf("%.5e\n", 
        ((double)(time_stop - time_start)) / p.num_runs);
    gconv_free(g);
    return 1;
}

void make_impulse_signal(num *signal, long N, long n0)
{
    long n;
    
    for(n = 0; n < N; n++)
        signal[n] = 0;
    
    signal[n0] = 1;
}

int accuracy_test(program_params p, num *output, num *input)
{
    double *error_sums = NULL;
    num *output0 = NULL;
    gconv *g0 = NULL, *g = NULL; 
    double linf_norm = 0.0;
    long m, n;
    int success = 0;
        
    if(!(error_sums = (double *)malloc(sizeof(double) * p.N))
        || !(output0 = (num *)malloc(sizeof(num) * p.N))
        || !(g0 = gconv_plan(output0, input, p.N, 1,
            "fir", p.sigma, p.K, 1e-15))
        || !(g = gconv_plan(output, input, p.N, 1,
            p.algo, p.sigma, p.K, p.tol)))
        goto fail;
    
    for(n = 0; n < p.N; n++)
        error_sums[n] = 0.0;
    
    for(n = 0; n < p.N; n++)
    {
        make_impulse_signal(input, p.N, n);
        gconv_execute(g0);
        gconv_execute(g);
        
        for(m = 0; m < p.N; m++)
            error_sums[m] += fabs((double)output0[m] - (double)output[m]);
    }
    
    for(n = 0; n < p.N; n++)
        if(error_sums[n] > linf_norm)
            linf_norm = error_sums[n];
    
    printf("%.8e\n", linf_norm);
    success = 1;
fail:
    gconv_free(g0);
    gconv_free(g);
    if(output0)
        free(output0);
    if(error_sums)
        free(error_sums);
    return success;
}

int write_output(const char *filename, num *output, num *exact, long N)
{
    FILE *f;
    long n;
    
    if(!(f = fopen(filename, "wt")))
        return 0;
    
    fprintf(f, "# n\toutput value\texact value\n");
    
    for(n = 0; n < N; n++)
        fprintf(f, "%ld\t%.16e\t%.16e\n", n, output[n], exact[n]);
    
    fclose(f);
    return 1;
}

int impulse_test(program_params p, num *output, num *input)
{
    num *exact = NULL;
    gconv *g = NULL, *g_exact = NULL; 
    int success = 0;
    
    if(!(exact = (num *)malloc(sizeof(num) * p.N))
        || !(g = gconv_plan(output, input, p.N, 1,
            p.algo, p.sigma, p.K, p.tol))
        || !(g_exact = gconv_plan(exact, input, p.N, 1,
            "fir", p.sigma, p.K, 1e-15)))
        goto fail;
    
    make_impulse_signal(input, p.N, p.n0);
    gconv_execute(g);
    gconv_execute(g_exact);
    
    if(!write_output(OUTPUT_FILE, output, exact, p.N))
        goto fail;
    
    success = 1;
fail:
    gconv_free(g_exact);
    gconv_free(g);
    if(exact)
        free(exact);
    return success;
}

int main(int argc, char **argv)
{
    program_params param;
    num *input = NULL;
    num *output = NULL;
    int success = 0;
    
    if(!parse_params(&param, argc, argv))
        return 1;
    
    /* Create input signal */
    if(!(input = (num *)malloc(sizeof(num) * param.N))
        || !(output = (num *)malloc(sizeof(num) * param.N)))
    {
        fprintf(stderr, "Out of memory\n");
        goto fail;
    }
    
    if(!strcmp(param.bench_type, "speed"))
    {
        if(!speed_test(param, output, input))
            goto fail;
    }
    else if(!strcmp(param.bench_type, "accuracy"))
    {
        if(!accuracy_test(param, output, input))
            goto fail;
    }
    else if(!strcmp(param.bench_type, "impulse"))
    {
        if(!impulse_test(param, output, input))
            goto fail;
    }
    else
    {
        fprintf(stderr, "Invalid bench type \"%s\"\n", param.bench_type);
        goto fail;
    }
    
    success = 1;
fail:
    if(!success)
        fprintf(stderr, "Bench %s failed\n", param.bench_type);

    if(output)
        free(output);
    if(input)
        free(input);
    return !success;
}

int parse_params(program_params *param, int argc, char **argv)
{
    static const char *default_algo = (const char *)"fir";
    char *option_string;
    char option_char;
    int i;

    if(argc < 2)
    {
        print_usage();
        return 0;
    }

    param->bench_type = argv[1];
    
    /* Set parameter defaults */
    param->N = 100;
    param->n0 = -1;
    param->num_runs = 10;
    param->sigma = 5;
    param->algo = default_algo;
    param->K = 3;
    param->tol = 1e-2;
    
    for(i = 2; i < argc;)
    {
        if(argv[i] && argv[i][0] == '-')
        {
            if((option_char = argv[i][1]) == 0)
            {
                fprintf(stderr, "Invalid parameter format.\n");
                return 0;
            }

            if(argv[i][2])
                option_string = &argv[i][2];
            else if(++i < argc)
                option_string = argv[i];
            else
            {
                fprintf(stderr, "Invalid parameter format.\n");
                return 0;
            }
            
            switch(option_char)
            {
            case 'a':   /* Read algorithm */
                param->algo = option_string;
                break;
            case 's':   /* Read sigma parameter */
                param->sigma = atof(option_string);
                
                if(param->sigma < 0)
                {
                    fprintf(stderr, "sigma must be positive.\n");
                    return 0;
                }
                break;
            case 'K':   /* Read number of steps */
                param->K = atoi(option_string);
                
                if(param->K <= 0)
                {
                    fprintf(stderr, "K must be positive.\n");
                    return 0;
                }
                break;
            case 't':   /* Read tolerance */
                param->tol = atof(option_string);
                
                if(param->tol < 0)
                {
                    fprintf(stderr, "Tolerance must be positive.\n");
                    return 0;
                }
                break;
            case 'N':   /* Input length */
                param->N = atoi(option_string);
                
                if(param->N < 0)
                {
                    fprintf(stderr, "Signal length must be positive.\n");
                    return 0;
                }
                break;
            case 'r':   /* Read number of runs */
                param->num_runs = (long)atof(option_string);
                
                if(param->num_runs <= 0)
                {
                    fprintf(stderr, "Number of runs must be positive.\n");
                    return 0;
                }
                break;
            case 'n':   /* Impulse position */
                param->n0 = atoi(option_string);
                break;
            case '-':
                print_usage();
                return 0;
            default:
                if(isprint(option_char))
                    fprintf(stderr, "Unknown option \"-%c\".\n", option_char);
                else
                    fprintf(stderr, "Unknown option.\n");
            
                return 0;
            }
            
            i++;
        }
        else
            i++;
    }
    
    if(param->n0 < 0 || param->n0 >= param->N)
        param->n0 = (param->N + 1) / 2;
    
    return 1;
}
