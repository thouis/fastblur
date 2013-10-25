# GCC makefile for Gaussian Convolution
# Pascal Getreuer 
# April 28, 2013

# The FFTW3 library (http://www.fftw.org) is required.
# Set the flags needed for linking with fftw3 and fftw3f.
LDFFTW3=-lfftw3 -lfftw3f

# Set this line to compute in single precision instead of double precision
NUM=-DNUM_SINGLE

# Make settings
SHELL=/bin/sh
CFLAGS=-O3 -g -ansi -pedantic -Wall -Wextra $(NUM)
LDLIBS=-lm $(LDFFTW3)

SOURCES=strategy_gaussian_conv.c \
gaussian_conv_fir.c gaussian_conv_dct.c gaussian_conv_am.c \
gaussian_conv_deriche.c gaussian_conv_vyv.c \
gaussian_conv_box.c gaussian_conv_ebox.c \
gaussian_conv_sii.c gaussian_short_conv.c \
filter_util.c erfc_cody.c inverfc_acklam.c invert_matrix.c basic.c

OBJECTS=$(SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip dist-xz

all: libgaussian.a

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

libgaussian.a: $(OBJECTS)
	ar cru $@ $^
	ranlib $@
