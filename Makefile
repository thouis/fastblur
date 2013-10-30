# GCC makefile for Gaussian Convolution
# Pascal Getreuer 
# April 28, 2013

# The FFTW3 library (http://www.fftw.org) is required.
# Set the flags needed for linking with fftw3 and fftw3f.
LDFFTW3=-lfftw3 -lfftw3f

# needed for building the fastblur testing program
OPENCV=/usr/local/opencv-clang

# paths for Ray's machine
#STDCPATH=/Users/thouis/homebrew/lib/llvm-3.3/lib/c++/v1

# Make settings
CFLAGS=-O3 -gdwarf-2 -ansi -pedantic -Wall -Wextra -fcilkplus 
LDLIBS=-lm $(LDFFTW3)

CC = gcc
CXX = g++

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

fastblur.o: fastblur.cpp
	$(CXX) -fcilkplus $(CFLAGS) -I$(OPENCV)/include -c $^ 

fastblur: fastblur.o libgaussian.a
	$(CXX) -o $@ $^ $(LDFLAGS) -L$(OPENCV)/lib -lgaussian -lopencv_highgui -lopencv_imgproc -lopencv_core -lcilkrts 

clean:
	rm *.o
	rm fastblur
