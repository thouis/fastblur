A Survey of Gaussian Convolution Algorithms
Pascal Getreuer, getreuer@gmail.com
Version 20130429 (April 29, 2013)

  *** Please cite IPOL article "A Survey of Gaussian Convolution ***
  *** Algorithms" if you publish results obtained with this      ***
  *** software.                                                  ***


== Overview ==

This C source code accompanies with Image Processing On Line (IPOL) article
"A Survey of Gaussian Convolution Algorithms" by Pascal Getreuer.  The 
article and online demo can be found at

    http://www.ipol.im

Future software releases and updates will be posted at 

    http://dev.ipol.im/~getreuer/code/


== License (GPL/BSD) ==

This program is free software: you can redistribute it and/or modify it
under, at your option, the terms of the GNU General Public License as 
published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version, or the terms of the 
simplified BSD license.

You should have received a copy of these licenses along with this program.
If not, see <http://www.gnu.org/licenses/> and
<http://www.opensource.org/licenses/bsd-license.html>.


== Program Usage ==

This source code includes two command line programs: gaussian_demo and
gaussian_bench.

    * gaussian_demo: performs Gaussian convolution on an image

    * gaussian_bench: performs speed, accuracy, and impulse tests
      for 1D Gaussian convolution


--- gaussian_demo ---

Usage: gaussian_conv_demo [options] <input> <output>

where "input" and "output" are BMP files (JPEG, PNG, or TIFF files can also 
be used if the program is compiled with libjpeg, libpng, and/or libtiff).  

Options:
   -a <algo>     algorithm to use, choices are
                 fir     FIR approximation, tol = kernel accuracy
                 dct     DCT-based convolution
                 box     box filtering, K = # passes
                 sii     stacked integral images, K = # boxes
                 am      Alvarez-Mazorra recursive filtering,
                         K = # passes, tol = boundary accuracy
                 deriche Deriche recursive filtering,
                         tol = boundary accuracy
                 vyv     Vliet-Young-Verbeek recursive filtering,
                         K = order, tol = boundary accuracy
   -s <number>   sigma, standard deviation of the Gaussian
   -K <number>   specifies number of steps (box, sii, am)
   -t <number>   accuracy tolerance (fir, am, deriche, yv)


--- gaussian_bench ---

Usage: gaussian_bench [bench type] [options] [output]

Bench type:
   speed         measure computation time
   accuracy      measure L^infty operator norm
   impulse       impulse response, written to impulse.txt

Options:
   -a <algo>     algorithm to use, choices are
                 fir     FIR approximation, tol = kernel accuracy
                 dct     DCT-based convolution
                 box     box filtering, K = # passes
                 sii     stacked integral images, K = # boxes
                 am      Alvarez-Mazorra recursive filtering,
                         K = # passes, tol = boundary accuracy
                 deriche Deriche recursive filtering,
                         K = order, tol = boundary accuracy
                 vyv     Vliet-Young-Verbeek recursive filtering,
                         K = order, tol = boundary accuracy
                 pyramid Burt-Adelson Gaussian pyramid
   -s <number>   sigma, standard deviation of the Gaussian
   -K <number>   specifies number of steps (box, sii, am)
   -t <number>   accuracy tolerance (fir, am, deriche, vyv)
   -N <number>   signal length

   -r <number>   (speed bench) number of runs
   -n <number>   (impulse bench) position of the impulse


Example usage:

    # Perform Gaussian convolution on the image einstein.png
    ./gaussian_demo -s 5 -a deriche -K 3 -t 1e-6 einstein.png blurred.png

    # Compute the L^infty operator error
    ./gaussian_bench accuracy -N 1000 -s 5 -a deriche -K 3 -t 1e-6

Further examples of using these programs are included in the shell scripts 
demo.sh and bench.sh.


== Compiling ==

Instructions are included below for compiling on Linux sytems with GCC, on
Windows with MinGW+MSYS, and on Windows with MSVC.

Compiling requires the FFTW3 Fourier transform library (http://www.fftw.org/).
For supporting additional image formats, the programs can optionally be
compiled with libjpeg, libpng, and/or libtiff.  Windows BMP images are always
supported.


== Compiling (Linux) ==

To compile this software under Linux, first install the development files for
libfftw, libjpeg, libpng, and libtiff.  On Ubuntu and other Debian-based 
systems, enter the following into a terminal:
    sudo apt-get install build-essential libfftw3-dev libjpeg8-dev libpng-dev libtiff-dev
On Redhat, Fedora, and CentOS, use
    sudo yum install make gcc libfftw-devel libjpeg-turbo-devel libpng-devel libtiff-devel

Then to compile the software, use make with makefile.gcc:

    tar -xf gaussian_20130429.tar.gz
    cd gaussian_20130429
    make -f makefile.gcc

This should produce two executables, gaussian_demo and gaussian_bench.

Source documentation can be generated with Doxygen (www.doxygen.org).

    make -f makefile.gcc srcdoc


== Compiling (Windows) ==

The MinGW+MSYS is a convenient toolchain for Linux-like development under 
Windows.  MinGW and MSYS can be obtained from

    http://downloads.sourceforge.net/mingw/

The FFTW3 library is needed to compile the programs.  FFTW3 can be 
obtained from

    http://www.fftw.org/

Instructions for building FFTW3 with MinGW+MSYS can be found at

    http://www.fftw.org/install/windows.html
    http://neuroimaging.scipy.org/doc/manual/html/devel/install/windows_scipy_build.html


--- Building with BMP only ---

The simplest way to build the programs is with support for only BMP images.
In this case, only the FFTW3 library is required.  Edit makefile.gcc and 
comment the LDLIB lines to disable use of libjpeg, libpng, and libtiff:

    #LDLIBJPEG=-ljpeg
    #LDLIBPNG=-lpng -lz
    #LDLIBTIFF=-ltiff

Then open an MSYS terminal and compile the program with 

    make CC=gcc -f makefile.gcc

This should produce the executables gaussian_demo and gaussian_bench.


--- Building with PNG, JPEG, and/or TIFF support ---

To use the programs with PNG, JPEG, and/or TIFF images, the following
libraries are needed.

    For PNG:    libpng and zlib
    For JPEG:   libjpeg 
    For TIFF:   libtiff

These libraries can be obtained at 
    
    http://www.libpng.org/pub/png/libpng.html
    http://www.zlib.net/
    http://www.ijg.org/
    http://www.remotesensing.org/libtiff/

It is not necessary to include support for all of these libraries, for 
example, you may choose to support only PNG by building zlib and libpng 
and commenting the LDLIBJPEG and LDLIBTIF lines in makefile.gcc.

Instructions for how to build the libraries with MinGW+MSYS are provided at

    http://permalink.gmane.org/gmane.comp.graphics.panotools.devel/103
    http://www.gaia-gis.it/spatialite-2.4.0/mingw_how_to.html

Once the libraries are installed, build the programs with the makefile.gcc
included in this archive.

    make CC=gcc -f makefile.gcc

This should produce two executables, gaussian_demo and gaussian_bench.


== Compiling (Mac OSX) ==

The following instructions are untested and may require adaptation, but
hopefully they provide something in the right direction.

First, install the XCode developer tools.  One way to do this is from 
the OSX install disc, in which there is a folder of optional installs 
including a package for XCode.

The program requires the free FFTW library to build.  Optionally, it can
also use the libpng, libjpeg, and libtiff libraries to support more image
formats.  These libraries can be obtained on Mac OSX from the Fink project 

    http://www.finkproject.org/

Go to the Download -> Quick Start page for instructions on how to get 
started with Fink.  The Fink Commander program may then be used to download
and install the packages fftw, libpng, libjpeg, and libtiff.  It may be 
necessary to install fftw-dev, libpng-dev, libjpeg-dev, and libtiff-dev as 
well.

Once the libraries are installed, compile the program using the included
makefile "makefile.gcc":

    make -f makefile.gcc

This should produce the executables gaussian_demo and gaussian_bench.
