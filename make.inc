############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _linux

#
#  The name of the libraries to be created/linked to
#


SuperLUroot	= /home/xxx/SuperLU_4.3
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu_4.3.a


TMGLIB       	= libtmglib.a

## BLASLIB   	= $(SuperLUroot)/lib/libblas.a

BLASLIB		= -L/usr/lib -lblas
LIBS		= $(SUPERLULIB) $(BLASLIB)


## ATLAS BLAS causes single-precision to fail
#BLASLIB   	= -L/usr/lib/atlas -lblas
## This BLAS causes single-precision to fail the test in SuperLU
#BLASLIB 	= -L/usr/lib -lblas

LIBS		= $(SUPERLULIB) $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#

ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = gcc -fPIC
CFLAGS       = -O3 -g -ggdb
NOOPTS       = 
FORTRAN	     = gfortran #g77
FFLAGS       = -O2 -g -fcheck=all -Wall  -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
LOADER       = $(CC)
LOADOPTS     = -fopenmp

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
#MATLAB	     = /usr/sww/matlab

