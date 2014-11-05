# ***************************************************************
# Copyright (C) 2008 Universidade do Minho
# All Rights Reserved
# ***************************************************************

# ***************************************************************
#
# User may want to change this settings to match their system
#
# ***************************************************************

# Compilers used
# The C compiler
CC = gcc
# The MPICC compiler
MPICC = /usr/local/openmpi/gcc4/bin/mpicc
# The fortran compiler
FC = g77
# The R compiler
RC = R CMD SHLIB

# Compiler flags
# R flags
#RFLAGS = -I/usr/local/lib/R/include
RFLAGS = -I/usr/share/R/include
# Python flags
#PYFLAGS = -I/usr/local/include/python2.5 -I/usr/local/lib/python2.5/site-packages/numpy/core/include
PYFLAGS = -I/usr/include/python2.5 -I/usr/lib/python2.5/site-packages/numpy/core/include
# C flags
CFLAGS = -O3 -fPIC # -g -Wall # -pg
# ar flags
ARFLAGS = ruv
# Fortran flags
FFLAGS = -04 -u -xl -fno-automatic -fugly-assumed -malign-double

# Additional flags that control compilation
# Linear constraints flag (includes calls to lapack functions)
LFLAGS = -DLINEAR
# Flags for MPI compilation
MPIOPTFLAGS= -DMPI  # -DMPE # Uncomment for MPE data production
# Flags for AMPL compilation
AMPLOPTFLAGS= -DAMPL -I./include


# ar and ranlib binary
AR = ar
RANLIB = ranlib

#
# External libraries: 
#
# Standard libraries. Math and Dynamic Library Load.
# Serial version
libs = -lm -ldl
# Parallel version
mpilibs= -lm -ldl # -lmpe # Uncommnet for MPE use
# AMPL interface library
ampllibs = ./libs/amplsolver.a
# BLAS and LaPACK libs
# Include the necessary libraries to link C and Fortran code
#lapacklibs = ./libs/lapack_LINUX.a ./libs/blas_LINUX.a -L /usr/lib/gcc/x86_64-redhat-linux/3.4.3 -lg2c
lapacklibs = -L/usr/lib -lblas -llapack


# ***************************************************************
#
# End of user settings
#
# ***************************************************************


# Source files
SRC = ./pattern.c ./pswarm.c ./pswarm_main.c ./cache.c ./user.c
HDRS = ./pattern.h ./pswarm.h ./pswarm_main.h
# The MVE interior point code for the Ellipsoid
SRC_LINEAR = ./mve_solver.c ./mve_presolve.c
# Main algorithm files
LIBSRC = ./pattern.c ./pswarm.c               
LIBHDRS = ./pattern.h ./pswarm.h
# R interface
RSRC = ./pswarm_r.c                           
RHDRS = ./pswarm_r.h
# Python interface
PYSRC = ./pswarm_py.c                         
PYHDRS = ./pswarm_py.h



# Default compilation for the parallel version with linear constraints and the PSwarm library
all: parallel_linear lib_linear

# Parallel version with linear constraints
parallel_linear: $(SRC) $(SRC_LINEAR) $(HDRS)
	$(MPICC) $(MPIOPTFLAGS) $(CFLAGS) $(LFLAGS) -o pswarm $(SRC) $(SRC_LINEAR) $(mpilibs) $(lapacklibs)

# Parallel version without linear constraints
parallel: $(SRC) $(HDRS)
	$(MPICC) $(MPIOPTFLAGS) $(CFLAGS) -o pswarm $(SRC) $(mpilibs)

# Serial standalone version with linear constraints
serial_linear: $(SRC) $(SRC_LINEAR) $(HDRS)
	$(CC) $(OPTFLAGS) $(CFLAGS) $(LFLAGS) -o pswarm $(SRC) $(SRC_LINEAR) $(libs) $(lapacklibs)

# Serial standalone version without linear constraints
serial: $(SRC) $(HDRS)
	$(CC) $(OPTFLAGS) $(CFLAGS) -o pswarm $(SRC) $(libs)

# Serial AMPL interface version with linear constraints
ampl_linear: $(SRC) $(SRC_LINEAR) $(HDRS)
	$(CC) $(AMPLOPTFLAGS) $(CFLAGS) $(LFLAGS) -o pswarm $(SRC) $(SRC_LINEAR) $(ampllibs) $(lapacklibs) $(libs)

# Serial AMPL interface version without linear constraints
ampl: $(SRC) $(HDRS)
	$(CC) $(AMPLOPTFLAGS) $(CFLAGS) -o pswarm $(SRC) $(ampllibs) $(libs)

# Serial R interface version with linear constraints
r_linear: lib_linear $(RSRC) $(RHDRS)
	$(CC) $(CFLAGS) $(RFLAGS) -c $(RSRC)
	x=`echo $(RSRC) | sed 's/\.[cs]/.o/g'` && $(RC) libpswarm.a $(lapacklibs) $$x -o pswarm_r.so

# Serial R interface version without linear constraints
r: lib $(RSRC) $(RHDRS)
	$(CC) $(CFLAGS) $(RFLAGS) -c $(RSRC)
	x=`echo $(RSRC) | sed 's/\.[cs]/.o/g'` && $(RC) libpswarm.a $$x -o pswarm_r.so

# Serial Python interface version with linear constraints
py_linear: lib_linear $(PYSRC) $(PYHDRS)
	$(CC) $(CFLAGS) $(PYFLAGS) -c $(PYSRC)
	x=`echo $(PYSRC) | sed 's/\.[cs]/.o/g'` && $(CC) -shared $(CFLAGS) $(PYFLAGS) $$x libpswarm.a $(lapacklibs) -o pswarm_py.so

# Serial Python interface version without linear constraints
py: lib $(PYSRC) $(PYHDRS)
	$(CC) $(CFLAGS) $(PYFLAGS) -c $(PYSRC)
	x=`echo $(PYSRC) | sed 's/\.[cs]/.o/g'` && $(CC) -shared $(CFLAGS) $(PYFLAGS) $$x libpswarm.a -o pswarm_py.so

# PSwarm library with linear constraints
lib_linear: $(LIBSRC) $(LIBHDRS) $(SRC_LINEAR)
	$(CC) $(CFLAGS) $(LFLAGS) -c $(LIBSRC) $(SRC_LINEAR)
	x=`echo $(LIBSRC) $(SRC_LINEAR) | sed 's/\.[cs]/.o/g'` && $(AR) $(ARFLAGS) libpswarm.a $$x
	$(RANLIB) libpswarm.a || true

# PSwarm library without linear constraints
lib: $(LIBSRC) $(LIBHDRS)
	$(CC) $(CFLAGS) -c $(LIBSRC)
	x=`echo $(LIBSRC) | sed 's/\.[cs]/.o/g'` && $(AR) $(ARFLAGS) libpswarm.a $$x
	$(RANLIB) libpswarm.a || true

# Compile the MVE interior point code. For debug purpose only. Uncomment the main function
# in the mve_presolve.c file
mve: mve_presolve.c
	$(CC) -g -o mve_presolve mve_presolve.c $(lapacklibs)

# Just clean up the mess
clean:
	rm -f *.o pswarm libpswarm.a pswarm_r.so pswarm_py.so
