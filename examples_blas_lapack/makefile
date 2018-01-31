##
# @file      makefile
# @author    Mitch Richling <https://www.mitchr.me/>
# @Copyright Copyright 1997,2014 by Mitch Richling.  All rights reserved.
# @brief     Simple make file to build the examples in this directory.@EOL
# @Keywords  
# @Std       GenericMake
#
#            
#            
TARMINE = example test_lapack ax_bx_c
TARGBC = blas1C blas2C blas3C blas2bC
TARGLC = slvSysC                 

# Put your compiler and basic compile options here
FC       = gfortran
CC       = gcc
FFLAGS   = -Wall
CFLAGS   = -Wall -g

#CFLAGS   = -Wall -Wunused-variable

# Put your Fortran 77 LAPACK and BLAS compile and link options here!
# What is below works on my stock Debian box, and works on my OS X
# 10.9.2 system when I comment out the first three and uncomment
# APLKRZ.
BLLIB  = -lblas
LPLIB  = -llapack
LPINC  = -I/usr/include/atlas
#APLKRZ  = -framework Accelerate 

# You will need a Fortran compiler for the TARG*F targets.
#TARGETS = $(TARGBC) $(TARGLC) $(TARMINE)
#TARGETS = $(TARGBC) $(TARGLC) 
TARGETS = $(TARMINE)
#
all : $(TARGETS)
	@echo Make Complete

test_lapack: test_lapack.cpp
	g++ -o test_lapack test_lapack.cpp -llapack

example_lapack: example_lapack.c
	(CC) $(CFLAGS) $(APLKRZ) $(LPINC) example_lapack.cpp ctimer.c -lm -lblas -llapack -o example_lapack

example : example.c
	$(CC) $(CFLAGS) $(APLKRZ) $(LPINC) example.c ctimer.c -lm $(LPLIB) $(BLLIB) -o example

ax_bx_c : ax_bx_c.c
	$(CC) $(CFLAGS) $(APLKRZ) $(LPINC) ax_bx_c.c ctimer.c -lm $(LPLIB) $(BLLIB) -o ax_bx_c


blas1C : blas1C.c blaio.c blaio.h
	$(CC) $(CFLAGS) $(APLKRZ) blas1C.c blaio.c -lm $(BLLIB) -o blas1C

blas2bC : blas2bC.c blaio.c blaio.h 
	$(CC) $(CFLAGS) $(APLKRZ) blas2bC.c blaio.c -lm  $(BLLIB) -o blas2bC

blas2C : blas2C.c blaio.c blaio.h 
	$(CC) $(CFLAGS) $(APLKRZ) blas2C.c blaio.c   -lm $(BLLIB) -o blas2C

blas3C : blas3C.c blaio.c blaio.h
	$(CC) $(CFLAGS) $(APLKRZ) blas3C.c blaio.c  -lm $(BLLIB) -o blas3C

slvSysC : slvSysC.c blaio.c blaio.h
	$(CC) $(CFLAGS) $(APLKRZ) $(LPINC) slvSysC.c blaio.c  -lm $(LPLIB) $(BLLIB) -o slvSysC

clean :
	rm -rf a.out *~ *.bak $(TARGBF) $(TARGBC) $(TARGLF) $(TARGLC) $(TARMINE)	