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

TARGBC = blas1C blas2C blas3C blas2bC
TARGLC = slvSysC                 

# Put your compiler and basic compile options here
FC       = gfortran
CC       = gcc
FFLAGS   = -Wall
CFLAGS   = -Wall

# Put your Fortran 77 LAPACK and BLAS compile and link options here!
# What is below works on my stock Debian box, and works on my OS X
# 10.9.2 system when I comment out the first three and uncomment
# APLKRZ.
BLLIB  = -lblas
LPLIB  = -llapack
LPINC  = -I/usr/include/atlas
#APLKRZ  = -framework Accelerate 

# You will need a Fortran compiler for the TARG*F targets.
TARGETS = $(TARGBC) $(TARGLC) example
#TARGETS = $(TARGBC) $(TARGLC) 

all : $(TARGETS)
	@echo Make Complete

example : example.c
	$(CC) $(CFLAGS) $(APLKRZ) $(LPINC) example.c ctimer.c -lm $(LPLIB) $(BLLIB) -o example

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
	rm -rf a.out *~ *.bak $(TARGBF) $(TARGBC) $(TARGLF) $(TARGLC)	