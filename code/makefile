PROG = program

CC       = gcc
ICC		 = icc
CFLAGS   = -Wall 
INCL=ctimer.c utils.c algebra.c
BLLIB  = -lblas
LPLIB  = -llapack
#MKL    = -O3 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
MKL    = -mkl

TARGETS = $(PROG)
#
all : $(TARGETS)
	@echo Make Complete

program : program.c
	$(CC) $(CFLAGS) program.c $(INCL) -lm $(LPLIB) $(BLLIB) -o program

optimized: program.c
	$(ICC) $(CFLAGS) program.c $(INCL) -O3 $(LPLIB) $(BLLIB) -lm $(MKL) -o program_opt_icc

optimized_gcc: program.c
	$(CC) $(CFLAGS) program.c $(INCL) $(MKL) $(LPLIB) $(BLLIB) -o program_opt_gcc

clean :
	rm -rf a.out *~ *.bak $(PROG)	