PROG = program

CC       = gcc
FFLAGS   = -Wall
CFLAGS   = -Wall
INCL=ctimer.c utils.c algebra.c
BLLIB  = -lblas
LPLIB  = -llapack
LPINC  = -I/usr/include/atlas

TARGETS = $(PROG)
#
all : $(TARGETS)
	@echo Make Complete

program : program.c
	$(CC) $(CFLAGS) $(APLKRZ) $(LPINC) program.c $(INCL) -lm $(LPLIB) $(BLLIB) -o program

clean :
	rm -rf a.out *~ *.bak $(PROG)	