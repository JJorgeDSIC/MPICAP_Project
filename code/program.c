/* -*- Mode:C; Coding:us-ascii-unix; fill-column:132 -*- */
/* ****************************************************************************************************************************** */
/**
   @file      program.c
   @author    Javier Jorge Cano <http://jjorge.es/>
   @Copyright Copyright 2018 by Javier Jorge Cano.  All rights reserved.
   @brief     Solving AX-XB=C with Schur Factorization. @EOL
   @Keywords  blas cblas C fortran numerical linear algebra vector matrix gemv ger
   @Std       C89

   Solving the equation: AX-XB=C, decomposing B as B=QTQ', that is the Schur factorization,
   results in the following problem that is simpler: AY-YT=D, where Y = XQ and D = CQ, from:

          AX-XB=C
          (B=QTQ')
          AX-XQTQ'=C
          (AX-XQTQ')Q=CQ
          AXQ-XQTQ'Q=CQ
          (as Q is orthogonal, QQ'=I)
          AXQ-XQT=CQ
          (Replacing)
          AY-YT=D

  Considering that T is triangular, it could be solved by substitution and solving 
  simpler Ax=b equations iteratively.
*/

/* ------------------------------------------------------------------------------------------------------------------------------ */

#include <math.h>               /* Math stuff      ISOC  */
#include <stdio.h>              /* I/O lib         ISOC  */
#include <stdlib.h>             /* Standard Lib    ISOC  */            
#include <float.h> 
#include "ctimer.h"
#include "utils.h"
#include "algebra.h"

int main(int argc, char **argv) {

  if ( argc < 2 ) {
    fprintf( stderr, "\nUso: %s <size> \n\n", argv[0] );
    return 0;
  }

  int n;
  sscanf(argv[1],"%d",&n);

  double t1,t2,tucpu,tscpu;

  srand(time(0));

  ctimer(&t1,&tucpu,&tscpu);

  double* A=malloc(n*n*sizeof(double));
  double* B=malloc(n*n*sizeof(double));
  double* X=malloc(n*n*sizeof(double));
  double* C=malloc(n*n*sizeof(double));

  generateProblem(A,B,X,C,n);

  double* Xcalc=calloc(n*n,sizeof(double));

  solveSystem(A,B,C,Xcalc,n);
  
  printf("Iteration ended\n");

  checkSolution(A,B,C,Xcalc,X,n);

  ctimer(&t2,&tucpu,&tscpu);
  printf("%f segundos \n",(float) (t2-t1));

  free(A);
  free(B);
  free(X);
  free(C);

  free(Xcalc);
  
  return 0;
} /* end func main */