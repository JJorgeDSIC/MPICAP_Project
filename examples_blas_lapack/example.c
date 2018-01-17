/* -*- Mode:C; Coding:us-ascii-unix; fill-column:132 -*- */
/* ****************************************************************************************************************************** */
/**
   @file      blas2C.c
   @author    Mitch Richling <https://www.mitchr.me/>
   @Copyright Copyright 1997 by Mitch Richling.  All rights reserved.
   @brief     Demonstrate several cblas (level 1) functions. @EOL
   @Keywords  blas cblas C fortran numerical linear algebra vector matrix gemv ger
   @Std       C89

   This is a simple program intended to illistrate how to make use of #gemv and #ger blas routines (as implimented in the cblas).
              
*/

/* ------------------------------------------------------------------------------------------------------------------------------ */

#include <math.h>               /* Math stuff      ISOC  */
#include <stdio.h>              /* I/O lib         ISOC  */
#include <stdlib.h>             /* Standard Lib    ISOC  */
#include <float.h> 
#include "ctimer.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>    /* The MacOS X blas/lapack stuff */
typedef __CLPK_integer       CLPKinteger;
typedef __CLPK_doublereal    CLPKdoublereal;
#else
#include <clapack.h>      /* C LAPACK         */
typedef int       CLPKinteger;
typedef double    CLPKdoublereal;
#endif
#define N 5
#define A(i,j)  A[i+j*N]
#define B(i,j)  B[i+j*N]
#define X(i,j)  X[i+j*N]
#define C(i,j)  C[i+j*N]  
#define AX(i,j) AX[i+j*N]
#define BX(i,j) BX[i+j*N]

int fillMatrix(double* A){

  int i,j;
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      //A(i,j)=rand() % N;
      //A[i+j*N]
      A(i,j)=i;
    }
  }
}


int fillMatrixWithValue(double* A, double a){

  int i,j;
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      //A(i,j)=rand() % N;
      A(i,j)=a;
    }
  }
}


int printMat(double* A){

  int i;
  int j;
 
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){

      printf(" %f", A(i,j));
    }
    printf("\n");
  }
}

void matmult(double *A, double *B, double *C, int n){
  const double one = 1.0;
  dgemm_( "N", "N", &n, &n, &n, &one, A, &n, B, &n, &one, C, &n );
}


void matsum(double *A, double *B, double *C){

  int i,j;
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      C(i,j) = A(i,j)+B(i,j);
    }
  }

}


int main(int argc, char **argv) {

  int i,j;
  int n = N;
  double t1,t2,tucpu,tscpu;
  const double one = 1.0;
  const double zero = 0.0;
  double *A;
  double *AX;
  double *B;
  double *BSchur;
  double *BX;
  double *X;
  double *C;

  A=malloc(N*N*sizeof(double));
  B=malloc(N*N*sizeof(double));
  BSchur=malloc(N*N*sizeof(double));
  X=malloc(N*N*sizeof(double));
  C=malloc(N*N*sizeof(double));
  AX=malloc(N*N*sizeof(double));
  BX=malloc(N*N*sizeof(double));

  srand(time(0));

  ctimer(&t1,&tucpu,&tscpu);
  
  printf("====Filling A====\n");
  fillMatrix(A);
  printf("====Filling B====\n");
  fillMatrix(B);
  fillMatrix(BSchur);
  printf("====Filling X====\n");
  fillMatrix(X);
  printf("====Filling C====\n");
  fillMatrixWithValue(C, 0.0);

  printf("====A====\n");
  printMat(A);
  printf("====B====\n");
  printMat(B);
  printf("====X====\n");
  printMat(X);
  printf("====C====\n");
  printMat(C);

  //AX=A*X
  matmult(A,X,AX,n);
  printf("====AX====\n");
  printMat(AX);
  //BX=B*X
  matmult(B,X,BX,n);

  printf("====BX====\n");
  printMat(BX);

  printf("========\n");

  /********************/
  /********************/
  /********************/


  matsum(AX, BX, C);

  /********************/
  /********************/
  /********************/

  printf("\n");
  printf("====AX====\n");
  printMat(AX);

  printf("====BX====\n");
  printMat(BX);

  printf("====C====\n");
  printMat(C);


  double *WI=malloc(N*sizeof(double));
  double *WR=malloc(N*sizeof(double));
  double *VS=malloc(N*sizeof(double));

  fillVectorWithValue(C, 0.0);
  fillVectorWithValue(C, 0.0);
  fillVectorWithValue(C, 0.0);

  //dgees( 'V', 'N', NULL, &n, BSchur, &n, &zero, WR, WI, VS, &n);

  ctimer(&t2,&tucpu,&tscpu);
  printf("Tiempo %f segundos \n",(float) (t2-t1));
  return 0;
} /* end func main */