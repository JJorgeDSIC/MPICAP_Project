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


#define A(i,j)  A[i+j*n]
#define DEBUG 1

void fillMatrix(double* A, int n){
  int i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){

      A(i,j)= drand48();
    //A[i+j*N]
    //A(i,j)=i;
    }
  }
}

void fillVectorWithValue(double* A, double a, int n){

  int i;
  for (i=0;i<n;i++){
    A[i]=a;
  }
}

void fillMatrixWithValue(double* A, double a, int n){

  int i,j;
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
    //A(i,j)=rand() % N;
      A(i,j)=a;
    }
  }
}

void printMat(double* A, int n){

  int i;
  int j;

    printf("A=[\n");
    for (i=0;i<n;i++){
      for (j=0;j<n;j++){
        if(j==n-1){
          printf(" %f", A(i,j));
        }else{
          printf(" %f,", A(i,j));
        }
      }
      printf("\n");
    }
    printf("];\n");
 
}

void printVect(double* A, int n){

  int i;

  for (i=0;i<n;i++){
    printf(" %f", A[i]);
  }
  printf("\n");
}

double frobeniusNorm(double* A, int n){

  int i;
  int j;
  double acum = 0;
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      acum+=A(i,j);
    }
  }
  return acum;
}

int schurFactorization(double* A, double* T, double* Q, int n){

    char JOBVS ='V';   //Schur vectors are computed.
    char SORT ='N';
    int LWORK = 3*n; 
    double *WORK =  (double *)malloc( LWORK*sizeof(double));
    int INFO;
    int SDIM;  
    double *WR=malloc(n*sizeof(double));
    double *WI=malloc(n*sizeof(double));
    int *BWORK=malloc(n*sizeof(int));

    const int oneI = 1;
    int nsquare = n*n;

    dcopy_(&nsquare, A, &oneI, T, &oneI);
    dgees_( &JOBVS, &SORT, 0, &n, T, &n, &SDIM, WR, WI, Q, &n, WORK, &LWORK, BWORK, &INFO);

    if (INFO!=0){
      printf("Error: dgeev returned error code %d ",INFO);
      return -1;
    }

    return 0;

}


int solveSystemWithY(double *A, double *T, double *D, double* Y, int n){

  double* Z=malloc(n*n*sizeof(double));
  double* Yaux=malloc(n*sizeof(double));

  //Non-triangular things
  double* Dk = malloc(n*sizeof(double));
  double* Ds = malloc(n*sizeof(double));
  double* P = malloc(n*n*sizeof(double));
  double* Ps = malloc(n*n*sizeof(double));

  double* R = malloc(n*n*sizeof(double));
  double* Rs = malloc(n*n*sizeof(double));
  double* Ry = malloc(n*sizeof(double));

  double* AR = malloc(n*n*sizeof(double));
  double* AP = malloc(n*n*sizeof(double));

  int* IPIV =malloc(n*sizeof(int));
  int INFO;
  int nsquare = n * n;
  double alpha;
  const double one = 1.0;
  const double minusOne = -1.0;
  const int zero = 0;
  const int oneI = 1;
  const int nPlusOne = n+1;


  int i = 0;
  int j = 0;

  while(i < n){
    //printf("Iter %d\n",i);
    if(i != n-1 & T[(i+1)+i*n]!=0){
      //printf("No es triangular\n");

      for (j = 0; j < i; ++j)
      {
        //printf("i=%d, j=%d\n",i,j);
        dcopy_(&n, &Y[j*n], &oneI, Yaux, &oneI);

        alpha = T[j+i*n];//T(j,i);

        dscal_(&n, &alpha, Yaux, &oneI);
        daxpy_(&n, &one, Yaux, &oneI, &D[i*n], &oneI);

        dcopy_(&n, &Y[j*n], &oneI, Yaux, &oneI);

        alpha = T[j+(i+1)*n];//T(j,i+1);
        dscal_(&n, &alpha, Yaux, &oneI);

        daxpy_(&n, &one, Yaux, &oneI, &D[(i+1)*n], &oneI);

        
    }
    
    // Dk <- D(:,i);
    dcopy_(&n, &D[i*n], &oneI, Dk, &oneI);
   

   
    // Ds <- D(:,i+1);
    dcopy_(&n, &D[(i+1)*n], &oneI, Ds, &oneI);


    //P<-Dk
    dcopy_(&n, Dk, &oneI, P, &oneI);
    // P = Dk/T(i+1,i);
    alpha = 1.0/T[(i+1)+i*n];
    dscal_(&n, &alpha, P, &oneI);

    // R = (A-T(i,i)*I)/T(i+1,i);
    //(A-T(i,i)*I)
    //R = A
    dcopy_(&nsquare, A, &oneI, R, &oneI);
    // R = R-(T(i,i)*I);
    daxpy_(&n, &minusOne, &T[i+i*n], &zero, R, &nPlusOne);

    //(A-T(i,i)*I)/T(i+1,i)
    alpha = 1.0/T[(i+1)+i*n];
    dscal_(&nsquare, &alpha, R, &oneI);

    // Z = (A*R-T(i,i+1)*I)-R*T(i+1,i+1);
    //A*R

    dgemm_( "N", "N", &n, &n, &n, &one, A, &n, R, &n, &zero, AR, &n );

    //printf("====AR====\n");
    //printMat(AR, n);
    //A*R-T(i,i+1)*I
    daxpy_(&n, &minusOne, &T[i+(i+1)*n], &zero, AR, &nPlusOne);
    //Rs <- R
    dcopy_(&nsquare, R, &oneI, Rs, &oneI);
    //Rs*T(i+1,i+1)
    alpha = T[(i+1)+(i+1)*n];
    dscal_(&nsquare, &alpha, Rs, &oneI);
    //Z = AR-Rs;
    daxpy_(&nsquare, &minusOne, Rs, &oneI, AR, &oneI);

    // W = Ds + A*P - P.*T(i+1,i+1);
    //AP=A*P
    dgemv_( "N", &n, &n, &one, A, &n, P, &oneI, &zero, AP, &oneI );
   
    //P*T(i+1,i+1);
    //Ps<-P
    dcopy_(&n, P, &oneI, Ps, &oneI);
    alpha = T[(i+1)+(i+1)*n];
    dscal_(&n, &alpha, Ps, &oneI);

    //Ds<-W
    //Ds<- Ds + A*P
    daxpy_(&n, &one, AP, &oneI, Ds, &oneI);
    //Ds<- Ds - Ps
    daxpy_(&n, &minusOne, Ps, &oneI, Ds, &oneI);

  
    // Y(:,i) = Z\W;
    dgesv_(&n, &oneI, AR, &n, IPIV, Ds, &n, &INFO);

    if (INFO!=0){
      printf("Error: dgeev returned error code %d ",INFO);
      return -1;
    }
    //Y(:,i+1)<-Ds
    dcopy_(&n, Ds, &oneI, &Y[i*n], &oneI);



    // Y(:,i+1) = R*Y(:,i) - P;
    //Ry<-R*Y(:,i)
    dgemv_( "N", &n, &n, &one, R, &n, &Y[i*n], &oneI, &zero, &Y[(i+1)*n], &oneI );

    //Y(:,i+1)<- Ry - P
    //Ry<- Ry - P
    daxpy_(&n, &minusOne, P, &oneI, &Y[(i+1)*n], &oneI);

    //Y(:,i+1)<-Ry
    //dcopy_(&n, Ry, &oneI, &Y[(i+1)*n], &oneI);

    i++;
    
  }else{

    //Z = A
    dcopy_(&nsquare, A, &oneI, Z, &oneI);
    // Z = Z-(T(i,i)*I);
    daxpy_(&n, &minusOne, &T[i+i*n], &zero, Z, &nPlusOne);


    for (j = 0; j < i; ++j)
    {
      //D(:,i) = D(:,i) + Y(:,j)*T(j,i);
      dcopy_(&n, &Y[j*n], &oneI, Yaux, &oneI);

      alpha = T[j+i*n];//T(j,i);
      dscal_(&n, &alpha, Yaux, &oneI);

      daxpy_(&n, &one, Yaux, &oneI, &D[i*n], &oneI);

    }
    
    // b = D(:,i); //D[i+j*n]
    dcopy_(&n, &D[i*n], &oneI, &Y[i*n], &oneI);
    // Y(:,i) = Z\b;
    dgesv_(&n, &oneI, Z, &n, IPIV, &Y[i*n], &n, &INFO);

    if (INFO!=0){
      printf("Error: dgeev returned error code %d ",INFO);
      return -1;
    }

  }

  i++;

  }

  return 0;

}

int generateProblem(double* A, double* B, double* X, double* C, int n){

  printf("====Filling A====\n");
  fillMatrix(A, n);
  printf("====Filling B====\n");
  fillMatrix(B, n);
  printf("====Filling X====\n");
  fillMatrix(X, n);
  printf("====Filling C====\n");
  //fillMatrix(&C, n);

  double* AX=calloc(n*n,sizeof(double));
  double* XB=calloc(n*n,sizeof(double));

  const double one = 1.0;
  const double minusOne = -1.0;
  const double zero = 0.0;
  const int oneI = 1;
  int nsquare = n * n;

  dgemm_( "N", "N", &n, &n, &n, &one, A, &n, X, &n, &zero, AX, &n );
  dgemm_( "N", "N", &n, &n, &n, &one, X, &n, B, &n, &zero, XB, &n );


  //printf("====AX====\n");
  //printMat(AX, n);

  //printf("====XB====\n");
  //printMat(XB, n);

  daxpy_(&nsquare, &minusOne, XB, &oneI, AX, &oneI);
  dcopy_(&nsquare, AX, &oneI, C, &oneI);

  //printf("====A====\n");
  //printMat(A, n);

  //printf("====B====\n");
  //printMat(B, n);

  //printf("====C====\n");
  //printMat(C, n);

  //printf("====X====\n");
  //printMat(X, n);

  free(AX);
  free(XB);

  printf("====Problem created====\n");

  return 0;
}

int solveSystem(double* A,double* B,double* C, double* Xcalc, int n){

  const double one = 1.0;
  const double zero = 0.0;

  //Triangular matrix
  double *T=malloc(n*n*sizeof(double));
  //Q orthogonal
  double *Q=malloc(n*n*sizeof(double));

  schurFactorization(B, T, Q, n);

  double* D=calloc(n*n, sizeof(double));

  dgemm_( "N", "N", &n, &n, &n, &one, C, &n, Q, &n, &zero, D, &n );

  //Solving the equivalent system with unknown Y 
  double* Y=calloc(n*n,sizeof(double));

  printf("====Solving the System====\n");
  solveSystemWithY(A, T, D, Y, n);

  //Obtaining X = Y*Q'
  dgemm_( "N", "T", &n, &n, &n, &one, Y, &n, Q, &n, &zero, Xcalc, &n );

  return 0;
}

int checkSolution(double* A, double* B, double* C, double* Xcalc, double* Xreal, int n){

  const double one = 1.0;
  const double minusOne = -1.0;
  const double zero = 0.0;
  const int oneI = 1;
  int nsquare = n*n;

  double* AX=calloc(n*n,sizeof(double));
  double* XB=calloc(n*n,sizeof(double));

  dgemm_( "N", "N", &n, &n, &n, &one, A, &n, Xcalc, &n, &zero, AX, &n );
  dgemm_( "N", "N", &n, &n, &n, &one, Xcalc, &n, B, &n, &zero, XB, &n );

  daxpy_(&nsquare, &minusOne, XB, &oneI, AX, &oneI);

  daxpy_(&nsquare, &minusOne, C, &oneI, AX, &oneI);

  daxpy_(&nsquare, &minusOne, Xcalc, &oneI, Xreal, &oneI);

  printf("====Tot====\n");
  printf("Frob. Norm with C=%.16f\n",frobeniusNorm(AX, n));

  printf("====Tot====\n");
  printf("Frob. Norm with real X=%.8f\n", frobeniusNorm(Xreal, n));

  return 0;
}


int main(int argc, char **argv) {

  int n = 100;
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
  printf("Tiempo %f segundos \n",(float) (t2-t1));
  return 0;
} /* end func main */