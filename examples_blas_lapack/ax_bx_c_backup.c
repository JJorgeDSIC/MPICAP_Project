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
#define B(i,j)  B[i+j*n]
#define X(i,j)  X[i+j*n]
#define C(i,j)  C[i+j*n]  
#define AX(i,j) AX[i+j*n]
#define BX(i,j) BX[i+j*n]
#define T(i,j) T[i+j*n]
#define D(i,j) D[i+j*n]
#define MAJORROW 1

   void fillMatrix(double* A, int n){

    int i,j;
    for (i=0;i<n;i++){
      for (j=0;j<n;j++){
        A(i,j)=rand() % n;
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

    //ROWS
    if(MAJORROW){
      for (i=0;i<n;i++){
        for (j=0;j<n;j++){
          printf(" %f", A(i,j));
        }
        printf("\n");
      }
    }else{
      //COLS
      for (j=0;j<n;j++){
        for (i=0;i<n;i++){
          printf(" %f", A(i,j));
        }
        printf("\n");
      }
    }
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

  void matmult(double *A, double *B, double *C, int n){
    const double one = 1.0;
    dgemm_( "N", "N", &n, &n, &n, &one, A, &n, B, &n, &one, C, &n );
  }

  void matsum(double *A, double *B, double *C, int n){
    const double one = 1.0;
    const int oneI = 1;
    int nsquare = n * n;
    memcpy( (void*)C, (void*)B, n * n * sizeof(double) );
    daxpy_(&nsquare, &one, A, &oneI, C, &oneI);
  }

  void QTQtranspose(double* Q, double *T, double *QTQt, int n){

    const double one = 1.0;
  //double *Qaux;
    double* QT=calloc(n*n,sizeof(double));
    dgemm_( "N", "N", &n, &n, &n, &one, Q, &n, T, &n, &one, QT, &n );
    dgemm_( "N", "T", &n, &n, &n, &one, QT, &n, Q, &n, &one, QTQt, &n );
  //printf("====QTQ'====\n");
  //printMat(QTQt, n);

  }

  void copyMatrix(double* dst, double* src, int n){
    memcpy( (void*)dst, (void*)src, n * n * sizeof(double) );
  }


  int schurFactorization(double* A, double* T, double* Q, int n){

  char JOBVS ='V';   //Schur vectors are computed.
  char SORT ='N';
  int LWORK = 3*n; 
  double *WORK =  (double *)malloc( LWORK*sizeof(double));
  int ZERO = 0;
  int INFO;
  int SDIM;  
  double *WR=malloc(n*sizeof(double));
  double *WI=malloc(n*sizeof(double));
  int *BWORK=malloc(n*sizeof(int));

  memcpy( (void*)T, (void*)A, n * n * sizeof(double) );

  dgees_( &JOBVS, &SORT, 0, &n, T, &n, &SDIM, WR, WI, Q, &n, WORK, &LWORK, BWORK, &INFO);

  if (INFO!=0){
    printf("Error: dgeev returned error code %d ",INFO);
    return -1;
  }

  return 0;

}


int solverSys(double *A, double *T, double *D, double* Y, int n){

  double* Z=malloc(n*n*sizeof(double));
  double* Yaux=malloc(n*sizeof(double));
  double* b=malloc(n*sizeof(double));


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

    if(i != n-1 & T(i+1,i)!=0){

     // for j=1:i
    //     D(:,i) = D(:,i) + Y(:,j)*T(j,i);
    //     D(:,i+1) = D(:,i+1) + Y(:,j)*T(j,i+1);
    // end
      printf("No es diagonal superior!!\n");
    // disp("No es triangular superior");
      for (j = 0; j < i; ++j)
      {
      /* code */
        dcopy_(&n, &Y[j*n], &oneI, Yaux, &oneI);

      alpha = T[j+i*n];//T(j,i);

      dscal_(&n, &alpha, Yaux, &oneI);
      daxpy_(&n, &one, Yaux, &oneI, &D[i*n], &oneI);


      dcopy_(&n, &Y[j*n], &oneI, Yaux, &oneI);

      alpha = T[j+(i+1)*n];//T(j,i+1);
      dscal_(&n, &alpha, Yaux, &oneI);

      daxpy_(&n, &one, Yaux, &oneI, &D[(i+1)*n], &oneI);
      
    }

    printf("====D====\n");
    printMat(D, n);
    
    // Dk <- D(:,i);
    dcopy_(&n, &D[i*n], &oneI, Dk, &oneI);
    printf("====Dk====\n");
    printVect(Dk, n);
    // Ds <- D(:,i+1);
    dcopy_(&n, &D[(i+1)*n], &oneI, Ds, &oneI);
    printf("====Ds====\n");
    printVect(Ds, n);
    //P<-Dk
    dcopy_(&n, Dk, &oneI, P, &oneI);
    // P = Dk/T(i+1,i);
    alpha = 1.0/T[(i+1)+i*n];
    dscal_(&n, &alpha, P, &oneI);

    printf("====P====\n");
    printVect(P, n);

    // R = (A-T(i,i)*I)/T(i+1,i);
    //(A-T(i,i)*I)
    //R = A
    dcopy_(&nsquare, A, &oneI, R, &oneI);
    // R = R-(T(i,i)*I);
    daxpy_(&n, &minusOne, &T(i,i), &zero, R, &nPlusOne);

    //(A-T(i,i)*I)/T(i+1,i)
    alpha = 1.0/T[(i+1)+i*n];
    dscal_(&nsquare, &alpha, R, &oneI);

    printf("====R====\n");
    printMat(R, n);
    // Z = (A*R-T(i,i+1)*I)-R*T(i+1,i+1);
    //A*R
    dgemm_( "N", "N", &n, &n, &n, &one, A, &n, R, &n, &one, AR, &n );
    //A*R-T(i,i+1)*I
    daxpy_(&n, &minusOne, &T[i+(i+1)*n], &zero, AR, &nPlusOne);
    //Rs <- R
    dcopy_(&nsquare, R, &oneI, Rs, &oneI);
    //Rs*T(i+1,i+1)
    alpha = T[(i+1)+(i+1)*n];
    dscal_(&nsquare, &alpha, Rs, &oneI);
    //Z = AR-Rs;
    daxpy_(&nsquare, &minusOne, Rs, &oneI, AR, &oneI);

    printf("====Z====\n");
    printMat(AR, n);

    // W = Ds + A*P - P.*T(i+1,i+1);
    //AP=A*P
    dgemv_( "N", &n, &n, &one, A, &n, P, &oneI, &one, AP, &oneI );
    printf("====AP====\n");
    printVect(AP, n);
    //P*T(i+1,i+1);
    //Ps<-P
    dcopy_(&n, P, &oneI, Ps, &oneI);
    alpha = T[(i+1)+(i+1)*n];
    dscal_(&n, &alpha, Ps, &oneI);
    printf("====P*T(i+1,i+1)====\n");
    printVect(P, n);
    //Ds<-W
    //Ds<- Ds + A*P
    daxpy_(&n, &one, AP, &oneI, Ds, &oneI);
    //Ds<- Ds - Ps
    daxpy_(&n, &minusOne, Ps, &oneI, Ds, &oneI);

    printf("====W====\n");
    printVect(Ds, n);
    // Y(:,i) = Z\W;
    dgesv_(&n, &oneI, AR, &n, IPIV, Ds, &n, &INFO);

    if (INFO!=0){
      printf("Error: dgeev returned error code %d ",INFO);
      return -1;
    }
    //Y(:,i+1)<-Ds
    dcopy_(&n, Ds, &oneI, &Y[i*n], &oneI);

    printf("====Y(:,i)====\n");
    printVect(&Y[i*n], n);

    // Y(:,i+1) = R*Y(:,i) - P;
    //Ry<-R*Y(:,i)
    printf("====P 2====\n");
    printVect(P, n);
    dgemv_( "N", &n, &n, &one, R, &n, &Y[i*n], &oneI, &one, Ry, &oneI );
    printf("====Ry====\n");
    printVect(Ry, n);
    //Y(:,i+1)<- Ry - P
    //Ry<- Ry - P
    daxpy_(&n, &minusOne, P, &oneI, Ry, &oneI);
    printf("====Ry====\n");
    printVect(Ry, n);
    //Y(:,i+1)<-Ry
    dcopy_(&n, Ry, &oneI, &Y[(i+1)*n], &oneI);

    printf("====Y(:,i+1)====\n");
    printVect(&Y[(i+1)*n], n);

    i++;
    //return 0;
  }else{

    //Z = A
    dcopy_(&nsquare, A, &oneI, Z, &oneI);
    // Z = Z-(T(i,i)*I);
    daxpy_(&n, &minusOne, &T(i,i), &zero, Z, &nPlusOne);
    // printf("====A====\n");
    // printMat(A, n);
    // printf("====T====\n");
    // printMat(T, n);
    printf("====Z====\n");
    printMat(Z, n);

    // for j=1:i
    //     D(:,i) = D(:,i) + Y(:,j)*T(j,i);
    // end
    for (j = 0; j < i; ++j)
    {

      dcopy_(&n, &Y[j*n], &oneI, Yaux, &oneI);
      // printf("====Y====\n");
      // printMat(Y, n);
      // printf("====Yaux====\n");
      // printVect(Yaux, n);
      alpha = T[j+i*n];//T(j,i);
      dscal_(&n, &alpha, Yaux, &oneI);
      //printf("====Yaux====\n");
      //printVect(Yaux, n);
      daxpy_(&n, &one, Yaux, &oneI, &D[i*n], &oneI);
      
      //printf("====Yaux====\n");
      //printVect(Yaux, n);
      printf("====D====\n");
      printMat(D, n);
    }
    

    // b = D(:,i); //D[i+j*n]
    dcopy_(&n, &D[i*n], &oneI, b, &oneI);
    // Y(:,i) = Z\b;
    dgesv_(&n, &oneI, Z, &n, IPIV, b, &n, &INFO);

    if (INFO!=0){
      printf("Error: dgeev returned error code %d ",INFO);
      return -1;
    }

    dcopy_(&n, b, &oneI, &Y[i*n], &oneI);
    
    //printf("====b====\n");
    //printVect(b, n);

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

  double* AX=calloc(n*n,sizeof(double));
  double* BX=calloc(n*n,sizeof(double));

  matmult(A,X,AX,n);
  matmult(B,X,BX,n);

  matsum(AX, BX, C, n);

  free(AX);
  free(BX);

  return 0;
}

int main(int argc, char **argv) {

  int i=0,j=0;
  int n = 5;
  double t1,t2,tucpu,tscpu;
  const double one = 1.0;
  const double minusOne = -1.0;
  const int zero = 0;
  const int oneI = 1;
  const int nPlusOne = n+1;
  int nsquare = n * n;

  srand(time(0));

  ctimer(&t1,&tucpu,&tscpu);


  double* A=malloc(n*n*sizeof(double));
  double* B=malloc(n*n*sizeof(double));
  double* X=malloc(n*n*sizeof(double));
  double* C=malloc(n*n*sizeof(double));

  generateProblem(A,B,X,C,n);
  // double A[] = {
  //    0.7577,0.7060,0.8235,0.4387,0.4898
  //   ,0.7431,0.0318,0.6948,0.3816,0.4456
  //   ,0.3922,0.2769,0.3171,0.7655,0.6463
  //   ,0.6555,0.0462,0.9502,0.7952,0.7094
  //   ,0.1712,0.0971,0.0344,0.1869,0.7547};

  // double B[] = {
  //    0.8147    ,0.9058   ,0.1270   ,0.9134  ,0.6324
  //   ,0.0975   ,0.2785   ,0.5469  ,0.9575    ,0.9649
  //   ,0.1576   ,0.9706   ,0.9572  ,0.4854    ,0.8003
  //   ,0.1419   ,0.4218   ,0.9157  ,0.7922    ,0.9595
  //   ,0.6557   ,0.0357   ,0.8491  ,0.9340    ,0.6787};

  // double C[] = {
  //   0.2760,0.4984,0.7513,0.9593,0.8407
  // ,0.6797,0.9597,0.2551,0.5472,0.2543
  // ,0.6551,0.3404,0.5060,0.1386,0.8143
  // ,0.1626,0.5853,0.6991,0.1493,0.2435
  // ,0.1190,0.2238,0.8909,0.2575,0.9293};

  double* Xcalc=calloc(n*n,sizeof(double));
  double* Y=calloc(n*n,sizeof(double));

  printf("====A====\n");
  printMat(A, n);

  printf("====B====\n");
  printMat(B, n);

  printf("====C====\n");
  printMat(C, n);

  double *T=malloc(n*n*sizeof(double));
  double *Q=malloc(n*n*sizeof(double));

  schurFactorization(B, T, Q, n);

  double *D=calloc(n*n, sizeof(double));

  //matmult(C,Q,D,n);
  dgemm_( "N", "N", &n, &n, &n, &one, C, &n, Q, &n, &one, D, &n );

  printf("====D====\n");
  printMat(D, n);

  // printf("====T====\n");
  // printMat(T, n);

  // printf("====Q====\n");
  // printMat(Q, n);

  //double* QTQt=calloc(n*n,sizeof(double));

  //QTQtranspose(Q,T,QTQt,n);

  // printf("====B====\n");
  // printMat(B, n);

  // printf("====QTQ'====\n");
  // printMat(QTQt, n);

  printf("====Solving the System====\n");
  solverSys(A, T, D, Y, n);

  printf("Iteration ended\n");
  printf("====Y====\n");
  printMat(Y, n);

  //matmult(C,Q,D,n);
  dgemm_( "N", "T", &n, &n, &n, &one, Y, &n, Q, &n, &one, Xcalc, &n );
  
  printf("====X====\n");
  printMat(Xcalc, n);

  double* AX=calloc(n*n,sizeof(double));
  double* BX=calloc(n*n,sizeof(double));

  dgemm_( "N", "N", &n, &n, &n, &one, A, &n, Xcalc, &n, &one, AX, &n );
  dgemm_( "N", "N", &n, &n, &n, &one, Xcalc, &n, B, &n, &one, BX, &n );

  daxpy_(&nsquare, &minusOne, BX, &oneI, AX, &oneI);

  printf("====C====\n");
  printMat(C, n);

  daxpy_(&nsquare, &minusOne, C, &oneI, AX, &oneI);

  printf("====X====\n");
  printMat(X, n);

  printf("====Xcalc====\n");
  printMat(Xcalc, n);

  daxpy_(&nsquare, &minusOne, Xcalc, &oneI, X, &oneI);

  printf("====Tot====\n");
  printf("Frob. Norm with C=%.8f\n",frobeniusNorm(AX, n));
  //printMat(AX, n);
  printf("====Tot====\n");
  printf("Frob. Norm with real X=%.8f\n", frobeniusNorm(X, n));

  

  ctimer(&t2,&tucpu,&tscpu);
  printf("Tiempo %f segundos \n",(float) (t2-t1));
  return 0;
} /* end func main */