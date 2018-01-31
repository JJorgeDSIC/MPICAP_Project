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
#define I(i,j) I[i+j*n]
#define T(i,j) T[i+j*n]
#define D(i,j) D[i+j*n]
#define MAJORROW 1

void fillMatrix(double* A, int n){

  int i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      //A(i,j)=rand() % N;
      //A[i+j*N]
      A(i,j)=i;
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


void trasposeMatrix(double* A, int n){

  int i;
  int j;
  double aux;
  printf("\n");
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      aux = A(i,j);
      A(i,j) = A(j,i);
      A(j,i) = aux;
      printf(" %f", A(i,j));
    }
    printf("\n");
  }
  printf("\n");
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
  }}else{
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

void matmult(double *A, double *B, double *C, int n){
  const double one = 1.0;
  dgemm_( "N", "N", &n, &n, &n, &one, A, &n, B, &n, &one, C, &n );
}

void QTQtranspose(double* Q, double *T, double *QTQt, int n){

  const double one = 1.0;
  //double *Qaux;
  //copyMatrix(Qaux, Q, n);
  double* QT=calloc(n*n,sizeof(double));
  dgemm_( "N", "N", &n, &n, &n, &one, Q, &n, T, &n, &one, QT, &n );
  dgemm_( "N", "T", &n, &n, &n, &one, QT, &n, Q, &n, &one, QTQt, &n );
  //printf("====QTQ'====\n");
  //printMat(QTQt, n);

}


void copyMatrix(double* dst, double* src, int n){
  memcpy( (void*)dst, (void*)src, n * n * sizeof(double) );
}

void matsum(double *A, double *B, double *C, int n){
  const double one = 1.0;
  const int oneI = 1;
  int nsquare = n * n;
  memcpy( (void*)C, (void*)B, n * n * sizeof(double) );
  daxpy_(&nsquare, &one, A, &oneI, C, &oneI);
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


  //printf("====DIAG====\n");
  //printVect(WR, n);

  //printf("====DIAG====\n");
  //printVect(WI, n);

  if (INFO!=0){
    printf("Error: dgeev returned error code %d ",INFO);
    return -1;
  }

  return 0;

}


int solverSys(double *A, double *T, double *D, int n){

//Y = zeros(N);
//Daux = zeros(N,1);

double* Y=calloc(n*n,sizeof(double));
double* Daux=calloc(n,sizeof(double));
double* Z=malloc(n*n*sizeof(double));
double* b=malloc(n*sizeof(double));
double* Yaux=malloc(n*sizeof(double));;

double* R=malloc(n*n*sizeof(double));
double* Rs=malloc(n*n*sizeof(double));
double* AP=malloc(n*n*sizeof(double));
double* AR=malloc(n*n*sizeof(double));
double* Dk = malloc(n*sizeof(double));
double* Ds = malloc(n*sizeof(double));

int* IPIV;
int INFO;
int nsquare = n * n;
double alpha;
const double one = 1.0;
const int zero = 0;
const int oneI = 1;
const double minusOne = -1.0;
const int nPlusOne = n+1;


int i = 0;
int j = 0;

while(i < 2){

  printf("======ITERATION: %d======\n", i);
  
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

    //printf("====D====\n");
    //printMat(D, n);
    
    // Dk = D(:,i);
    //dcopy_(&n, &D[i*n], &oneI, Dk, &oneI);
    //printf("====Dk, col=(%d)====\n",i);
    //printVect(Dk, n);
    // P = Dk/T(i+1,i);
    //dscal_(&n, 1.0/T[(i+1)+i*n], Dk, &oneI);
    //printf("====Dk, col=(%d)====\n",i);
    //printVect(Dk, n);
    // Ds = D(:,i+1);
    //dcopy_(&n, &D[(i+1)*n], &oneI, Ds, &oneI);
    //printf("====Ds, col=(%d)====\n",i);
    //printVect(Ds, n);


    //1)Z =(A-T(i,i)*I)
    //Z = A
    //dcopy_(&nsquare, A, &oneI, Z, &oneI);
    // Z = Z-(T(i,i)*I);
    //daxpy_(&n, &minusOne, &T(i,i), &zero, Z, &nPlusOne);

    //2)/T(i+1,i);
    // R = Z/T(i+1,i);
    //dscal_(&nsquare, 1.0/T[(i+1)+i*n], Z, &oneI);

    //A*R
    //dgemm_( "N", "N", &n, &n, &n, &one, A, &n, Z, &n, &one, AR, &n );
    //(A*R-T(i,i+1)*I)
    //daxpy_(&n, &minusOne, &T[i+(i+i)*n], &zero, AR, &nPlusOne);
    //dcopy_(&nsquare, Z, &oneI, Rs, &oneI);
    //R*T(i+1,i+1)
    //dscal_(&nsquare, T[(i+1)+(i+1)*n], Rs, &oneI);

    // J = (A*R-T(i,i+1)*I)-R*T(i+1,i+1);
    // (same as) J = AR-Rs
    //daxpy_(&n, &minusOne, Rs, &one, AR, &one);

    // W = Ds + A*P - P.*T(i+1,i+1);
    //(same as)
    // W = Ds + A*Dk - Dk*T(i+1,i+1);
    //AP=A*Dk
    //dgemm_( "N", "N", &n, &n, &n, &one, A, &n, Dk, &n, &one, AP, &n );
    //Dk*T(i+1,i+1);
    //dscal_(&n, T[(i+1)+(i+1)*n], Dk, &oneI);

    //daxpy_(&n, &one, Ds, &one, AR, &one);
    // %Zx = W
    // %x = J\W

    // Y(:,i) = J\W;
    //dgesv_(&n, &oneI, AR, &n, IPIV, b, &n, &INFO);



    // Y(:,i+1) = R*Y(:,i) - P;
    // i=i+1;
    i++;
    return 0;
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
    for (j = 0; j <i; ++j)
    {
      //Yaux <- Y(:,j) 
      dcopy_(&n, &Y[j*n], &oneI, Yaux, &oneI);
      //alpha <- T(j,i)
      alpha = T[j+i*n];//T(j,i);
      //Yaux*T(j,i)
      dscal_(&n, &alpha, Yaux, &oneI);
      //D(:,i)<-D(:,i) + Yaux*T(j,i)
      daxpy_(&n, &one, Yaux, &oneI, &D[i*n], &oneI);
      
      printf("====Yaux====\n");
      printVect(Yaux, n);
      printf("====D====\n");
      printMat(D, n);
    }
    

    // b = D(:,i); //D[i+j*n]
    dcopy_(&n, &D[i*n], &oneI, b, &oneI);
    //printf("====D====\n");
    //printMat(D, n);
    printf("====b====\n");
    printVect(b, n);
    // %x = Y(:i);
    // %Zx = b
    // %x = Z\b

    // sgesv  ( integer   N,
    //   integer   NRHS,
    //   real, dimension( lda, * )   A,
    //   integer   LDA,
    //   integer, dimension( * )   IPIV,
    //   real, dimension( ldb, * )   B,
    //   integer   LDB,
    //   integer   INFO 
    // ) 
    dgesv_(&n, &oneI, Z, &n, IPIV, b, &n, &INFO);

    //if (INFO!=0){
    //  printf("Error: dgeev returned error code %d ",INFO);
    // return -1;
    //}
    //dcopy_(&n, b, &oneI, &Y[i*n], &oneI);
    //printf("====Y====\n");
    //printMat(Y, n);
    printf("====b====\n");
    printVect(b, n);

    // Y(:,i) = Z\b;

  }
  i++;

}

//I = eye(N);
//i = 1;


return 0;

}


// void matsum(double *A, double *B, double *C, int n){

//   int i,j;
//   for (i=0;i<n;i++){
//     for (j=0;j<n;j++){
//       C(i,j) = A(i,j)+B(i,j);
//     }
//   }

// } 


int main(int argc, char **argv) {

  int i=0,j=0;
  int n = 5;
  double t1,t2,tucpu,tscpu;
  const double one = 1.0;
  int zero = 0;
  double *A;
  double *AX;
  double *B;
  double *BSchur;
  double *BX;
  double *X;
  double *C;
  double *I;

  A=malloc(n*n*sizeof(double));
  B=malloc(n*n*sizeof(double));
  BSchur=malloc(n*n*sizeof(double));
  X=malloc(n*n*sizeof(double));
  C=calloc(n*n,sizeof(double));
  AX=malloc(n*n*sizeof(double));
  BX=malloc(n*n*sizeof(double));
  I=malloc(n*n*sizeof(double));

  for(i = 0; i < n; i++)
    I(i,i) = 1.0;

  srand(time(0));

  ctimer(&t1,&tucpu,&tscpu);
  
  printf("====Filling A====\n");
  fillMatrix(A, n);
  printf("====Filling B====\n");
  fillMatrix(B, n);
  fillMatrix(BSchur, n);
  printf("====Filling X====\n");
  fillMatrix(X, n);
  printf("====Filling C====\n");
  //fillMatrixWithValue(C, 0.0, n);

  // printf("====A====\n");
  // printMat(A, n);
  // printf("====B====\n");
  // printMat(B, n);
  // printf("====X====\n");
  // printMat(X, n);
  // printf("====C====\n");
  // printMat(C, n);
  // printf("====I====\n");
  // printMat(I, n);

  //AX=A*X
  matmult(A,X,AX,n);
  //printf("====AX====\n");
  //printMat(AX, n);
  //BX=B*X
  matmult(B,X,BX,n);

  //printf("====BX====\n");
  //printMat(BX, n);

  //printf("========\n");

  /********************/
  /********************/
  /********************/


  matsum(AX, BX, C, n);

  /********************/
  /********************/
  /********************/

  // printf("\n");
  // printf("====AX====\n");
  // printMat(AX, n);

  // printf("====BX====\n");
  // printMat(BX, n);

  // printf("====C====\n");
  // printMat(C, n);

  // double Btest[] = {
  //   0.8147,0.9058,0.1270,0.9134,0.6324,
  //   0.0975,0.2785,0.5469,0.9575,0.9649,
  //   0.1576,0.9706,0.9572,0.4854,0.8003,
  //   0.1419,0.4218,0.9157,0.7922,0.9595,
  //   0.6557,0.0357,0.8491,0.9340,0.6787};

  double Atest[] = {
     0.7577,0.7060,0.8235,0.4387,0.4898
    ,0.7431,0.0318,0.6948,0.3816,0.4456
    ,0.3922,0.2769,0.3171,0.7655,0.6463
    ,0.6555,0.0462,0.9502,0.7952,0.7094
    ,0.1712,0.0971,0.0344,0.1869,0.7547};


    /**
    0.8147    0.0975    0.1576    0.1419    0.6557
    0.9058    0.2785    0.9706    0.4218    0.0357
    0.1270    0.5469    0.9572    0.9157    0.8491
    0.9134    0.9575    0.4854    0.7922    0.9340
    0.6324    0.9649    0.8003    0.9595    0.6787
    **/

  double Btest[] = {
     0.8147    ,0.9058   ,0.1270   ,0.9134  ,0.6324
    ,0.0975   ,0.2785   ,0.5469  ,0.9575    ,0.9649
    ,0.1576   ,0.9706   ,0.9572  ,0.4854    ,0.8003
    ,0.1419   ,0.4218   ,0.9157  ,0.7922    ,0.9595
    ,0.6557   ,0.0357   ,0.8491  ,0.9340    ,0.6787};

  double Ctest[] = {
    0.2760,0.4984,0.7513,0.9593,0.8407
  ,0.6797,0.9597,0.2551,0.5472,0.2543
  ,0.6551,0.3404,0.5060,0.1386,0.8143
  ,0.1626,0.5853,0.6991,0.1493,0.2435
  ,0.1190,0.2238,0.8909,0.2575,0.9293};

  printf("====Atest====\n");
  printMat(Atest, n);

  printf("====Btest====\n");
  printMat(Btest, n);

  printf("====Ctest====\n");
  printMat(Ctest, n);

  double *T=malloc(n*n*sizeof(double));
  double *Q=malloc(n*n*sizeof(double));

  schurFactorization(Btest, T, Q, n);

  double *D=calloc(n*n, sizeof(double));

  matmult(Ctest,Q,D,n);

  //printf("====D====\n");
  //printMat(D, n);

  // printf("====T====\n");
  // printMat(T, n);

  // printf("====Q====\n");
  // printMat(Q, n);

  double* QTQt=calloc(n*n,sizeof(double));

  QTQtranspose(Q,T,QTQt,n);

  // printf("====Btest====\n");
  // printMat(Btest, n);

  // printf("====QTQ'====\n");
  // printMat(QTQt, n);

  solverSys(Atest, T, D, n);
  // check for errors
 
  // double *WORK=malloc(n*sizeof(double));
  // double *VS=malloc(n*sizeof(double));

  // fillVectorWithValue(WI, 0.0, n);
  // fillVectorWithValue(WR, 0.0, n);
  // fillVectorWithValue(WORK, 0.0, n);
  // fillVectorWithValue(VS, 0.0, n);

  // int n3 = 3*n;

  // int info = -1;

  // dgees_( 'V', 'N', &n, &n, BSchur, &n, &zero, WR, WI, VS, &n, WORK, &n, NULL, &info);

  ctimer(&t2,&tucpu,&tscpu);
  printf("Tiempo %f segundos \n",(float) (t2-t1));
  return 0;
} /* end func main */