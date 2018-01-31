#include <iostream>
#include <fstream>
#include <math.h>               /* Math stuff      ISOC  */
#include <stdio.h>              /* I/O lib         ISOC  */
#include <stdlib.h> 
#include "ctimer.h"
using namespace std;

#define A(i,j)  A[i+j*n]

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
extern int dgemm_( char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

void fillMatrixWithValue(double* A, double a, int n);
void fillVectorWithValue(double* A, double a, int n);
void fillMatrix(double* A);

void fillVectorWithValue(double* A, double a, int n){

  int i;
  for (i=0;i<n;i++){
      A[i]=a;
  }
}

void fillMatrix(double* A, int n){

  int i,j;
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      //A(i,j)=rand() % n;
      //A[i+j*n]
      A(i,j)=i;
    }
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
 
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){

      printf(" %f", A(i,j));
    }
    printf("\n");
  }
}


void matmult(double *A, double *B, double *C, int n){
  const double one = 1.0;
  dgemm( "N", "N", &n, &n, &n, &one, A, &n, B, &n, &one, C, &n );
}


int main(int argc, char** argv){

  // check for an argument
  //if (argc<2){
  //  cout << "Usage: " << argv[0] << " " << " filename" << endl;
  //  return -1;
  //}

  int i,j;
  int n = 4;
  double t1,t2,tucpu,tscpu;
  const double one = 1.0;
  const int zero = 0.0;

  double *A;
  double *AX;
  double *B;
  double *BX;
  double *X;
  double *C;

  double *BSchur;

  A=new double[n*n];
  B=new double[n*n];
  X=new double[n*n];
  C=new double[n*n];
  AX=new double[n*n];
  BX=new double[n*n];
  BSchur=new double[n*n];


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
  fillMatrixWithValue(C, 0.0, n);

  printf("====A====\n");
  printMat(A,n);
  printf("====B====\n");
  printMat(B,n);
  printf("====X====\n");
  printMat(X,n);
  printf("====C====\n");
  printMat(C,n);

  //AX=A*X
  matmult(A,X,AX,n);
  printf("====AX====\n");
  printMat(AX,n);
  //BX=B*X
  matmult(B,X,BX,n);

  printf("====BX====\n");
  printMat(BX,n);

  printf("========\n");

  // // allocate data
  // char Nchar='N';
  // double *eigReal=new double[n];
  // double *eigImag=new double[n];
  // double *vl,*vr;
  // int one=1;
  // int lwork=6*n;
  // double *work=new double[lwork];
  // int info;

  // // calculate eigenvalues using the DGEEV subroutine
  // dgeev_(&Nchar,&Nchar,&n,data,&n,eigReal,eigImag,
  //       vl,&one,vr,&one,
  //       work,&lwork,&info);



  // // check for errors
  // if (info!=0){
  //   cout << "Error: dgeev returned error code " << info << endl;
  //   return -1;
  // }

  // // output eigenvalues to stdout
  // cout << "--- Eigenvalues ---" << endl;
  // for (int i=0;i<n;i++){
  //   cout << "( " << eigReal[i] << " , " << eigImag[i] << " )\n";
  // }
  // cout << endl;

  // // deallocate
  // delete [] data;
  // delete [] eigReal;
  // delete [] eigImag;
  // delete [] work;

  delete[] A;
  delete[] B;
  delete[] BSchur;
  delete[] X;
  delete[] C;
  delete[] AX;
  delete[] BX;

  return 0;
} 