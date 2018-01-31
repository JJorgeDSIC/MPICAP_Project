#include <stdio.h>              /* I/O lib         ISOC  */
#include <stdlib.h>             /* Standard Lib    ISOC  */            
#define A(i,j)  A[i+j*n]
#define DEBUG 1

void fillMatrix(double* A, int n){
  int i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      A(i,j)= drand48();
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