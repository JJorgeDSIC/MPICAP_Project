#include <stdlib.h> 
#include <stdio.h>              /* I/O lib         ISOC  */
#include "utils.h"

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

    free(WR);
    free(WI);
    free(BWORK);

    return 0;

}


int solveSystemWithY(double *A, double *T, double *D, double* Y, int n){

  double* Z=malloc(n*n*sizeof(double));
  double* Yaux=malloc(n*sizeof(double));

  //Non-triangular things
  double* Dk = malloc(n*sizeof(double));
  double* Ds = malloc(n*sizeof(double));

  double* R = malloc(n*n*sizeof(double));

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
    
    if(i != n-1 & T[(i+1)+i*n]!=0){

      for (j = 0; j < i; ++j)
      {

        alpha = T[j+i*n];//T(j,i);

        daxpy_(&n, &alpha, &Y[j*n], &oneI, &D[i*n], &oneI);

        alpha = T[j+(i+1)*n];//T(j,i+1);
       
        daxpy_(&n, &alpha, &Y[j*n], &oneI, &D[(i+1)*n], &oneI);
      }
    
    // Dk <- D(:,i);
    dcopy_(&n, &D[i*n], &oneI, Dk, &oneI);
   
    // Ds <- D(:,i+1);
    dcopy_(&n, &D[(i+1)*n], &oneI, Ds, &oneI);

    // Dk = Dk/T(i+1,i);
    alpha = 1.0/T[(i+1)+i*n];
    dscal_(&n, &alpha, Dk, &oneI);

    // R = (A-T(i,i)*I)/T(i+1,i);
    //(A-T(i,i)*I)
    //R = A
    dcopy_(&nsquare, A, &oneI, R, &oneI);

    // R = R-(T(i,i)*I);
    daxpy_(&n, &minusOne, &T[i+i*n], &zero, R, &nPlusOne);
    //(R-T(i,i)*I)/T(i+1,i)
    alpha = 1.0/T[(i+1)+i*n];
    dscal_(&nsquare, &alpha, R, &oneI);

    // Z = (A*R-T(i,i+1)*I)-R*T(i+1,i+1);
    //A*R
    dgemm_( "N", "N", &n, &n, &n, &one, A, &n, R, &n, &zero, AR, &n );
    //A*R-T(i,i+1)*I
    daxpy_(&n, &minusOne, &T[i+(i+1)*n], &zero, AR, &nPlusOne);

    //Z = AR-R*T(i+1,i+1);
    alpha = -T[(i+1)+(i+1)*n];
    daxpy_(&nsquare, &alpha, R, &oneI, AR, &oneI);

    // W = Ds + A*Dk - Dk*T(i+1,i+1);
    //AP=A*P
    dgemv_( "N", &n, &n, &one, A, &n, Dk, &oneI, &zero, AP, &oneI );
    //Ds<- Ds + A*Dk
    daxpy_(&n, &one, AP, &oneI, Ds, &oneI);
    //Ds<- Ds - Dk*T(i+1,i+1);
    alpha = -T[(i+1)+(i+1)*n];
    daxpy_(&n, &alpha, Dk, &oneI, Ds, &oneI);

    // Y(:,i) = Z\W;
    dgesv_(&n, &oneI, AR, &n, IPIV, Ds, &n, &INFO);

    if (INFO!=0){
      printf("Error: dgeev returned error code %d ",INFO);
      return -1;
    }
    //Y(:,i+1)<-Ds
    dcopy_(&n, Ds, &oneI, &Y[i*n], &oneI);

    //  Y(:,i+1) = R*Y(:,i) - Dk;
    //  Y(:,i+1)<-R*Y(:,i)
    dgemv_( "N", &n, &n, &one, R, &n, &Y[i*n], &oneI, &zero, &Y[(i+1)*n], &oneI );

    //Y(:,i+1)<- Y(:,i+1) - Dk
    daxpy_(&n, &minusOne, Dk, &oneI, &Y[(i+1)*n], &oneI);

    i++;
    
  }else{

    //Z = A
    dcopy_(&nsquare, A, &oneI, Z, &oneI);
    // Z = Z-(T(i,i)*I);
    daxpy_(&n, &minusOne, &T[i+i*n], &zero, Z, &nPlusOne);

    for (j = 0; j < i; ++j)
    {
      //D(:,i) = D(:,i) + Y(:,j)*T(j,i);
      alpha = T[j+i*n];//T(j,i);
      daxpy_(&n, &alpha, &Y[j*n], &oneI, &D[i*n], &oneI);
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


  free(Z);
  free(Yaux);
  free(Dk);
  free(Ds);
  free(R);
  free(AR);
  free(AP);
  free(IPIV);

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

  //printf("====Solving the System====\n");
  solveSystemWithY(A, T, D, Y, n);

  //Obtaining X = Y*Q'
  dgemm_( "N", "T", &n, &n, &n, &one, Y, &n, Q, &n, &zero, Xcalc, &n );

  free(T);
  free(Q);
  free(D);

  return 0;
}

double frobeniusNorm(double* A, int n){

  int i;
  int j;
  double acum = 0;
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      acum+=A[i+j*n];
    }
  }
  return acum;
}


int generateProblem(double* A, double* B, double* X, double* C, int n){

  //printf("====Filling A====\n");
  fillMatrix(A, n);
  //printf("====Filling B====\n");
  fillMatrix(B, n);
  //printf("====Filling X====\n");
  fillMatrix(X, n);
  //printf("====Filling C====\n");
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

  daxpy_(&nsquare, &minusOne, XB, &oneI, AX, &oneI);
  dcopy_(&nsquare, AX, &oneI, C, &oneI);

  free(AX);
  free(XB);

  //printf("====Problem created====\n");

  return 0;
}


int checkSolution(double* A, double* B, double* C, double* Xcalc, double* Xreal, int n, float diff){

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

  //printf("====Tot====\n");
  //printf("Frob. Norm with C=%.16f\n",frobeniusNorm(AX, n));

  daxpy_(&nsquare, &minusOne, Xcalc, &oneI, Xreal, &oneI);

  //printf("====Tot====\n");
  //printf("Frob. Norm with real X=%.16f\n", frobeniusNorm(Xreal, n));

  printf("%d %.6f %.16f %.16f\n",n,(float) diff, frobeniusNorm(AX, n), frobeniusNorm(Xreal, n));

  free(AX);
  free(XB);

  return 0;
}