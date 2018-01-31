int schurFactorization(double* A, double* T, double* Q, int n);
int solveSystemWithY(double *A, double *T, double *D, double* Y, int n);
int solveSystem(double* A,double* B,double* C, double* Xcalc, int n);
double frobeniusNorm(double* A, int n);
int generateProblem(double* A, double* B, double* X, double* C, int n);
int checkSolution(double* A, double* B, double* C, double* Xcalc, double* Xreal, int n);