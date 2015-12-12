#ifndef BLAS_H_INCLUDED__
#define BLAS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../start.h"

extern void 
Display_Mes(char *mes);

extern void 
Display_Err(char *err);

extern double
gettimeofday_sec();

extern double *
Double1Malloc(int ndata);

extern double **
Double2Malloc(int ndatax, 
    int ndatay);

extern void 
Double1Free(double *ptr);

extern void 
Double2Free(double **ptr, 
    int ndatay);

extern int *
Intger1Malloc(int ndata);

extern int **
Intger2Malloc(int ndatax, 
    int ndatay);

extern void 
Intger1Free(int *ptr);
extern void 
Intger2Free(int **ptr, 
    int ndata);

extern void 
DoubleVecAdd(double *out, 
    double *x, 
    double *y, 
    int ndata);
extern void 
DoubleVecSub(double *out, 
    double *x, 
    double *y, 
    int ndata);
extern void 
DoubleVecMul(double *out, 
    double *x, 
    double *y, 
    int ndata);
extern void 
DoubleScalar(double *out, 
    double a, 
    double *x, 
    int ndata);
extern void 
DoubleScalarxpy(double *out, 
    double a, 
    double *x, 
    double *y, 
    int ndata);
extern double 
DoubleDot(double *x, 
    double *y, 
    int ndata);
extern void 
DoubleMVMCSR(double *out, 
    double *val, 
    int *col, 
    int *ptr, 
    double *vec, 
    int ndata);
extern double 
Double1Norm(double *x, 
    int ndata);
extern double 
Double2Norm(double *x, 
    int ndata);
extern double 
DoubleMaxNorm(double *x, 
    int ndata);
extern double 
error_check_CRS(double *val, 
    const int *col, 
    const int *ptr, 
    double *b, 
    double *x_new ,
    double *x_0, 
    int N);
extern void 
DoubleVecCopy(double *a,
    double *b, 
    int ndata);
extern void 
DoubleVecInit(double *vec, 
    double val, 
    int ndata);
extern void 
Double2VecInit(double **vec, 
    double val, 
    int ndatax, 
    int ndatay);
extern void 
DoubleCalArApKCG(double *Ar, 
    double *Ap, 
    double *val, 
    int *col, 
    int *ptr, 
    double *rvec, 
    double *pvec, 
    int ndata, 
    int kskip);
extern void 
DoubleCalDeltaEtaZetaKCG(double *delta, 
    double *eta, 
    double *zeta, 
    double *Ar, 
    double *Ap, 
    double *rvec, 
    double *pvec, 
    int ndata, 
    int kskip);
extern void 
DoubleCalArApKCR(double *Ar, 
    double *Ap, 
    double *val, 
    int *col, 
    int *ptr, 
    double *rvec, 
    double *pvec, 
    int ndata, 
    int kskip);
extern void 
DoubleCalDeltaEtaZetaKCR(double *delta, 
    double *eta, 
    double *zeta, 
    double *Ar, 
    double *Ap, 
    double *rvec, 
    int ndata, 
    int kskip);


#endif //BLAS_H_INCLUDED__

