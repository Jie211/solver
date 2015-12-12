#ifndef KSKIPCG_H_INCLUDED__
#define KSKIPCG_H_INCLUDED__

#include "../functions/io.h"
#include "../functions/blas.h"
#include "../start.h"


extern void 
KSKIPCG_Init(double *v1, 
    double *v2, 
    double *v3, 
    double *v4, 
    double *v5, 
    double *v6, 
    double *v7, 
    double *v8, 
    double *v9, 
    int ndata, 
    int kskip);
extern int 
KSKIPCG_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    int nnz,
    double eps, 
    int i_max, 
    int kskip, 
    int fix);
#endif //KSKIPCG_H_INCLUDED__

