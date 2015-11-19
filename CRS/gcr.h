#ifndef GCR_H_INCLUDED__
#define GCR_H_INCLUDED__

#include "../functions/blas.h"
#include "../functions/io.h"
#include "../start.h"

extern void 
GCR_Init(double *v1, 
    double *v2, 
    double *v3, 
    double **v4, 
    double **v5, 
    double *x, 
    int ndata, 
    int restart);

extern int 
GCR_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max, 
    int restart);

#endif //GCR_H_INCLUDED__

