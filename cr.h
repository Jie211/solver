#ifndef CR_H_INCLUDED__
#define CR_H_INCLUDED__

#include "blas.h"
#include "io.h"
#include "start.h"

extern int 
CR_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max);

extern void 
CR_Init(double *v1, 
    double *v2, 
    double *v3, 
    double *v4, 
    double *x, 
    double ndata);
#endif //CR_H_INCLUDED__

