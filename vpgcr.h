#ifndef VPGCR_H_INCLUDED__
#define VPGCR_H_INCLUDED__

#include "io.h"
#include "blas.h"
#include "innersolvers.h"
#include "start.h"

extern void 
VPGCR_Init(double *v1, 
    double *v2, 
    double *v3, 
    double *v4, 
    double **v5, 
    double **v6, 
    double *x, 
    int ndata,
    int restart);


extern int 
VPGCR_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max, 
    int restart);



#endif //VPGCR_H_INCLUDED__

