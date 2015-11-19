#ifndef VPCG_H_INCLUDED__
#define VPCG_H_INCLUDED__

#include "../functions/io.h"
#include "../functions/blas.h"
#include "../innersolvers.h"
#include "../start.h"


extern void 
VPCG_Init(double *v1, 
    double *v2, 
    double *v3, 
    double *v4, 
    double *x, 
    double ndata);
extern int 
VPCG_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max);

#endif //VPCG_H_INCLUDED__

