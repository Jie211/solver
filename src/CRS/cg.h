#ifndef CG_H_INCLUDED__
#define CG_H_INCLUDED__

#include "../functions/blas.h"
#include "../functions/io.h"
#include "../start.h"

extern int 
CG_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max);

extern void 
CG_Init(double *v1, 
    double *v2, 
    double *v3, 
    double *x, 
    double ndata);

#endif //CG_H_INCLUDED__

