#ifndef VPGMRES_H_INCLUDED__
#define VPGMRES_H_INCLUDED__

#include "../functions/io.h"
#include "../functions/blas.h"
#include "../innersolvers.h"
#include "../start.h"


extern int 
VPGMRES_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    int nnz, 
    double eps, 
    int i_max, 
    int rs);

#endif //VPGMRES_H_INCLUDED__

