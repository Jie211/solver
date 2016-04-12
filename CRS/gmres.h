#ifndef GMRES_H_INCLUDED__
#define GMRES_H_INCLUDED__

#include "../functions/blas.h"
#include "../functions/io.h"
#include "../functions/cudafunc.cuh"
#include "../start.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>


extern int 
GMRES_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    int nnz, 
    double eps, 
    int i_max, 
    int rs);

extern void 
GMRES_Init(double *rvec, 
    double *Av, 
    double *vvec, 
    double *vmtx, 
    double *evec, 
    double *hmtx, 
    double *yvec, 
    double *wvec, 
    double *cvec, 
    double *svec, 
    double *x0vec, 
    int ndata, 
    int rs);
#endif //GMRES_H_INCLUDED__

