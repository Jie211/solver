#ifndef CG_H_INCLUDED__
#define CG_H_INCLUDED__

#include "../functions/blas.h"
#include "../functions/io.h"
#include "../functions/cudafunc.cuh"
#include "../start.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>


extern int 
CG_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    int nnz,
    double eps, 
    int i_max);

extern void 
CG_Init(double *v1, 
    double *v2, 
    double *v3, 
    double *x, 
    double ndata);

#endif //CG_H_INCLUDED__

