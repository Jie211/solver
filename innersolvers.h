#ifndef INNERSOLVERS_H_INCLUDED__
#define INNERSOLVERS_H_INCLUDED__

#include "./functions/blas.h"
#include "./CRS/cg.h"
#include "./CRS/cr.h"
#include "./CRS/kskipcg.h"
#include "./CRS/kskipcr.h"
#include "./CRS/gcr.h"
#include "start.h"

extern int
InnerSolverSelecter(double *val, 
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
#endif //INNERSOLVERS_H_INCLUDED__

