#ifndef INNERSOLVERS_H_INCLUDED__
#define INNERSOLVERS_H_INCLUDED__

#include "blas.h"
#include "cg.h"
#include "cr.h"
#include "kskipcg.h"
#include "kskipcr.h"
#include "gcr.h"

extern int
InnerSolverSelecter(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max, 
    int kskip, 
    int fix);
#endif //INNERSOLVERS_H_INCLUDED__

