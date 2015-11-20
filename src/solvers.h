#ifndef SOLVERS_H_INCLUDED__
#define SOLVERS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "start.h"
#include "./functions/blas.h"
#include "./CRS/cg.h"
#include "./CRS/cr.h"
#include "./CRS/kskipcg.h"
#include "./CRS/kskipcr.h"
#include "./CRS/vpcg.h"
#include "./CRS/vpcr.h"
#include "./CRS/vpgcr.h"
#include "./CRS/gcr.h"

extern int 
SolverSelecter(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max,
    int kskip,
    int fix);

#endif //SOLVERS_H_INCLUDED__

