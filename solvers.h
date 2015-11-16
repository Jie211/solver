#ifndef SOLVERS_H_INCLUDED__
#define SOLVERS_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "cg.h"
#include "cr.h"
#include "kskipcg.h"
#include "kskipcr.h"
#include "vpcg.h"
#include "vpcr.h"
#include "vpgcr.h"
#include "start.h"

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
    int fix,
    int restart);

#endif //SOLVERS_H_INCLUDED__

