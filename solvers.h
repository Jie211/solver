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

