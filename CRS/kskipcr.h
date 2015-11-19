#ifndef KSKIPCR_H_INCLUDED__
#define KSKIPCR_H_INCLUDED__

#include "../functions/blas.h"
#include "../functions/io.h"
#include "../start.h"


extern void 
  KSKIPCR_Init(double *v1, 
      double *v2, 
      double *v3, 
      double *v4, 
      double *v5, 
      double *v6, 
      double *v7, 
      double *v8, 
      double *v9, 
      int ndata, 
      int kskip);
extern int 
KSKIPCR_CRS(double *val, 
    int *col, 
    int *ptr, 
    double *bvec, 
    double *xvec, 
    int ndata, 
    double eps, 
    int i_max, 
    int kskip, 
    int fix);
#endif //KSKIPCR_H_INCLUDED__

