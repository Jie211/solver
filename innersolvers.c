#include "innersolvers.h"

int InnerSolverSelecter(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max, int kskip, int fix){
  int error=0;
#ifdef IS_CG
  error=CG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif IS_CR
  error=CR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif IK_CG
  error=KSKIPCG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
#elif IK_CR
  error=KSKIPCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
#else
  printf("no inner solver selected\n");
  return -1;
#endif
  if(error==-1){
    printf("error in innersolvers\n");
    return -1;
  }
  return 0;
}
