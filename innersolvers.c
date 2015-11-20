#include "innersolvers.h"

int InnerSolverSelecter(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max, int kskip, int fix){
  int error=0;
if(IS_CG){
  error=CG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
}else if (IS_CR){
  error=CR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
}else if (IK_CG){
  error=KSKIPCG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
}else if (IK_CR){
  error=KSKIPCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
}else if (IS_GCR){
  error=GCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, restart_inner);
}
  if(error==-1){
    Display_Err("error in innersolvers");
    return -1;
  }
  return 0;
}
