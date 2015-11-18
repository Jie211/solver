#include "solvers.h"

int SolverSelecter(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max, int kskip, int fix){
  int error=0;

#ifdef S_CG
  error=CG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif S_CR
  error=CR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif K_CG
  error=KSKIPCG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
#elif K_CR
  error=KSKIPCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
#elif VP_CG
  error=VPCG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif VP_CR
  error=VPCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif VP_GCR
  error=VPGCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, RESTART);
#elif S_GCR
  error=GCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, RESTART);
#endif

  if(error==1){
    Display_Mes("good");
  }else if(error==2){
    Display_Mes("bad");
  }else if(error==-1){
    Display_Err("erro in solver.c");
    return -1;
  }
  return 0;
}
