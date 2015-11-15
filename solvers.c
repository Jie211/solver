#include "solvers.h"

int SolverSelecter(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max, int kskip, int fix){
  int error=0;
#if (defined (VP_CG) || defined (VP_CR))  &&  (!defined (IS_CG) && !defined (IS_CR) && !defined (IK_CG) && !defined (IK_CR))
  printf("---- if use VP solver, you need to select some Inner solver else ----\n");
  return -1;
#endif

#ifdef S_CG
  printf("---- CG selected ----\n");
  error=CG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif S_CR
  printf("---- CR selected ----\n");
  error=CR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif K_CG
  printf("---- Kskip-CG selected ----\n");
  error=KSKIPCG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
#elif K_CR
  printf("---- Kskip-CR selected ----\n");
  error=KSKIPCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max, kskip, fix);
#elif VP_CG
  printf("---- VPCG selected ----\n");
  error=VPCG_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#elif VP_CR
  printf("---- VPCR selected ----\n");
  error=VPCR_CRS(val, col, ptr, bvec, xvec, ndata, eps, i_max);
#else
  printf("---- no solver selected ----\n");
#endif
  if(error==1){
    printf("---- good ----\n");
  }else if(error==2){
    printf("---- bad ----\n");
  }else if(error==-1){
    printf("** error in solver.c **\n");
    return -1;
  }
  return 0;
}
