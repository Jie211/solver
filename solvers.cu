#include "solvers.h"

int SolverSelecter(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, int nnz, double eps, int i_max, int kskip, int fix){

  int error=0;

  if(S_CG){
    error=CG_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max);
  }else if(S_CR){
    error=CR_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max);
  }else if(K_CG){
    error=KSKIPCG_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max, kskip, fix);
  }else if(K_CR){
    error=KSKIPCR_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max, kskip, fix);
  }else if(VP_CG){
    error=VPCG_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max);
  }else if(VP_CR){
    error=VPCR_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max);
  }else if(VP_GCR){
    error=VPGCR_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max, restart_outer);
  }else if(S_GCR){
    error=GCR_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max, restart_outer);
  }else if(S_GMRES){
    error=GMRES_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max, restart_outer);
  }else if(VP_GMRES){
    error=VPGMRES_CRS(val, col, ptr, bvec, xvec, ndata, nnz, eps, i_max, restart_outer);
  }else{
    error=-1;
  }

  if(error==1){
    Display_Mes("good");
  }else if(error==2){
    Display_Mes("bad");
  }else if(error==-1){
    Display_Err("error in solver.c");
    return -1;
  }
  return 0;
}
