#include "cg.h"

void CG_Init(double *v1, double *v2, double *v3, double *x, double ndata){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,ndata);
  DoubleVecInit(x,0.0,ndata);
}

int CG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max){
  /* int i, j, k, n; */
  int loop;
  
  double *rvec, *pvec, *Av, *x_0, error=0.0;
  double alpha, beta, bnorm, rnorm;
  double rr, rr2;

  bool flag=false;

#ifndef INNER
  double t_error=0.0;
  FILE *p_x, *p_his;
  p_x=FileInit("./output/CG_x.txt", "w");
  p_his=FileInit("./output/CG_his.txt", "w");
#endif

  Av=Double1Malloc(ndata);
  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  //init vectory
  CG_Init(rvec, pvec, Av, xvec, ndata);

  DoubleVecCopy(x_0, xvec, ndata);

  // b 2norm
  bnorm = Double2Norm(bvec, ndata);

  //Ax
  DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);

  //r=b-Ax
  DoubleVecSub(rvec, bvec, Av, ndata);

  //p=r
  DoubleVecCopy(pvec, rvec, ndata);

  /* // r 2norm */
  /* rnorm = Double2Norm(rvec, ndata); */
  // (r,r)
  rr=DoubleDot(rvec, rvec, ndata); 

  for(loop=0;loop<i_max;loop++){
    //rnorm
    rnorm = Double2Norm(rvec, ndata);
    error=rnorm/bnorm;
#ifndef INNER
    printf("%d %.12e\n",loop, error);
    fprintf(p_his,"%d %.12e\n",loop, error);
#endif
    if(error <= eps){
      flag=true;
      break;
    }

    //Ap
    DoubleMVMCSR(Av, val, col, ptr, pvec, ndata);

    //alpha=(r,r)/(p,ap)
    alpha = rr / DoubleDot(pvec, Av, ndata);

    //x=x+alpha*pvec
    DoubleScalarxpy(xvec, alpha, pvec, xvec, ndata);

    //r=r-alpha*Ap
    DoubleScalarxpy(rvec, -alpha, Av, rvec, ndata);

    //(r,r)
    rr2=DoubleDot(rvec, rvec, ndata);

    //beta=(r_new,r_new)/(r,r)
    beta = rr2/rr;

    rr=rr2;

    //p=r+beta*p
    DoubleScalarxpy(pvec, beta, pvec, rvec, ndata);

  }
#ifndef INNER
  FileOutPutVec(p_x, xvec, ndata);
  t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);
#endif
#ifdef INNER
  printf("Inner %d %.12e\n", loop, error);
#endif
  Double1Free(Av);
  Double1Free(rvec);
  Double1Free(pvec);
  Double1Free(x_0);
#ifndef INNER
  FileClose(p_x);
  FileClose(p_his);
#endif
  if(flag){
    return 1;
  }else{
    return 2;
  }
  return 0;
}
