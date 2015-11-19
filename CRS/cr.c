#include "cr.h"

void CR_Init(double *v1, double *v2, double *v3, double *v4, double *x, double ndata){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,ndata);
  DoubleVecInit(v4,0.0,ndata);
  DoubleVecInit(x,0.0,ndata);
}

int CR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max){
  /* int i, j, k, n; */
  int loop;
  
  double *rvec, *pvec, *qvec, *svec, *x_0, error=0.0;
  double alpha, beta, bnorm, rnorm;
  double rs, rs2;

  bool flag=false;

#ifndef INNER
  double t_error=0.0;
  FILE *p_x, *p_his;
  p_x=FileInit("./output/CR_x.txt", "w");
  p_his=FileInit("./output/CR_his.txt", "w");
#endif

  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  qvec=Double1Malloc(ndata);
  svec=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  //init vectory
  CR_Init(rvec, pvec, qvec, svec, xvec, ndata);

  DoubleVecCopy(x_0, xvec, ndata);

  // b 2norm
  bnorm = Double2Norm(bvec, ndata);

  //Ax
  DoubleMVMCSR(qvec, val, col, ptr, xvec, ndata);

  //r=b-Ax
  DoubleVecSub(rvec, bvec, qvec, ndata);

  //p=r
  DoubleVecCopy(pvec, rvec, ndata);

  //q=Ap
  DoubleMVMCSR(qvec, val, col, ptr, pvec, ndata);

  //s=q
  DoubleVecCopy(svec, qvec, ndata);

  // (r,s)
  rs=DoubleDot(rvec, svec, ndata); 

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


    //alpha=(r,s)/(q,q)
    alpha = rs / DoubleDot(qvec, qvec, ndata);

    //x=x+alpha*pvec
    DoubleScalarxpy(xvec, alpha, pvec, xvec, ndata);

    //r=r-alpha*qvec
    DoubleScalarxpy(rvec, -alpha, qvec, rvec, ndata);

    //s=Ar
    DoubleMVMCSR(svec, val, col, ptr, rvec, ndata);

    //(r,s)
    rs2=DoubleDot(rvec, svec, ndata);

    //beta=(r_new,s_new)/(r,s)
    beta = rs2/rs;

    rs=rs2;

    //p=r+beta*p
    DoubleScalarxpy(pvec, beta, pvec, rvec, ndata);

    //q=s+beta*q
    DoubleScalarxpy(qvec, beta, qvec, svec, ndata);

  }
#ifndef INNER
  FileOutPutVec(p_x, xvec, ndata);
  t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);
#endif
#ifdef INNER
  printf("Inner %d %.12e\n", loop, error);
#endif

  Double1Free(rvec);
  Double1Free(pvec);
  Double1Free(qvec);
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
