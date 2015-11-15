#include "vpcr.h"

void VPCR_Init(double *v1, double *v2, double *v3, double *v4, double *v5, double *v6, double *x, double ndata){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,ndata);
  DoubleVecInit(v4,0.0,ndata);
  DoubleVecInit(v5,0.0,ndata);
  DoubleVecInit(v6,0.0,ndata);
  DoubleVecInit(x,0.0,ndata);
}

int VPCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max){
  /* int i, j, k, n; */
  int loop;
  
  double *rvec, *pvec, *qvec, *svec, *zvec, *wvec, *x_0, error=0.0;
  double alpha, beta, bnorm, rnorm;
  double zs, zs2;

  bool flag=false;
  double t_error=0.0;
  double error_message;

  FILE *p_x, *p_his;
  p_x=FileInit("./output/VPCR_x.txt", "w");
  p_his=FileInit("./output/VPCR_his.txt", "w");

  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  qvec=Double1Malloc(ndata);
  svec=Double1Malloc(ndata);
  zvec=Double1Malloc(ndata);
  wvec=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  //init vectory
  VPCR_Init(rvec, pvec, qvec, svec, zvec, wvec, xvec, ndata);

  DoubleVecCopy(x_0, xvec, ndata);

  // b 2norm
  bnorm = Double2Norm(bvec, ndata);

  //Ax
  DoubleMVMCSR(qvec, val, col, ptr, xvec, ndata);

  //r=b-Ax
  DoubleVecSub(rvec, bvec, qvec, ndata);


  //solve z by Az=r
  error_message=InnerSolverSelecter(val, col, ptr, rvec, zvec, ndata, I_EPS, I_I_MAX, I_KSKIP, I_FIX);
  if(error_message!=0){
    printf("error in vpcr\n");
    return -1;
  }

  //p=z
  DoubleVecCopy(pvec, zvec, ndata);

  //q=Ap
  DoubleMVMCSR(qvec, val, col, ptr, pvec, ndata);

  //s=q
  DoubleVecCopy(svec, qvec, ndata);

  // (z,s)
  zs=DoubleDot(zvec, svec, ndata); 

  for(loop=0;loop<i_max;loop++){
    //rnorm
    rnorm = Double2Norm(rvec, ndata);
    error=rnorm/bnorm;
    printf("Outer %d %.12e\n",loop, error);
    fprintf(p_his,"%d %.12e\n",loop, error);
    if(error <= eps){
      flag=true;
      break;
    }

    //init wvec
    DoubleVecInit(wvec, 0.0, ndata);
    
    //solve w by Aw=q
    error_message=InnerSolverSelecter(val, col, ptr, qvec, wvec, ndata, I_EPS, I_I_MAX, I_KSKIP, I_FIX);
    if(error_message!=0){
      printf("error in vpcr\n");
      return -1;
    }

    //alpha=(z,s)/(w,q)
    alpha = zs / DoubleDot(wvec, qvec, ndata);

    //x=x+alpha*pvec
    DoubleScalarxpy(xvec, alpha, pvec, xvec, ndata);

    //r=r-alpha*qvec
    DoubleScalarxpy(rvec, -alpha, qvec, rvec, ndata);

    //z=z-alpha*wvec
    DoubleScalarxpy(zvec, -alpha, wvec, zvec, ndata);

    //s=Az
    DoubleMVMCSR(svec, val, col, ptr, zvec, ndata);

    //(z,s)
    zs2=DoubleDot(zvec, svec, ndata);

    //beta=(z_new,s_new)/(z,s)
    beta = zs2/zs;

    zs=zs2;

    //p=z+beta*p
    DoubleScalarxpy(pvec, beta, pvec, zvec, ndata);

    //q=s+beta*q
    DoubleScalarxpy(qvec, beta, qvec, svec, ndata);

  }
  FileOutPutVec(p_x, xvec, ndata);
  t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);

  Double1Free(rvec);
  Double1Free(pvec);
  Double1Free(qvec);
  Double1Free(zvec);
  Double1Free(wvec);
  Double1Free(x_0);
  FileClose(p_x);
  FileClose(p_his);

  if(flag){
    return 1;
  }else{
    return 2;
  }
  return 0;
}
