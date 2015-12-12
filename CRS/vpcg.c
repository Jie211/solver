#include "vpcg.h"
void VPCG_Init(double *v1, double *v2, double *v3, double *v4, double *x, double ndata){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,ndata);
  DoubleVecInit(v4,0.0,ndata);
  DoubleVecInit(x,0.0,ndata);
}

int VPCG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, int nnz, double eps, int i_max){
  /* int i, j, k, n; */
  int loop;
  
  double *rvec, *pvec, *zvec, *Av, *x_0, error=0.0;
  double alpha, beta, bnorm, rnorm;
  double rz, rz2;

  bool flag=false;
  double t_error=0.0;
  int error_message;

  double st, et, t1;

  FILE *p_x, *p_his;
  p_x=FileInit("./output/VPCG_x.txt", "w");
  p_his=FileInit("./output/VPCG_his.txt", "w");

  st=gettimeofday_sec();

  Av=Double1Malloc(ndata);
  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  zvec=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  //init vectory
  VPCG_Init(rvec, pvec, zvec, Av, xvec, ndata);

  DoubleVecCopy(x_0, xvec, ndata);

  // b 2norm
  bnorm = Double2Norm(bvec, ndata);

  //Ax
  DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);

  //r=b-Ax
  DoubleVecSub(rvec, bvec, Av, ndata);

  //solve z  by  Az=r
  error_message=InnerSolverSelecter(val, col, ptr, rvec, zvec, ndata, nnz, eps_inner, loop_inner, kskip_inner, fix_inner);
  if(error_message!=0){
    printf("error in vpcg.c\n");
    return -1;
  }

  //p=z
  DoubleVecCopy(pvec, zvec, ndata);

  // (r,z)
  rz=DoubleDot(rvec, zvec, ndata); 

  for(loop=0;loop<i_max;loop++){
    //rnorm
    rnorm = Double2Norm(rvec, ndata);
    error=rnorm/bnorm;
    if(verbose){
      printf("Outer %d %.12e\n",loop, error);
    }
    fprintf(p_his,"%d %.12e\n",loop, error);
    if(error <= eps){
      flag=true;
      break;
    }

    //Ap
    DoubleMVMCSR(Av, val, col, ptr, pvec, ndata);

    //alpha=(r,z)/(p,ap)
    alpha = rz / DoubleDot(pvec, Av, ndata);

    //x=x+alpha*pvec
    DoubleScalarxpy(xvec, alpha, pvec, xvec, ndata);

    //r=r-alpha*Ap
    DoubleScalarxpy(rvec, -alpha, Av, rvec, ndata);
    
    //init zvec
    DoubleVecInit(zvec,0.0,ndata);
    
    //solve new_z Az=r
    error_message=InnerSolverSelecter(val, col, ptr, rvec, zvec, ndata, nnz, eps_inner, loop_inner, kskip_inner, fix_inner);
    if(error_message!=0){
      return -1;
    }

    //(r,z)
    rz2=DoubleDot(rvec, zvec, ndata);

    //beta=(r_new,r_new)/(r,r)
    beta = rz2/rz;

    rz=rz2;

    //p=z+beta*p
    DoubleScalarxpy(pvec, beta, pvec, zvec, ndata);

  }

  et=gettimeofday_sec();
  t1=et-st;

  FileOutPutVec(p_x, xvec, ndata);
  t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);
  printf("Execution Time=%lf s\n", t1);

  Double1Free(Av);
  Double1Free(rvec);
  Double1Free(pvec);
  Double1Free(zvec);
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
