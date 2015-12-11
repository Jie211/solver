#include "vpgcr.h"

void VPGCR_Init(double *v1, double *v2, double *v3, double *v4, double **v5, double **v6, double *x, int ndata, int restart){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,ndata);
  DoubleVecInit(v4,0.0,restart);
  Double2VecInit(v5,0.0,ndata, restart);
  Double2VecInit(v6,0.0,ndata, restart);
  DoubleVecInit(x,0.0,ndata);
}

int VPGCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max, int restart){

  int loop=0, kloop, iloop;

  double *rvec, *zvec, *Av;
  double **qvec, **pvec, *qq, *x_0;
  double alpha, beta, rnorm, bnorm;
  bool flag=false;
  bool out=false;
  double t_error=0.0;
  int error_message;
  double error=0.0;

  double st, et, t1;

  FILE *p_x, *p_his;
  p_x=FileInit("./output/VPGCR_x.txt", "w");
  p_his=FileInit("./output/VPGCR_his.txt", "w");

  st=gettimeofday_sec();

  rvec=Double1Malloc(ndata);
  zvec=Double1Malloc(ndata);
  Av=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  qq=Double1Malloc(restart);
  if(verbose){
    printf("---- restart set to %d ----\n", restart);
  }
  qvec=Double2Malloc(ndata, restart);
  pvec=Double2Malloc(ndata, restart);

  VPGCR_Init(rvec, zvec, Av, qq, qvec, pvec, xvec, ndata, restart);

  DoubleVecCopy(x_0, xvec, ndata);

  //b 2norm
  bnorm=Double2Norm(bvec, ndata);

  while(loop<i_max){
    //Ax
    DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);

    //r=b-Ax
    DoubleVecSub(rvec, bvec, Av, ndata);

    //init pvec[0]
    DoubleVecInit(pvec[0], 0.0, ndata);
    //solve p by Ap[0]=r
    error_message=InnerSolverSelecter(val, col, ptr, rvec, pvec[0], ndata, eps_inner, loop_inner, kskip_inner, fix_inner);
    if(error_message!=0){
      printf("** error in vpgcr.c **\n");
      return -1;
    }

    //q[0*ndata+x]=A*p[0*ndata+x]
    DoubleMVMCSR(qvec[0], val, col, ptr, pvec[0], ndata);

    for(kloop=0;kloop<restart;kloop++){
      rnorm=Double2Norm(rvec, ndata);
      error=rnorm/bnorm;
      if(verbose){
        printf("Outer %d %.12e\n",loop, error);
      }
      fprintf(p_his,"%d %.12e\n",loop, error);
      if(error <= eps){
        flag=true;
        out=true;
        break;
      }else if(loop >= i_max){
        out=true;
        break;
      }
      loop++;

      //(q,q)
      qq[kloop]=DoubleDot(qvec[kloop], qvec[kloop], ndata);

      //alpha=(r, q)/(q, q)
      alpha=DoubleDot(rvec, qvec[kloop], ndata) / qq[kloop];

      //x=alpha*pvec[k]+xvec
      DoubleScalarxpy(xvec, alpha, pvec[kloop], xvec, ndata);
      if(kloop==restart-1){
        break;
      }

      //r=-alpha*qvec+rvec
      DoubleScalarxpy(rvec, -alpha, qvec[kloop], rvec, ndata);

      //init z
      DoubleVecInit(zvec, 0.0, ndata);
      //Az=r
      error_message=InnerSolverSelecter(val, col, ptr, rvec, zvec, ndata, eps_inner, loop_inner, kskip_inner, fix_inner);
      if(error_message!=0){
        printf("** error in vpgcr **\n");
        return -1;
      }

      //Az
      DoubleMVMCSR(Av, val, col, ptr, zvec, ndata);

      //init p[k+1] q[k+1]
      DoubleVecInit(pvec[kloop+1], 0.0, ndata);
      DoubleVecInit(qvec[kloop+1], 0.0, ndata);

      for(iloop=0; iloop<=kloop; iloop++){
        //beta=-(Az,qvec)/(q,q)
        beta= -(DoubleDot(Av, qvec[iloop], ndata)/qq[iloop]);
        //pvec[k+1]=beta*pvec[i]+pvec[k+1]
        DoubleScalarxpy(pvec[kloop+1], beta, pvec[iloop], pvec[kloop+1], ndata);
        //qvec[k+1]=beta*qvec[i]+qvec[k+1]
        DoubleScalarxpy(qvec[kloop+1], beta, qvec[iloop], qvec[kloop+1], ndata);
      }

      //p[k]=z+p[k]
      DoubleVecAdd(pvec[kloop+1], zvec, pvec[kloop+1], ndata);
      //q[k]=Az+q[k]
      DoubleVecAdd(qvec[kloop+1], Av, qvec[kloop+1], ndata);
    }
    if(out){
      break;
    }
  }

  et=gettimeofday_sec();
  t1=et-st;

  FileOutPutVec(p_x, xvec, ndata);
  t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);
  printf("Execution Time=%lf", t1);

  Double1Free(rvec);
  Double1Free(zvec);
  Double1Free(Av);
  Double1Free(x_0);
  Double1Free(qq);
  Double2Free(qvec, restart);
  Double2Free(pvec, restart);
  
  if(flag){
    return 1;
  }else{
    return 2;
  }
  return 0;

}
