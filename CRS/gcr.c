#include "gcr.h"

void GCR_Init(double *v1, double *v2, double *v3, double **v4, double **v5, double *x, int ndata, int restart){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,restart);
  Double2VecInit(v4,0.0,ndata, restart);
  Double2VecInit(v5,0.0,ndata, restart);
  DoubleVecInit(x,0.0,ndata);
}

int GCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max, int restart){

  int loop=0, kloop, iloop;

  double *rvec, *Av;
  double **qvec, **pvec, *qq , *x_0;
  double alpha, beta, rnorm, bnorm;
  bool flag=false;
  bool out=false;
  double error=0.0;

  double t_error=0.0;
  FILE *p_x=NULL, *p_his=NULL;

  double st, et, t1;

  if(!INNER){
    p_x=FileInit("./output/GCR_x.txt", "w");
    p_his=FileInit("./output/GCR_his.txt", "w");
  }

  st=gettimeofday_sec();

  rvec=Double1Malloc(ndata);
  Av=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  qq=Double1Malloc(restart);
  if(!INNER && verbose){
    printf("---- restart set to %d ----\n", restart);
  }
  qvec=Double2Malloc(ndata, restart);
  pvec=Double2Malloc(ndata, restart);

  GCR_Init(rvec, Av, qq, qvec, pvec, xvec, ndata, restart);

  DoubleVecCopy(x_0, xvec, ndata);

  //b 2norm
  bnorm=Double2Norm(bvec, ndata);

  while(loop<i_max){
    //Ax
    DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);

    //r=b-Ax
    DoubleVecSub(rvec, bvec, Av, ndata);

    //p=r
    DoubleVecCopy(pvec[0], rvec, ndata);

    //Ap
    DoubleMVMCSR(qvec[0], val, col, ptr, pvec[0], ndata);

    for(kloop=0;kloop<restart;kloop++){
      rnorm=Double2Norm(rvec, ndata);
      error=rnorm/bnorm;
      if(!INNER ){
        if(verbose){
          printf("Outer %d %.12e\n",loop, error);
        }
        fprintf(p_his,"%d %.12e\n",loop, error);
      }
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

      //Ar
      DoubleMVMCSR(Av, val, col, ptr, rvec, ndata);

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
      DoubleVecAdd(pvec[kloop+1], rvec, pvec[kloop+1], ndata);
      //q[k]=Az+q[k]
      DoubleVecAdd(qvec[kloop+1], Av, qvec[kloop+1], ndata);
    }
    if(out){
      break;
    }
  }
  et=gettimeofday_sec();
  t1=et-st;
  if(!INNER){
    FileOutPutVec(p_x, xvec, ndata);
    t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("Execution Time=%lf s\n", t1);
  }
  if(INNER && verbose){
    printf("Inner %d %.12e\n", loop, error);
  }
  Double1Free(rvec);
  Double1Free(Av);
  Double1Free(x_0);
  Double1Free(qq);
  Double2Free(qvec, restart);
  Double2Free(pvec, restart);
  if(!INNER){
    FileClose(p_x);
    FileClose(p_his);
  }

  if(flag){
    return 1;
  }else{
    return 2;
  }
  return 0;

}
