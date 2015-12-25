#include "kskipcr.h"
void KSKIPCR_Init(double *v1, double *v2, double *v3, double *v4, double *v5, double *v6, double *v7, double *v8, double *v9, int ndata, int kskip)
{
  DoubleVecInit(v1, 0.0, (2*kskip+1)*ndata);
  DoubleVecInit(v2, 0.0, (2*kskip+1)*ndata);
  DoubleVecInit(v3, 0.0, (2*kskip+1)*ndata);
  DoubleVecInit(v4, 0.0, (2*kskip+1)*ndata);
  DoubleVecInit(v5, 0.0, (2*kskip+1)*ndata);

  DoubleVecInit(v6, 0.0, ndata);
  DoubleVecInit(v7, 0.0, ndata);
  DoubleVecInit(v8, 0.0, ndata);
  DoubleVecInit(v9, 0.0, ndata);
}

int KSKIPCR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, int nnz, double eps, int i_max, int kskip, int fix)
{
  int nloop, iloop, jloop;

  double *Ar, *Ap, *delta, *eta, *zeta;
  double *rvec, *pvec, *Av, *x_0;
  double rnorm, bnorm, alpha, beta, error=0.0;
  bool flag=false;

  double t_error;
  FILE *p_x=NULL, *p_his=NULL;

  double st, et, t1;

  if(!INNER){
    p_x=FileInit("./output/KskipCR_x.txt", "w");
    p_his=FileInit("./output/KskipCR_his.txt", "w");
  }

  st=gettimeofday_sec();
  Ar=Double1Malloc((2*kskip+1)*ndata);
  Ap=Double1Malloc((2*kskip+1)*ndata);
  delta=Double1Malloc((2*kskip+1)*ndata);
  eta=Double1Malloc((2*kskip+1)*ndata);
  zeta=Double1Malloc((2*kskip+1)*ndata);

  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  Av=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  KSKIPCR_Init(Ar, Ap, delta, eta, zeta, rvec, pvec, Av, xvec, ndata, kskip);

  DoubleVecCopy(x_0, xvec, ndata);

  //Ax
  DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);

  //r=b-Ax
  DoubleVecSub(rvec, bvec, Av, ndata);

  //p=r
  DoubleVecCopy(pvec, rvec, ndata);

  bnorm=Double2Norm(bvec, ndata);

  for(nloop=0; nloop<i_max; nloop+=(kskip+1))
  {
    rnorm=Double2Norm(rvec, ndata);
    error=rnorm/bnorm;
    if(!INNER){
      if(verbose){
        printf("%d %.12e\n", nloop, error);
      }
      fprintf(p_his, "%d %.12e\n", nloop, error);
    }
    if(error<=eps){
      flag=true;
      break;
    }

    //Ar-> Ar^2k+1
    //Ap-> Ap^2k+1
    DoubleCalArApKCR(Ar, Ap, val, col, ptr, rvec, pvec, ndata, kskip);

    //delta=(r,Ar)
    //eta=(A1p,Ap)
    //zeta=(r,Ap)
    DoubleCalDeltaEtaZetaKCR(delta, eta, zeta, Ar, Ap, rvec, ndata, kskip);

    for(iloop=nloop; iloop<=nloop+kskip; iloop++)
    {
      //alpha=delta_1/eta_1
      alpha=delta[0]/eta[0];

      //beta (delta_1 - 2*alpha*zeta_2 + alpha^2*eta[2])/delta_1
      beta=(delta[0] - 2*alpha*zeta[1] + alpha*alpha*eta[1]) / delta[0];

      //update delta eta zeta
      for(jloop=0; jloop<2*kskip-2*(iloop-nloop); jloop++)
      {
        double delta2=0.0;
        delta[jloop] = delta[jloop] - 2*alpha*zeta[jloop+1] + alpha*alpha*eta[jloop+1];
        delta2=delta[jloop] - 2*alpha*zeta[jloop+1] + alpha*alpha*eta[jloop+1];
        /* eta[jloop] = delta[jloop+1] + 2*beta*(zeta[jloop+1]-alpha*eta[jloop+1]) + beta*beta*eta[jloop]; */
        eta[jloop] = delta2 + 2*beta*(zeta[jloop+1]-alpha*eta[jloop+1]) + beta*beta*eta[jloop];
        zeta[jloop] = delta[jloop] - alpha*zeta[jloop+1] - alpha*(zeta[jloop+1]-alpha*eta[jloop+1]) + beta*zeta[jloop] - alpha*beta*eta[jloop];
      }

      //new Ap
      DoubleMVMCSR(Av, val, col, ptr, pvec, ndata);

      //x=x+alpha*p
      DoubleScalarxpy(xvec, alpha, pvec, xvec, ndata);

      //r=r-alpha*Ap
      DoubleScalarxpy(rvec, -alpha, Av, rvec, ndata);

      //p=r+beta*p
      DoubleScalarxpy(pvec, beta, pvec, rvec, ndata);
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
  if(INNER){
    printf("Inner %d %.12e\n",nloop,error);
  }
  Double1Free(Ar);
  Double1Free(Ap);
  Double1Free(delta);
  Double1Free(eta);
  Double1Free(zeta);
  Double1Free(rvec);
  Double1Free(pvec);
  Double1Free(Av);
  Double1Free(x_0);
  if(!INNER){
    FileClose(p_x);
    FileClose(p_his);
  }
  if(flag){
    return 1;
  }
  return 2;
}
