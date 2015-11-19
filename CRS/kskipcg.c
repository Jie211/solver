#include "kskipcg.h"

void KSKIPCG_Init(double *v1, double *v2, double *v3, double *v4, double *v5, double *v6, double *v7, double *v8, double *v9, int ndata, int kskip)
{
  DoubleVecInit(v1, 0.0, 2*kskip*ndata);
  DoubleVecInit(v2, 0.0, (2*kskip+2)*ndata);
  DoubleVecInit(v3, 0.0, 2*kskip);
  DoubleVecInit(v4, 0.0, 2*kskip+1);
  DoubleVecInit(v5, 0.0, 2*kskip+2);

  DoubleVecInit(v6, 0.0, ndata);
  DoubleVecInit(v7, 0.0, ndata);
  DoubleVecInit(v8, 0.0, ndata);
  DoubleVecInit(v9, 0.0, ndata);
}

int KSKIPCG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, double eps, int i_max, int kskip, int fix)
{
  int nloop, iloop, jloop;

  double *Ar, *Ap, *delta, *eta, *zeta;
  double *rvec, *pvec, *Av, *x_0;
  double rnorm, bnorm, alpha, beta, gamma;
  bool flag=false;
  double error=0.0;

#ifndef INNER
  double t_error;
  FILE *p_x, *p_his;
  p_x=FileInit("./output/KskipCG_x.txt", "w");
  p_his=FileInit("./output/KskipCG_his.txt", "w");
#endif

  Ar=Double1Malloc(2*kskip*ndata);
  Ap=Double1Malloc((2*kskip+2)*ndata);
  delta=Double1Malloc(2*kskip);
  eta=Double1Malloc(2*kskip+1);
  zeta=Double1Malloc(2*kskip+2);

  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  Av=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  KSKIPCG_Init(Ar, Ap, delta, eta, zeta, rvec, pvec, Av, xvec, ndata, kskip);

  DoubleVecCopy(x_0, xvec, ndata);

  //Ax
  DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);
  
  //r=b-Ax
  DoubleVecSub(rvec, bvec, Av, ndata);

  //p=r
  DoubleVecCopy(pvec, rvec, ndata);
  
  //b 2norm
  bnorm=Double2Norm(bvec, ndata);

  for(nloop=0; nloop<i_max; nloop+=(kskip+1))
  {
    rnorm=Double2Norm(rvec, ndata);
    error=rnorm/bnorm;
#ifndef INNER
    printf("%d %.12e\n", nloop, error);
    fprintf(p_his,"%d %.12e\n", nloop, error);
#endif
    if(error<=eps){
      flag=true;
      break;
    }

    //Ar-> Ar^2k
    //Ap-> Ap^2k+2
    DoubleCalArApKCG(Ar, Ap, val, col, ptr, rvec, pvec, ndata, kskip);

    //gamma=(r,r)
    gamma=DoubleDot(rvec, rvec, ndata);

    //delta=(r,Ar)
    //eta=(r,Ap)
    //zeta=(p,Ap)
    DoubleCalDeltaEtaZetaKCG(delta, eta, zeta, Ar, Ap, rvec, pvec, ndata, kskip);

    for(iloop=nloop; iloop<=nloop+kskip; iloop++)
    {
      //alpha = gamma/zeta_1
      alpha=gamma/zeta[0];
      //beta = (alpha*zeta_2/zeta_1) -1
      beta=alpha*zeta[1]/zeta[0] - 1.0;
      if(fix==1){
        gamma=beta*gamma;
      }else if(fix==2){
        /* double tmp1 = alpha*alpha*zeta[1]; */
        /* double tmp2 = gamma; */
        /* gamma=tmp1 - tmp2; */
        double tmp0 = gamma - alpha * eta[0];
        double tmp1 = eta[0] - alpha*zeta[1];
        gamma = tmp0 - alpha*tmp1;
      }

      //update delta eta zeta
      for(jloop=0; jloop<2*kskip-2*(iloop-nloop);jloop++){
        delta[jloop] = delta[jloop] - 2*alpha*eta[jloop+1] + alpha*alpha*eta[jloop+2];
        double eta_old=eta[jloop];
        eta[jloop] = delta[jloop] + beta*zeta[jloop+1] - alpha*beta*zeta[jloop+1];
        zeta[jloop] = eta[jloop+1] + beta*eta_old + beta*beta*zeta[jloop] - alpha*beta*zeta[jloop+1];
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
#ifndef INNER
  FileOutPutVec(p_x, xvec, ndata);
  t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);
#endif
#ifdef INNER
  printf("Inner %d %.12e\n",nloop,error);
#endif
  Double1Free(Ar);
  Double1Free(Ap);
  Double1Free(delta);
  Double1Free(eta);
  Double1Free(zeta);
  Double1Free(rvec);
  Double1Free(pvec);
  Double1Free(Av);
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
