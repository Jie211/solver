#include "kskipcg.h"

/* void KSKIPCG_Init(double *v1, double *v2, double *v3, double *v4, double *v5, double *v6, double *v7, double *v8, double *v9, int ndata, int kskip) */
/* void KSKIPCG_Init(double *v1, double **v2, double *v3, double *v4, double *v5, double *v6, double *v7, double *v8, double *v9, int ndata, int kskip) */
void KSKIPCG_Init(double **v1, double **v2, double *v3, double *v4, double *v5, double *v6, double *v7, double *v8, double *v9, int ndata, int kskip)
{
  /* DoubleVecInit(v1, 0.0, 2*kskip*ndata); */
  Double2VecInit(v1, 0.0, ndata, 2*kskip+1);
  /* DoubleVecInit(v2, 0.0, (2*kskip+2)*ndata); */
  Double2VecInit(v2, 0.0, ndata, (2*kskip+2));
  DoubleVecInit(v3, 0.0, 2*kskip);
  DoubleVecInit(v4, 0.0, 2*kskip+1);
  DoubleVecInit(v5, 0.0, 2*kskip+2);

  DoubleVecInit(v6, 0.0, ndata);
  DoubleVecInit(v7, 0.0, ndata);
  DoubleVecInit(v8, 0.0, ndata);
  DoubleVecInit(v9, 0.0, ndata);
}

int KSKIPCG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, int nnz, double eps, int i_max, int kskip, int fix)
{

  if(verbose)
    printf("get %d loop!\n", i_max);
  int nloop, iloop, jloop;
  int i,  ii;

  double *delta, *eta, *zeta;
  double **Ap;
  /* double *Ap; */
  double **Ar;
  /* double *Ar; */
  double *rvec, *pvec, *Av, *x_0;
  double rnorm, bnorm, alpha, beta, gamma;
  bool flag=false;
  double error=0.0;

  double st, et, t1=0.0;

  double cp_s, cp_e, cp_t=0;
  double mv_s, mv_e, mv_t=0;
  double dot_s, dot_e, dot_t=0;

  double t_error;
  FILE *p_x=NULL, *p_his=NULL;

  double dot=0.0, tmp1=0.0, tmp2=0.0, tmp3=0.0;
  if(!INNER){
    p_x=FileInit("./output/KskipCG_x.txt", "w");
    p_his=FileInit("./output/KskipCG_his.txt", "w");
  }


  /* Ar=Double1Malloc(2*kskip*ndata); */
  Ar=Double2Malloc(ndata, 2*kskip+1);
  /* Ap=Double1Malloc((2*kskip+2)*ndata); */
  Ap=Double2Malloc(ndata, (2*kskip+2));

  delta=Double1Malloc(2*kskip);
  eta=Double1Malloc(2*kskip+1);
  zeta=Double1Malloc(2*kskip+2);

  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  Av=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);


  double *d_Ar, *d_Ap, *d_pvec, *d_rvec, *d_tmp, *d_xvec, *d_Av;
  int ThreadPerBlock=128;
  int BlockPerGrid=(ndata-1)/(ThreadPerBlock/32)+1;
  /* BlockPerGrid=ceil((double)ndata/(double)ThreadPerBlock); */

  int DotThreadPerBlock=128;
  int DotBlockPerGrid=(ndata-1)/(DotThreadPerBlock)+1;

  if(cuda){
    checkCudaErrors( cudaMalloc((void **)&d_Ar, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_Ap, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_pvec, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_rvec, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_tmp, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_xvec, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_Av, sizeof(double)*ndata)  );

  }

  st=gettimeofday_sec();

  KSKIPCG_Init(Ar, Ap, delta, eta, zeta, rvec, pvec, Av, xvec, ndata, kskip);

  DoubleVecCopy(x_0, xvec, ndata);

  //Ax
  if(cuda){
    cp_s=gettimeofday_sec();

    checkCudaErrors( cudaMemcpy(d_xvec, xvec, sizeof(double)*ndata, cudaMemcpyHostToDevice)  );
    checkCudaErrors( cudaMemset(d_Av, 0, sizeof(double)*ndata)  );

    cp_e=gettimeofday_sec();
    cp_t+=cp_e-cp_s;

    mv_s=gettimeofday_sec();
    DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_xvec, d_Av);
    mv_e=gettimeofday_sec();
    mv_t+=mv_e-mv_s;

    checkCudaErrors( cudaPeekAtLastError() );

    cp_s=gettimeofday_sec();
    checkCudaErrors( cudaMemcpy(Av, d_Av, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
    cp_e=gettimeofday_sec();
    cp_t+=cp_e-cp_s;


  }else{
    DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);
  }

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
    if(!INNER){
      if(verbose){
        printf("%d %.12e\n", nloop, error);
      }
      fprintf(p_his,"%d %.12e\n", nloop, error);
    }
    if(error<=eps){
      flag=true;
      break;
    }

    //Ar-> Ar^2k
    //Ap-> Ap^2k+2
    /* DoubleCalArApKCG(Ar, Ap, val, col, ptr, rvec, pvec, ndata, kskip); */

    if(cuda){
      cp_s=gettimeofday_sec();

      checkCudaErrors(cudaMemcpy(d_pvec, pvec, sizeof(double)*ndata, cudaMemcpyHostToDevice));
      checkCudaErrors(cudaMemcpy(d_rvec, rvec, sizeof(double)*ndata, cudaMemcpyHostToDevice));
      checkCudaErrors(cudaMemset(d_Ar, 0, sizeof(double)*ndata));
      checkCudaErrors(cudaMemset(d_Ap, 0, sizeof(double)*ndata));

      cp_e=gettimeofday_sec();
      cp_t+=cp_e-cp_s;

      mv_s=gettimeofday_sec();
      DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_rvec, d_Ar);
      mv_e=gettimeofday_sec();
      mv_t+=mv_e-mv_s;

      checkCudaErrors(cudaPeekAtLastError());

      cp_s=gettimeofday_sec();
      checkCudaErrors(cudaMemcpy(Ar[0], d_Ar, sizeof(double)*ndata, cudaMemcpyDeviceToHost));
      cp_e=gettimeofday_sec();
      cp_t+=cp_e-cp_s;


      mv_s=gettimeofday_sec();
      DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_pvec, d_Ap);
      mv_e=gettimeofday_sec();
      mv_t+=mv_e-mv_s;

      checkCudaErrors(cudaPeekAtLastError());

      cp_s=gettimeofday_sec();
      checkCudaErrors(cudaMemcpy(Ap[0], d_Ap, sizeof(double)*ndata, cudaMemcpyDeviceToHost));
      cp_e=gettimeofday_sec();
      cp_t+=cp_e-cp_s;

      for(ii=1;ii<2*kskip;ii++){
        cp_s=gettimeofday_sec();
        checkCudaErrors( cudaMemcpy(d_Ar, Ar[ii-1], sizeof(double)*ndata, cudaMemcpyHostToDevice)  );
        checkCudaErrors( cudaMemset(d_tmp, 0, sizeof(double)*ndata)  );
        cp_e=gettimeofday_sec();
        cp_t+=cp_e-cp_s;



        mv_s=gettimeofday_sec();
        DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_Ar, d_tmp);
        mv_e=gettimeofday_sec();
        mv_t+=mv_e-mv_s;

        checkCudaErrors( cudaPeekAtLastError() );

        cp_s=gettimeofday_sec();
        checkCudaErrors( cudaMemcpy(Ar[ii], d_tmp, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
        cp_e=gettimeofday_sec();
        cp_t+=cp_e-cp_s;

      }

      for(ii=1;ii<2*kskip+2;ii++){

        cp_s=gettimeofday_sec();
        checkCudaErrors( cudaMemcpy(d_Ap, Ap[ii-1], sizeof(double)*ndata, cudaMemcpyHostToDevice)  );
        checkCudaErrors( cudaMemset(d_tmp, 0, sizeof(double)*ndata)  );
        cp_e=gettimeofday_sec();
        cp_t+=cp_e-cp_s;

        mv_s=gettimeofday_sec();
        DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_Ap, d_tmp);
        mv_e=gettimeofday_sec();
        mv_t+=mv_e-mv_s;

        checkCudaErrors( cudaPeekAtLastError() );

        cp_s=gettimeofday_sec();
        checkCudaErrors( cudaMemcpy(Ap[ii], d_tmp, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
        cp_e=gettimeofday_sec();
        cp_t+=cp_e-cp_s;

      }

    }else{
      DoubleCalArApKCG(Ar, Ap, val, col, ptr, rvec, pvec, ndata, kskip);
    }


    //gamma=(r,r)
    if(cuda){

      dot_s=gettimeofday_sec();
      dot=DoubleCudaDot_Host(ndata, d_rvec, d_rvec, DotBlockPerGrid, DotThreadPerBlock);
      dot_e=gettimeofday_sec();
      dot_t+=dot_e-dot_s;

      gamma=dot;
    }else{
      gamma=DoubleDot(rvec, rvec, ndata);
    }

    //delta=(r,Ar)
    //eta=(r,Ap)
    //zeta=(p,Ap)
    if(cuda){


      for(i=0;i<2*kskip+2;i++){
        cp_s=gettimeofday_sec();
        checkCudaErrors( cudaMemcpy(d_Ap, Ap[i], sizeof(double)*ndata, cudaMemcpyHostToDevice)  );
        cp_e=gettimeofday_sec();
        cp_t+=cp_e-cp_s;
        if(i<2*kskip){

          cp_s=gettimeofday_sec();
          checkCudaErrors( cudaMemcpy(d_Ar, Ar[i], sizeof(double)*ndata, cudaMemcpyHostToDevice)  );
          cp_e=gettimeofday_sec();
          cp_t+=cp_e-cp_s;

          dot_s=gettimeofday_sec();
          tmp1=DoubleCudaDot_Host(ndata, d_rvec, d_Ar, DotBlockPerGrid, DotThreadPerBlock);
          dot_e=gettimeofday_sec();
          dot_t+=dot_e-dot_s;
          delta[i]=tmp1;
        }
        if(i<2*kskip+1){
          dot_s=gettimeofday_sec();
          tmp2=DoubleCudaDot_Host(ndata, d_rvec, d_Ap, DotBlockPerGrid, DotThreadPerBlock);
          dot_e=gettimeofday_sec();
          dot_t+=dot_e-dot_s;
          eta[i]=tmp2;
        }
        dot_s=gettimeofday_sec();
        tmp3=DoubleCudaDot_Host(ndata, d_pvec, d_Ap, DotBlockPerGrid, DotThreadPerBlock);
        dot_e=gettimeofday_sec();
        dot_t+=dot_e-dot_s;
        zeta[i]=tmp3;
      }

    }else{
      DoubleCalDeltaEtaZetaKCG(delta, eta, zeta, Ar, Ap, rvec, pvec, ndata, kskip);
    }


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
      if(cuda){
        cp_s=gettimeofday_sec();
        checkCudaErrors( cudaMemcpy(d_pvec, pvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMemset(d_Av, 0, sizeof(double)*ndata)  );
        cp_e=gettimeofday_sec();
        cp_t+=cp_e-cp_s;

        mv_s=gettimeofday_sec();
        DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_pvec, d_Av);
        mv_e=gettimeofday_sec();
        mv_t+=mv_e-mv_s;

        checkCudaErrors( cudaPeekAtLastError() );

        cp_s=gettimeofday_sec();
        checkCudaErrors( cudaMemcpy(Av, d_Av, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
        cp_e=gettimeofday_sec();
        cp_t+=cp_e-cp_s;


      }else{
        DoubleMVMCSR(Av, val, col, ptr, pvec, ndata);
      }

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
    if(cuda){
      printf("copy = %lf s\n", cp_t);
      printf("mv = %lf s\n", mv_t);
      printf("dot = %lf s\n", dot_t);
    }
  }
  if(INNER && verbose){
    printf("Inner %d %.12e\n",nloop,error);
  }

  if(cuda){

    checkCudaErrors( cudaFree(d_Ar)  );
    checkCudaErrors( cudaFree(d_Ap)  );
    checkCudaErrors( cudaFree(d_pvec)  );
    checkCudaErrors( cudaFree(d_rvec)  );
    checkCudaErrors( cudaFree(d_tmp)  );
    checkCudaErrors( cudaFree(d_xvec)  );
    checkCudaErrors( cudaFree(d_Av)  );
  }

  /* Double1Free(Ar); */
  Double2Free(Ar,2*kskip+1);
  /* Double1Free(Ap); */
  Double2Free(Ap,2*kskip+2);
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
