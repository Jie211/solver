#include "cr.h"

void CR_Init(double *v1, double *v2, double *v3, double *v4, double *x, double ndata){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,ndata);
  DoubleVecInit(v4,0.0,ndata);
  DoubleVecInit(x,0.0,ndata);
}

int CR_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, int nnz, double eps, int i_max){
  /* int i, j, k, n; */
  int loop;
  
  double *rvec, *pvec, *qvec, *svec, *x_0, error=0.0;
  double alpha, beta, bnorm, rnorm;
  double rs, rs2;

  bool flag=false;
  double t_error=0.0;
  FILE *p_x=NULL, *p_his=NULL;

  double st, et, t1=0.0;
  double st2, et2, t2=0.0;

  if(!INNER){
    p_x=FileInit("./output/CR_x.txt", "w");
    p_his=FileInit("./output/CR_his.txt", "w");
  }

  st=gettimeofday_sec();

  rvec=Double1Malloc(ndata);
  pvec=Double1Malloc(ndata);
  qvec=Double1Malloc(ndata);
  svec=Double1Malloc(ndata);
  x_0=Double1Malloc(ndata);

  double *d_qvec, *d_xvec, *d_pvec, *d_svec, *d_rvec;

  int ThreadPerBlock=128;
  int BlockPerGrid=(ndata-1)/(ThreadPerBlock/32)+1;
  /* BlockPerGrid=ceil((double)ndata/(double)ThreadPerBlock); */

  int DotThreadPerBlock=128;
  int DotBlockPerGrid=ceil((double)ndata/(double)ThreadPerBlock);

  if(cuda){
    checkCudaErrors( cudaMalloc((void **)&d_qvec, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_xvec, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_pvec, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_svec, sizeof(double)*ndata)  );
    checkCudaErrors( cudaMalloc((void **)&d_rvec, sizeof(double)*ndata)  );
  }

  //init vectory
  CR_Init(rvec, pvec, qvec, svec, xvec, ndata);

  DoubleVecCopy(x_0, xvec, ndata);

  // b 2norm
  bnorm = Double2Norm(bvec, ndata);

  //Ax
  if(cuda){
    /* checkCudaErrors( cudaDeviceSynchronize() ); */
    st2=gettimeofday_sec();

    checkCudaErrors( cudaMemcpy(d_xvec, xvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemset(d_qvec, 0, sizeof(double)*ndata)  );
    
    DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_xvec, d_qvec);

    checkCudaErrors( cudaPeekAtLastError()  );

    checkCudaErrors( cudaMemcpy(qvec, d_qvec, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
    /* checkCudaErrors(cudaDeviceSynchronize()); */
    et2=gettimeofday_sec();
    t2+=et2-st2;
  }else{
    DoubleMVMCSR(qvec, val, col, ptr, xvec, ndata);
  }

  //r=b-Ax
  DoubleVecSub(rvec, bvec, qvec, ndata);

  //p=r
  DoubleVecCopy(pvec, rvec, ndata);

  //q=Ap
  if(cuda){
    /* checkCudaErrors( cudaDeviceSynchronize() ); */
    st2=gettimeofday_sec();

    checkCudaErrors( cudaMemcpy(d_pvec, pvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemset(d_qvec, 0, sizeof(double)*ndata)  );

    DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_pvec, d_qvec);

    checkCudaErrors( cudaPeekAtLastError()  );

    checkCudaErrors( cudaMemcpy(qvec, d_qvec, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
    /* checkCudaErrors(cudaDeviceSynchronize()); */
    et2=gettimeofday_sec();
    t2+=et2-st2;
  }else{
    DoubleMVMCSR(qvec, val, col, ptr, pvec, ndata);
  }

  //s=q
  DoubleVecCopy(svec, qvec, ndata);

  // (r,s)
  if(cuda){
    /* checkCudaErrors( cudaDeviceSynchronize()  ); */
    st2=gettimeofday_sec();
    checkCudaErrors( cudaMemcpy(d_rvec, rvec, sizeof(double)*ndata, cudaMemcpyHostToDevice)  );
    checkCudaErrors( cudaMemcpy(d_svec, svec, sizeof(double)*ndata, cudaMemcpyHostToDevice)  );
    rs=DoubleCudaDot_Host(ndata, d_rvec, d_svec, DotBlockPerGrid, DotThreadPerBlock);
    /* checkCudaErrors( cudaDeviceSynchronize()  ); */
    et2=gettimeofday_sec();
    t2+=et2-st2;
  }else{
    rs=DoubleDot(rvec, svec, ndata); 
  }

  for(loop=0;loop<i_max;loop++){
    //rnorm
    rnorm = Double2Norm(rvec, ndata);
    error=rnorm/bnorm;
    if(!INNER){
      if(verbose){
        printf("%d %.12e\n",loop, error);
      }
      fprintf(p_his,"%d %.12e\n",loop, error);
    }
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
    if(cuda){
    /* checkCudaErrors( cudaDeviceSynchronize() ); */
    st2=gettimeofday_sec();

    checkCudaErrors( cudaMemcpy(d_rvec, rvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemset(d_svec, 0, sizeof(double)*ndata)  );

    DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_rvec, d_svec);

    checkCudaErrors( cudaPeekAtLastError()  );

    checkCudaErrors( cudaMemcpy(svec, d_svec, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
    /* checkCudaErrors(cudaDeviceSynchronize()); */
    et2=gettimeofday_sec();
    t2+=et2-st2;

    }else{
      DoubleMVMCSR(svec, val, col, ptr, rvec, ndata);
    }

    //(r,s)
    if(cuda){
      /* checkCudaErrors( cudaDeviceSynchronize()  ); */
      st2=gettimeofday_sec();
      rs2=DoubleCudaDot_Host(ndata, d_rvec, d_svec, DotBlockPerGrid, DotThreadPerBlock);
      /* checkCudaErrors( cudaDeviceSynchronize()  ); */
      et2=gettimeofday_sec();
      t2+=et2-st2;

    }else{
      rs2=DoubleDot(rvec, svec, ndata);
    }

    //beta=(r_new,s_new)/(r,s)
    beta = rs2/rs;

    rs=rs2;

    //p=r+beta*p
    DoubleScalarxpy(pvec, beta, pvec, rvec, ndata);

    //q=s+beta*q
    DoubleScalarxpy(qvec, beta, qvec, svec, ndata);

  }
  et=gettimeofday_sec();
  t1=et-st;

  if(!INNER){
    FileOutPutVec(p_x, xvec, ndata);
    t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    printf("Execution Time=%lf s\n", t1);
    if(cuda){
      printf("->Copy Time=%lf s\n", t2);
    }
  }
  if(INNER && verbose){
    printf("Inner %d %.12e\n", loop, error);
  }

  if(cuda){
    checkCudaErrors(cudaFree(d_xvec));
    checkCudaErrors(cudaFree(d_qvec));
    checkCudaErrors(cudaFree(d_pvec));
    checkCudaErrors(cudaFree(d_svec));
    checkCudaErrors(cudaFree(d_rvec));
  }

  Double1Free(rvec);
  Double1Free(pvec);
  Double1Free(qvec);
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
