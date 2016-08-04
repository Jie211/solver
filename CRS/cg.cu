#include "cg.h"

void CG_Init(double *v1, double *v2, double *v3, double *x, double ndata){
  DoubleVecInit(v1,0.0,ndata);
  DoubleVecInit(v2,0.0,ndata);
  DoubleVecInit(v3,0.0,ndata);
  DoubleVecInit(x,0.0,ndata);
}

int CG_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, int nnz, double eps, int i_max){
  /* int i, j, k, n; */
  int loop;

  double *rvec, *pvec, *Av, *x_0, dot, error=0.0;
  double alpha, beta, bnorm, rnorm;
  double rr, rr2;

  bool flag=false;

  double t_error=0.0;
  FILE *p_x=NULL, *p_his=NULL;

  double ct_start, ct_end, copytime=0.0;
  double dt_start, dt_end, dottime=0.0;
  double mvt_start, mvt_end, mvtime=0.0;

  if(!INNER){
    p_x=FileInit("./output/CG_x.txt", "w");
    p_his=FileInit("./output/CG_his.txt", "w");
  }



  /* Av=Double1Malloc(ndata); */
  /* rvec=Double1Malloc(ndata); */
  /* pvec=Double1Malloc(ndata); */
  /* x_0=Double1Malloc(ndata); */

  checkCudaErrors( cudaMallocHost((void **)&Av, sizeof(double)*ndata) );
  checkCudaErrors( cudaMallocHost((void **)&rvec, sizeof(double)*ndata) );
  checkCudaErrors( cudaMallocHost((void **)&pvec, sizeof(double)*ndata) );
  checkCudaErrors( cudaMallocHost((void **)&x_0, sizeof(double)*ndata) );


  double *d_Av=NULL, *d_xvec=NULL, *d_pvec=NULL;
  double *d_rvec=NULL;
  /* double *d_val=NULL; */
  /* int *d_col=NULL, *d_ptr=NULL; */

  int ThreadPerBlock=128;
  int BlockPerGrid=(ndata-1)/(ThreadPerBlock/32)+1;
  /* BlockPerGrid=ceil((double)ndata/(double)ThreadPerBlock); */

  int DotThreadPerBlock=128;
  int DotBlockPerGrid=ceil((double)ndata/(double)ThreadPerBlock);

  if(cuda){

  checkCudaErrors( cudaMalloc((void **)&d_Av, sizeof(double)*ndata) );
  checkCudaErrors( cudaMalloc((void **)&d_pvec, sizeof(double)*ndata) );
  checkCudaErrors( cudaMalloc((void **)&d_xvec, sizeof(double)*ndata) );
  checkCudaErrors( cudaMalloc((void **)&d_rvec, sizeof(double)*ndata) );

  }

  //init vectory
  CG_Init(rvec, pvec, Av, xvec, ndata);

  DoubleVecCopy(x_0, xvec, ndata);

  // b 2norm
  bnorm = Double2Norm(bvec, ndata);



  //Ax
  if(cuda){

    /* checkCudaErrors(cudaDeviceSynchronize()); */


    ct_start=gettimeofday_sec();
    checkCudaErrors( cudaMemcpy(d_xvec, xvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemset(d_Av, 0, sizeof(double)*ndata) );
    ct_end=gettimeofday_sec();
    copytime+=ct_end-ct_start;

    mvt_start=gettimeofday_sec();
    DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_xvec, d_Av);
    mvt_end=gettimeofday_sec();
    mvtime+=mvt_end-mvt_start;

    checkCudaErrors( cudaPeekAtLastError() );


    ct_start=gettimeofday_sec();
    checkCudaErrors( cudaMemcpy(Av, d_Av, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
    ct_end=gettimeofday_sec();
    copytime+=ct_end-ct_start;


    /* checkCudaErrors(cudaDeviceSynchronize()); */

  }else{
    DoubleMVMCSR(Av, val, col, ptr, xvec, ndata);
  }

  //r=b-Ax
  DoubleVecSub(rvec, bvec, Av, ndata);
  //p=r
  DoubleVecCopy(pvec, rvec, ndata);

  /* // r 2norm */
  /* rnorm = Double2Norm(rvec, ndata); */
  // (r,r)
  if(cuda){
    ct_start=gettimeofday_sec();
    checkCudaErrors( cudaMemcpy(d_rvec, rvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
    ct_end=gettimeofday_sec();
    copytime+=ct_end-ct_start;

    dt_start=gettimeofday_sec();
    dot=DoubleCudaDot_Host(ndata, d_rvec, d_rvec, DotBlockPerGrid, DotThreadPerBlock);
    dt_end=gettimeofday_sec();
    dottime+=dt_end-dt_start;

    rr=dot;

    /* rr=DoubleDot(rvec, rvec, ndata);  */
  }else{
    rr=DoubleDot(rvec, rvec, ndata); 
  }

  for(loop=0;loop<i_max;loop++){
    //rnorm
    rnorm = Double2Norm(rvec, ndata);
    error=rnorm/bnorm;
    if(!INNER){
      if(verbose){
        printf("%d %.12e\n",loop+1, error);
      }
      fprintf(p_his,"%d %.12e\n",loop+1, error);
    }
    if(error <= eps){
      flag=true;
      break;
    }

    //Ap
    if(cuda){
      /* checkCudaErrors(cudaDeviceSynchronize()); */


      ct_start=gettimeofday_sec();
      checkCudaErrors( cudaMemcpy(d_pvec, pvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
      checkCudaErrors( cudaMemset(d_Av, 0, sizeof(double)*ndata) );
      ct_end=gettimeofday_sec();
      copytime+=ct_end-ct_start;

      mvt_start=gettimeofday_sec();
      DoubleCudaMVMCSR<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*(ThreadPerBlock+16)>>>(ndata, d_val, d_col, d_ptr, d_pvec, d_Av);
      mvt_end=gettimeofday_sec();
      mvtime+=mvt_end-mvt_start;

      checkCudaErrors( cudaPeekAtLastError() );

      ct_start=gettimeofday_sec();
      checkCudaErrors( cudaMemcpy(Av, d_Av, sizeof(double)*ndata, cudaMemcpyDeviceToHost) );
      ct_end=gettimeofday_sec();
      copytime+=ct_end-ct_start;


      /* checkCudaErrors(cudaDeviceSynchronize()); */

    }else{
      DoubleMVMCSR(Av, val, col, ptr, pvec, ndata);
    }

    //alpha=(r,r)/(p,ap)
    if(cuda){
      
      /* checkCudaErrors(cudaDeviceSynchronize()); */

      dt_start=gettimeofday_sec();
      dot=DoubleCudaDot_Host(ndata, d_pvec, d_Av, DotBlockPerGrid, DotThreadPerBlock);
      dt_end=gettimeofday_sec();
      dottime+=dt_end-dt_start;

      /* checkCudaErrors(cudaDeviceSynchronize()); */

      /* dot=DoubleDot(pvec, Av, ndata); */
    }else{
      dot=DoubleDot(pvec, Av, ndata);
    }
    alpha = rr / dot;

    //x=x+alpha*pvec
    DoubleScalarxpy(xvec, alpha, pvec, xvec, ndata);

    //r=r-alpha*Ap
    DoubleScalarxpy(rvec, -alpha, Av, rvec, ndata);

    //(r,r)
    if(cuda){

      /* checkCudaErrors(cudaDeviceSynchronize()); */


      ct_start=gettimeofday_sec();
      checkCudaErrors( cudaMemcpy(d_rvec, rvec, sizeof(double)*ndata, cudaMemcpyHostToDevice) );
      ct_end=gettimeofday_sec();
      copytime+=ct_end-ct_start;

      dt_start=gettimeofday_sec();
      rr2=DoubleCudaDot_Host(ndata, d_rvec, d_rvec, DotBlockPerGrid, DotThreadPerBlock);
      dt_end=gettimeofday_sec();
      dottime+=dt_end-dt_start;


      /* checkCudaErrors(cudaDeviceSynchronize()); */

      /* rr2=DoubleDot(rvec, rvec, ndata); */
    }else{
      rr2=DoubleDot(rvec, rvec, ndata);
    }

    //beta=(r_new,r_new)/(r,r)
    beta = rr2/rr;

    rr=rr2;

    //p=r+beta*p
    DoubleScalarxpy(pvec, beta, pvec, rvec, ndata);

  }


  if(!INNER){
    FileOutPutVec(p_x, xvec, ndata);
    t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
    printf("|b-ax|2/|b|2=%.1f\n", t_error);
    /* printf("Execution Time=%lf s\n", t1); */
    if(cuda){
      /* printf("->Copy Time=%lf s\n", t2); */
      printf("Copy=%lf Dot=%lf MV=%lf\n", copytime, dottime, mvtime);
    }
  }
  if(INNER && verbose)
    printf("Inner %d %.12e\n", loop, error);


  if(cuda){

    /* checkCudaErrors( cudaFree(d_val) ); */
    /* checkCudaErrors( cudaFree(d_col) ); */
    /* checkCudaErrors( cudaFree(d_ptr) ); */
    checkCudaErrors( cudaFree(d_xvec) );
    checkCudaErrors( cudaFree(d_pvec) );
    checkCudaErrors( cudaFree(d_Av) );
    checkCudaErrors( cudaFree(d_rvec) );
  }

  /* cudaFree(Av); */
  /* cudaFree(rvec); */
  /* cudaFree(pvec); */
  /* cudaFree(x_0); */

  /* Double1Free(Av); */
  /* Double1Free(rvec); */
  /* Double1Free(pvec); */
  /* Double1Free(x_0); */

  checkCudaErrors(cudaFreeHost(Av));
  checkCudaErrors(cudaFreeHost(rvec));
  checkCudaErrors(cudaFreeHost(pvec));
  checkCudaErrors(cudaFreeHost(x_0));



  if(!INNER){
    FileClose(p_x);
    FileClose(p_his);
  }
  if(flag){
    return 1;
  }
  return 2;
}
