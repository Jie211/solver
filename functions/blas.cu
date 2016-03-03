#include "blas.h"
void Display_Mes(const char *mes){
  printf("---- %s ----\n", mes);
}
void Display_Err(const char *err){
  printf("** %s **\n", err);
}
double gettimeofday_sec(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}
double *Double1Malloc(int ndata)
{
  double *tmp;
  tmp=(double *)malloc(sizeof(double)*ndata);
  if(!tmp){
    printf("double 1dim vec malloc error\n");
    exit(-1);
  }
  return (tmp);
}
double **Double2Malloc(int ndatax, int ndatay)
{
  int i;
  double **tmp;
  tmp=(double **)malloc(sizeof(double*)*ndatay);
  if(!tmp){
    printf("double 2dim vec malloc error\n");
    exit(-1);
  }
  for(i=0;i<ndatay;i++){
    tmp[i]=(double *)malloc(sizeof(double)*ndatax);
    if(!tmp[i]){
      printf("double 2dim vec malloc error\n");
      exit(-1);
    }
  }
  return (tmp);
}
void Double1Free(double *ptr)
{
  free(ptr);
}
void Double2Free(double **ptr, int ndatay)
{
  int i;
  for(i=0;i<ndatay;i++){
    free(ptr[i]);
  }
  free(ptr);
}
int *Intger1Malloc(int ndata)
{
  int *tmp;
  tmp=(int *)malloc(sizeof(int)*ndata);
  if(!tmp){
    printf("int 1dim vec malloc error\n");
    exit(-1);
  }
  return (tmp);
}
int **Intger2Malloc(int ndatax, int ndatay)
{
  int i;
  int **tmp;
  tmp=(int **)malloc(sizeof(int*)*ndatay);
  if(!tmp){
    printf("int 2dim vec malloc error\n");
    exit(-1);
  }
  for(i=0;i<ndatay;i++){
    tmp[i]=(int *)malloc(sizeof(int)*ndatax);
    if(!tmp[i]){
      printf("int 2dim vec malloc error\n");
      exit(-1);
    }
  }
  return (tmp);
}
void Intger1Free(int *ptr)
{
  free(ptr);
}
void Intger2Free(int **ptr, int ndata)
{
  int i;
  for(i=0;i<ndata;i++){
    free(ptr[i]);
  }
  free(ptr);
}
void DoubleVecAdd(double *out, double *x, double *y, int ndata)
{
  int i;
  for(i=0;i<ndata;i++){
    out[i]=x[i]+y[i];
  }
}
void DoubleVecSub(double *out, double *x, double *y, int ndata)
{
  int i;
  for(i=0;i<ndata;i++){
    out[i]=x[i]-y[i];
  }
}
void DoubleVecMul(double *out, double *x, double *y, int ndata)
{
  int i;
  for(i=0;i<ndata;i++){
    out[i]=x[i]*y[i];
  }
}
void DoubleScalar(double *out, double a, double *x, int ndata)
{
  int i;
  for(i=0;i<ndata;i++){
    out[i]=a*x[i];
  }
}
void DoubleScalarxpy(double *out, double a, double *x, double *y, int ndata)
{
  int i;
  double tmp;
  for(i=0;i<ndata;i++){
    tmp=y[i];
    out[i]=(a*x[i])+tmp;
  }
}
double DoubleDot(double *x, double *y, int ndata)
{
  double tmp=0.0;
  int i;
  for(i=0;i<ndata;i++){
    tmp += x[i]*y[i];
  }
  return tmp; 
}
double Double1Norm(double *x, int ndata)
{
  int i;
  double tmp=0.0;
  for(i=0;i<ndata;i++){
    tmp+=fabs(x[i]);
  }
  return tmp;
}
double Double2Norm(double *x, int ndata)
{
  int i;
  double tmp=0.0;
  for(i=0;i<ndata;i++){
    tmp+=x[i]*x[i];
  }
  return sqrt(tmp);
}
void DoubleMVMCSR(double *out, double *val, int *col, int *ptr, double *vec, int ndata)
{
  int i,j;
  double tmp=0.0;
#pragma omp parallel for private(j) reduction(+:tmp) schedule(static) firstprivate(out, val, vec) lastprivate(out)
  for(i=0;i<ndata;i++){
    tmp=0.0;
    for(j=ptr[i];j<ptr[i+1];j++){
      tmp+=val[j] * vec[col[j]];
    }
    out[i]=tmp;
  }
}
double DoubleMaxNorm(double *x, int ndata)
{
  int i;
  double max=fabs(x[0]);
  double m=0.0;
  for(i=1;i<ndata;i++){
    m=fabs(x[i]);
    if(m>max){
      max=m;
    }
  }
  return max;
}
void DoubleVecCopy(double *a, double *b, int ndata)
{
  int i;
  for(i=0;i<ndata;i++){
    a[i]=b[i];
  }
}
double error_check_CRS(double *val, const int *col, const int *ptr, double *b, double *x_new ,double *x_0, int N)
{
  double tmp1=0.0;
  double tmp2=0.0;
  double tmp=0.0;
  double *Ax, *Ax_0;
  Ax=(double *)malloc(sizeof(double)*N);
  Ax_0=(double *)malloc(sizeof(double)*N);

  int i, j;

  for(i=0;i<N;i++){
    tmp=0.0;
    for(j=ptr[i];j<ptr[i+1];j++){
      tmp+=val[j]*x_new[col[j]];
    }
    Ax[i]=b[i]-tmp;
  }
  for(i=0;i<N;i++){
    tmp=0.0;
    for(j=ptr[i];j<ptr[i+1];j++){
      tmp+=val[j]*x_0[col[j]];
    }
    Ax_0[i]=b[i]-tmp;
  }
  tmp1=Double2Norm(Ax, N);
  tmp2=Double2Norm(Ax_0, N);
  tmp=log10(tmp1/tmp2);
  free(Ax);
  free(Ax_0);
  return tmp;
}
void DoubleVecInit(double *vec, double val, int ndata)
{
  int i;
  for(i=0;i<ndata;i++){
    vec[i]=val;
  }
}
void Double2VecInit(double **vec, double val, int ndatax, int ndatay)
{
  int i, j;
  for(i=0;i<ndatay;i++){
    for(j=0;j<ndatax;j++){
      vec[i][j]=val;
    }
  }
}
/* void DoubleCalArApKCG(double *Ar, double *Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, int ndata, int kskip) */
/* void DoubleCalArApKCG(double *Ar, double **Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, int ndata, int kskip) */
void DoubleCalArApKCG(double **Ar, double **Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, int ndata, int kskip)
{
  int i, j, ii;
  double tmp1=0.0;
  double tmp2=0.0;
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val, pvec, rvec) lastprivate(Ar, Ap)
for(i=0;i<ndata;i++){
    tmp1=0.0;
    tmp2=0.0;
    for(j=ptr[i];j<ptr[i+1];j++){
      tmp1 += val[j]*rvec[col[j]];
      tmp2 += val[j]*pvec[col[j]];
    }
    /* Ar[0*ndata+i]=tmp1; */
    Ar[0][i]=tmp1;
    /* Ap[0*ndata+i]=tmp2; */
    Ap[0][i]=tmp2;
  }
  for(ii=1;ii<2*kskip+2;ii++){
#pragma omp parallel for private(i, j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val) lastprivate(Ar, Ap)
    for(i=0;i<ndata;i++){
      tmp1=0.0;
      tmp2=0.0;
      for(j=ptr[i];j<ptr[i+1];j++){
        if(ii<2*kskip){
          /* tmp1 += val[j]*Ar[(ii-1)*ndata+col[j]]; */
          tmp1 += val[j]*Ar[(ii-1)][col[j]];
        }
        /* tmp2 += val[j]*Ap[(ii-1)*ndata+col[j]]; */
        tmp2 += val[j]*Ap[(ii-1)][col[j]];
      }
      if(ii<2*kskip){
        /* Ar[ii*ndata+i]=tmp1; */
        Ar[ii][i]=tmp1;
      }
      /* Ap[ii*ndata+i]=tmp2; */
      Ap[ii][i]=tmp2;
    }
  }
}
/* void DoubleCalDeltaEtaZetaKCG(double *delta, double *eta, double *zeta, double *Ar, double *Ap, double *rvec, double *pvec, int ndata, int kskip) */
/* void DoubleCalDeltaEtaZetaKCG(double *delta, double *eta, double *zeta, double *Ar, double **Ap, double *rvec, double *pvec, int ndata, int kskip) */
void DoubleCalDeltaEtaZetaKCG(double *delta, double *eta, double *zeta, double **Ar, double **Ap, double *rvec, double *pvec, int ndata, int kskip)
{
  int i, j;
  double tmp1=0.0;
  double tmp2=0.0;
  double tmp3=0.0;
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2, tmp3) schedule(static) firstprivate(delta, eta, zeta, Ar, rvec, Ap, pvec) lastprivate(delta, eta, zeta)
  for(i=0;i<2*kskip+2;i++){
    tmp1=0.0;
    tmp2=0.0;
    tmp3=0.0;
    for(j=0;j<ndata;j++){
      if(i<2*kskip){
        /* tmp1 += rvec[j]*Ar[i*ndata+j]; */
        tmp1 += rvec[j]*Ar[i][j];
      }
      if(i<2*kskip+1){
        /* tmp2 += rvec[j]*Ap[i*ndata+j]; */
        tmp2 += rvec[j]*Ap[i][j];
      }
      /* tmp3 += pvec[j]*Ap[i*ndata+j]; */
      tmp3 += pvec[j]*Ap[i][j];
    }
    if(i<2*kskip){
      delta[i]=tmp1;
    }
    if(i<2*kskip+1){
      eta[i]=tmp2;
    }
    zeta[i]=tmp3;
  }
}
void DoubleCalArApKCR(double *Ar, double *Ap, double *val, int *col, int *ptr, double *rvec, double *pvec, int ndata, int kskip)
{
  int i, j, ii;
    double tmp1=0.0;
    double tmp2=0.0;

#pragma omp parallel for private(j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val, pvec, rvec) lastprivate(Ar, Ap)
  for(i=0;i<ndata;i++)
  {
    tmp1=0.0;
    tmp2=0.0;
    for(j=ptr[i];j<ptr[i+1];j++)
    {
      tmp1 += val[j] * rvec[col[j]];
      tmp2 += val[j] * pvec[col[j]];
    }
    Ar[i]=tmp1;
    Ap[i]=tmp2;
  }
 
#pragma omp parallel for private(i, j) reduction(+:tmp1, tmp2) schedule(static) firstprivate(Ar, Ap, val) lastprivate(Ar, Ap)
  for(ii=1;ii<2*kskip+1;ii++){
    for(i=0;i<ndata;i++){
      tmp1=0.0;
      tmp2=0.0;
      for(j=ptr[i];j<ptr[i+1];j++){
        tmp1 += val[j]*Ar[(ii-1)*ndata+col[j]];
        tmp2 += val[j]*Ap[(ii-1)*ndata+col[j]];
      }
      Ar[ii*ndata+i]=tmp1;
      Ap[ii*ndata+i]=tmp2;
    }
  }
 
}
void DoubleCalDeltaEtaZetaKCR(double *delta, double *eta, double *zeta, double *Ar, double *Ap, double *rvec, int ndata, int kskip)
{
  int i, j;
  double tmp1=0.0;
  double tmp2=0.0;
  double tmp3=0.0;
#pragma omp parallel for private(j) reduction(+:tmp1, tmp2, tmp3) schedule(static) firstprivate(delta, eta, zeta, Ar, rvec, Ap) lastprivate(delta, eta, zeta)
  for(i=0;i<2*kskip+1;i++){
    tmp1=0.0;
    tmp2=0.0;
    tmp3=0.0;
    for(j=0;j<ndata;j++){
      tmp1 += rvec[j]*Ar[i*ndata+j];
      tmp2 += Ap[0*ndata+j]*Ap[i*ndata+j];
      tmp3 += rvec[j]*Ap[i*ndata+j];
    }
    delta[i]=tmp1;
    eta[i]=tmp2;
    zeta[i]=tmp3;
  }
}

double DoubleCudaDot_Host(int N, double *a, double *b, int BlockPerGrid, int ThreadPerBlock){
  /* int i; */
  /* double sum=0.0; */
  /*  */
  /* double *tmp_h, *tmp_d; */
  /*  */
  /* tmp_h=(double *)malloc(sizeof(double)*BlockPerGrid); */
  /* checkCudaErrors( cudaMalloc((void **)&tmp_d, sizeof(double)*BlockPerGrid) ); */
  /*  */
  /* DoubleCudaDot<<<BlockPerGrid, ThreadPerBlock, sizeof(double)*ThreadPerBlock+16>>>(N, ThreadPerBlock, a, b, tmp_d); */
  /*  */
  /* checkCudaErrors( cudaMemcpy(tmp_h, tmp_d, sizeof(double)*BlockPerGrid, cudaMemcpyDeviceToHost) ); */
  /*  */
  /*  */
  /* for(i=0;i<BlockPerGrid;i++){ */
  /*   sum+=tmp_h[i]; */
  /* } */
  /*  */
  /* free(tmp_h); */
  /* checkCudaErrors( cudaFree(tmp_d) ); */
  /*  */
  /* return sum; */

  double *h_out=0, *d_out=0, sum=0.0;
  int i;
  checkCudaErrors(cudaMalloc((void **)&d_out, sizeof(double)*BlockPerGrid));
  h_out=(double *)malloc(sizeof(double)*BlockPerGrid);

  DoubleCudaDot<<<BlockPerGrid, ThreadPerBlock>>>(N, ThreadPerBlock, a, b, d_out);

  checkCudaErrors(cudaMemcpy(h_out, d_out, sizeof(double)*BlockPerGrid, cudaMemcpyDeviceToHost));

  for(i=0;i<BlockPerGrid;i++){
    sum+=h_out[i];
  }

  free(h_out);
  cudaFree(d_out);
  return sum;

}

