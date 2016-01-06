#include "cudafunc.cuh"

void DoubleCudaMalloc(double *ptr, int size){
  cudaMalloc((void **)&ptr, sizeof(double)*size);
}
void IntgerCudaMalloc(int *ptr, int size){
  cudaMalloc((void **)&ptr, sizeof(int)*size);
}
void DoubleMemCpyH2D(double *host, double *device, int size){
  cudaMemcpy(device, host, sizeof(double)*size, cudaMemcpyHostToDevice);
}
void IntgerMemCpyH2D(int *host, int *device, int size){
  cudaMemcpy(device, host, sizeof(int)*size, cudaMemcpyHostToDevice);
}
void DoubleMemCpyD2H(double *host, double *device, int size){
  cudaMemcpy(host, device, sizeof(double)*size, cudaMemcpyDeviceToHost);
}
void IntgerMemCpyD2H(int *host, int *device, int size){
  cudaMemcpy(host, device, sizeof(int)*size, cudaMemcpyDeviceToHost);
}



__global__ void
DoubleCudaMVMCSR(int n, double *val, int *col, int *ptr, double *b, double *c){
  extern __shared__ double vals[];

  int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int warp_id = thread_id/32;
  int lane = thread_id & (32 - 1);

  int row = warp_id;
  if(row<n)
  {
    int row_start = ptr[row];
    int row_end = ptr[row+1];

    vals[threadIdx.x] = 0.0;

    for(int jj = row_start+lane; jj<row_end; jj+=32)
    { 
      vals[threadIdx.x]+=val[jj] * b[col[jj]];
    }

    if(lane <16)
      vals[threadIdx.x] += vals[threadIdx.x +16];
    if(lane<8)
      vals[threadIdx.x] += vals[threadIdx.x + 8];
    if(lane<4)
      vals[threadIdx.x] += vals[threadIdx.x + 4];
    if(lane<2)
      vals[threadIdx.x] += vals[threadIdx.x + 2];
    if(lane<1)
      vals[threadIdx.x] += vals[threadIdx.x + 1];

    if(lane == 0){
      c[row] += vals[threadIdx.x];
    }
  }
}
__global__ void
DoubleCudaMVMCSR2(int n, double *val, int *col, int *ptr, double *b, double *c)
{
  long row=blockDim.x * blockIdx.x + threadIdx.x;
  long int i;
  if(row<n){
    double tmp=0.0;
    long int row_start=ptr[row];
    long int row_end=ptr[row+1];
    for(i=row_start;i<row_end;i++){
      tmp+=val[i]*b[col[i]];
    }
    /* __syncthreads(); */
    c[row]=tmp;
    /* printf("%d %.12e\n", row, c[row]); */
  }
  /* __syncthreads(); */
}

__global__ void
DoubleCudaDot(int n, int ThreadPerBlock, double *a, double *b, double *c) {
  extern __shared__ double share[];

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j;

    if(i < n)
      share[threadIdx.x] = a[i] * b[i];
    else
      share[threadIdx.x] = 0.0;
    __syncthreads();

    for(j=ThreadPerBlock/2; j>31; j>>=1){
      if(threadIdx.x < j)
        share[threadIdx.x] += share[threadIdx.x + j];
      __syncthreads();
    }
    if(threadIdx.x < 16){
      share[threadIdx.x] += share[threadIdx.x + 16];
      __syncthreads();
      share[threadIdx.x] += share[threadIdx.x + 8];
      __syncthreads();
      share[threadIdx.x] += share[threadIdx.x + 4];
      __syncthreads();
      share[threadIdx.x] += share[threadIdx.x + 2];
      __syncthreads();
      share[threadIdx.x] += share[threadIdx.x + 1];
    }
    __syncthreads();

    if(threadIdx.x == 0)
      c[blockIdx.x] = share[0];
}
