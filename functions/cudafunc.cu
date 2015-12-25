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
DoubleCudaDot(double *out, double *x, double *y, int ndata){
  extern __shared__ double tmp[];
  int t = blockDim.x * blockIdx.x + threadIdx.x;
  int loc_t = threadIdx.x;
  int block_sz = blockDim.x;

  if (t < ndata) tmp[loc_t] = x[t]*y[t];
  __syncthreads();

  for (int stride = blockDim.x/2; stride > 32; stride /= 2) {
    if (loc_t < stride)
      tmp[loc_t] += tmp[loc_t + stride];
    __syncthreads();
  }

  if (loc_t < 32) {
    if (block_sz >= 64) tmp[loc_t] += tmp[loc_t + 32];
    if (block_sz >= 32) tmp[loc_t] += tmp[loc_t + 16];
    if (block_sz >= 16) tmp[loc_t] += tmp[loc_t + 8];
    if (block_sz >= 8) tmp[loc_t] += tmp[loc_t + 4];
    if (block_sz >= 4) tmp[loc_t] += tmp[loc_t + 2];
    if (block_sz >= 2) tmp[loc_t] += tmp[loc_t + 1];
  }

  if (threadIdx.x == 0) {
    out[blockIdx.x] = tmp[0];
  }
}
