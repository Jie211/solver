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
  /* extern __shared__ double share[]; */
  /*  */
  /*   int i = blockDim.x * blockIdx.x + threadIdx.x; */
  /*   int j; */
  /*  */
  /*   if(i < n) */
  /*     share[threadIdx.x] = a[i] * b[i]; */
  /*   else */
  /*     share[threadIdx.x] = 0.0; */
  /*   __syncthreads(); */
  /*  */
  /*   for(j=ThreadPerBlock/2; j>31; j>>=1){ */
  /*     if(threadIdx.x < j) */
  /*       share[threadIdx.x] += share[threadIdx.x + j]; */
  /*     __syncthreads(); */
  /*   } */
  /*   if(threadIdx.x < 16){ */
  /*     share[threadIdx.x] += share[threadIdx.x + 16]; */
  /*     __syncthreads(); */
  /*     share[threadIdx.x] += share[threadIdx.x + 8]; */
  /*     __syncthreads(); */
  /*     share[threadIdx.x] += share[threadIdx.x + 4]; */
  /*     __syncthreads(); */
  /*     share[threadIdx.x] += share[threadIdx.x + 2]; */
  /*     __syncthreads(); */
  /*     share[threadIdx.x] += share[threadIdx.x + 1]; */
  /*   } */
  /*   __syncthreads(); */
  /*  */
  /*   if(threadIdx.x == 0) */
  /*     c[blockIdx.x] = share[0]; */
  __shared__ double cache[128];
  int i=blockIdx.x*(blockDim.x*2)+threadIdx.x;
  int tid=threadIdx.x;
  double mySum;
  if(i<n){
    mySum=a[i]*b[i];
  }else{
    mySum=0.0;
  }
  if(i+blockDim.x<n){
    mySum+=a[i+blockDim.x]*b[i+blockDim.x];
  }
  cache[tid]=mySum;
  __syncthreads();
  for(int s=blockDim.x/2;s>32;s>>=1){
    if(tid<s){
      cache[tid]=mySum=mySum+cache[tid+s];
    }
    __syncthreads();
  }
  if(tid<32){
    if(blockDim.x>=64)
      mySum+=cache[tid+32];
    for(int offset=32/2;offset>0;offset/=2){
      mySum+=__shfl_down(mySum,offset);
    }
  }
  if(tid==0){
    c[blockIdx.x]=mySum;
  }
}
__device__ double partial_dot( const double* v1, const double* v2, int N, double* out  ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if( i >= N ) return double( 0 );
    cache[ threadIdx.x ] = 0.f;
      while( i < N ) {
        cache[ threadIdx.x ] += v1[ i ] * v2[ i ];
          i += gridDim.x * blockDim.x;
      }
  __syncthreads();
  /* i = BLOCK_SIZE / 2;\ */
  i = 16 / 2;
    while( i > 0 ) {
      if( threadIdx.x < i ) cache[ threadIdx.x ] += cache[ threadIdx.x + i ];
        __syncthreads();
          i /= 2; 
    }
  return cache[ 0 ];
}

__device__ double sum( const double* v ) {
  cache[ threadIdx.x ] = 0.f;
  int i = threadIdx.x;
  while( i < gridDim.x ) {
    cache[ threadIdx.x ] += v[ i ];
      i += blockDim.x;
  }
  __syncthreads();
  /* i = BLOCK_SIZE / 2;\ */
  i = 16 / 2;
    while( i > 0 ) {
      if( threadIdx.x < i ) cache[ threadIdx.x ] += cache[ threadIdx.x + i ];
        __syncthreads();
          i /= 2; 
    }
  return cache[ 0 ];
}

__global__ void full_dot( const double* v1, const double* v2, int N, double* out ) {
  __shared__ bool lastBlock;
  double r = partial_dot( v1, v2, N, out );
  if( threadIdx.x == 0 ) {
    out[ blockIdx.x ] = r;
    __threadfence();
    const unsigned int v = atomicInc( &count, gridDim.x );
    lastBlock = ( v == gridDim.x - 1 );
  }
    __syncthreads();

  if( lastBlock ) {
    r = sum( out );
    if( threadIdx.x == 0 ) {
      out[ 0 ] =  r;
      count = 0;
    }
  }   
}


