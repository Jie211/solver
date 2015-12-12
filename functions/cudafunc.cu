#include "cudafunc.h"

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
