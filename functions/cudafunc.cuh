#ifndef CUDAFUNC_H_INCLUDED__
#define CUDAFUNC_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include "../start.h"

extern void 
DoubleCudaMalloc(double *ptr, 
    int size);
extern void 
IntgerCudaMalloc(int *ptr, 
    int size);
extern void DoubleMemCpyH2D(double *host, 
    double *device, 
    int size);
extern void IntgerMemCpyH2D(int *host, 
    int *device, 
    int size);
extern void DoubleMemCpyD2H(double *host, 
    double *device, 
    int size);
extern void IntgerMemCpyD2H(int *host, 
    int *device, 
    int size);
extern __global__ void
DoubleCudaMVMCSR(int n, double *val, int *col, int *ptr, double *b, double *c);

extern __global__ void
DoubleCudaMVMCSR2(int n, double *val, int *col, int *ptr, double *b, double *c);

extern __global__ void
DoubleCudaDot(int n, int ThreadPerBlock, double *a, double *b, double *c);

extern __device__ double 
partial_dot( const double* v1, const double* v2, int N, double* out  );

extern __device__ double
sum( const double* v );

extern __global__ void
full_dot( const double* v1, const double* v2, int N, double* out );




#endif //CUDAFUNC_H_INCLUDED__

