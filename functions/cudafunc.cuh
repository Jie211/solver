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
DoubleCudaDot(double *out, double *x, double *y, int ndata);

#endif //CUDAFUNC_H_INCLUDED__

