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
#endif //CUDAFUNC_H_INCLUDED__

