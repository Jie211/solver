#include <stdio.h>
#include <stdlib.h>
#include "./functions/blas.h"
#include "./functions/io.h"
#include "solvers.h"


#define THREADS 8

#define EPS 1e-8
#define I_MAX 8000

#define KSKIP 2
#define FIX 2

#define RESTART 1000
//########################//
#define I_EPS 1e-1
#define I_I_MAX 5

#define I_KSKIP 2
#define I_FIX 2

#define I_RESTART 1000

#if defined(VP_CG) || defined(VP_CR) || defined(VP_GCR)
#define INNER
#endif

char bx_path[512];
char ptr_path[512];
char col_path[512];

extern int 
CSR_start(int argc, 
    char const* argv[]);

extern int 
DisplaySolver(void);
