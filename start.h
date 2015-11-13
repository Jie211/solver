#include <stdio.h>
#include <omp.h>

#define EPS 1e-8
#define I_MAX 10000
#define KSKIP 2
#define FIX 2

#define I_EPS 1e-1
#define I_I_MAX 5
#define I_KSKIP 2
#define I_FIX 2

#if defined(VP_CG) || defined(VP_CR)
#define INNER
#endif

// #define BX_PATH ""
// #define PTR_PATH ""
// #define COL_PATH ""
//


extern 
int CSR_start(int argc, 
    char const* argv[]);
