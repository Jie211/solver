#include "start.h"
#include "io.h"
#include "blas.h"
#include "solvers.h"


int CSR_start(int argc, char const* argv[]){
  int N, NNZ;
 
  double *bvec,*xvec, *val;
  int *col, *ptr;
  
  omp_set_num_threads(THREADS);
  printf("---- OpenMP set to %d ----\n", THREADS);
  
  int error = UsageCheck(argc, argv);
  if(error!=0){
    printf("** error in start **\n");
    return -1;
  }
  
  GetHead(col_path, ptr_path, bx_path, &N, &NNZ);
  
  bvec=Double1Malloc(N);
  xvec=Double1Malloc(N);

  val=Double1Malloc(NNZ);
  col=Intger1Malloc(NNZ);
  ptr=Intger1Malloc(N+1);

  GetData(col_path, ptr_path, bx_path, col, ptr, val, bvec, xvec, N, NNZ);


  error = SolverSelecter(val, col, ptr, bvec, xvec, N, EPS, I_MAX, KSKIP, FIX, RESTART);

  if(error!=0){
    printf("** error in start **\n");
    return(-1);
  }

  Double1Free(bvec);
  Double1Free(xvec);
  Double1Free(val);
  Intger1Free(col);
  Intger1Free(ptr);

  return 0;
}
