#include "start.h"
#include "io.h"
#include "blas.h"
#include "solvers.h"


int CSR_start(int argc, char const* argv[]){
  int N, NNZ;
  char basepath[512]="";
  char bxpath[512]="";
  char colpath[512]="";
  char ptrpath[512]="";

  double *bvec,*xvec, *val;
  int *col, *ptr;
  
  omp_set_num_threads(8);
  
  int error = UsageCheck(argc, argv);
  if(error!=0){
    printf("error in start\n");
    return -1;
  }

  strcat(basepath, "./Matrix/CSR/");
  strcat(basepath, argv[1]);
  strcat(basepath, "/");
  strcat(bxpath,basepath);
  strcat(colpath,basepath);
  strcat(ptrpath,basepath);
  strcat(bxpath, "bx.txt");
  strcat(colpath, "ColVal.txt");
  strcat(ptrpath, "Ptr.txt");

  GetHead(colpath, ptrpath, bxpath, &N, &NNZ);
  
  bvec=Double1Malloc(N);
  xvec=Double1Malloc(N);

  val=Double1Malloc(NNZ);
  col=Intger1Malloc(NNZ);
  ptr=Intger1Malloc(N+1);

  GetData(colpath, ptrpath, bxpath, col, ptr, val, bvec, xvec, N, NNZ);


  error = SolverSelecter(val, col, ptr, bvec, xvec, N, EPS, I_MAX, KSKIP, FIX);
  if(error!=0){
    printf("error in start\n");
    return(-1);
  }

  Double1Free(bvec);
  Double1Free(xvec);
  Double1Free(val);
  Intger1Free(col);
  Intger1Free(ptr);

  return 0;
}
