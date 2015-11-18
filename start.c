#include "start.h"
#include "io.h"
#include "blas.h"
#include "solvers.h"

int DisplaySolver(void){

#if (defined (VP_CG) || defined (VP_CR) || defined (VP_GCR))  &&  (!defined (IS_CG) && !defined (IS_CR) && !defined (IK_CG) && !defined (IK_CR) && !defined(IS_GCR))
  Display_Err("if use VP solver, you need to select some Inner solver else")
  return -1;
#endif

#ifdef S_CG
  Display_Mes("CG selected");
#elif S_CR
  Display_Mes("CG selected");
#elif K_CG
  Display_Mes("Kskip-CG selected");
#elif K_CR
  Display_Mes("Kskip-CR selected");
#elif VP_CG
  Display_Mes("VPCG selected");
#elif VP_CR
  Display_Mes("VPCR selected");
#elif VP_GCR
  Display_Mes("VPGCR selected");
#elif S_GCR
  Display_Mes("GCR selected");
#else
  Display_Err("no solver selected");
  return -1;
#endif

#ifdef IS_CG
  Display_Mes("Inner CG selected");
#elif IS_CR
  Display_Mes("Inner CR selected");
#elif IK_CG
  Display_Mes("Inner Kskip-CG selected");
#elif IK_CR
  Display_Mes("Inner Kskip-CR selected");
#elif IS_GCR
  Display_Mes("Inner GCR selected");
#else
  Display_Err("no inner solver selected");
  return -1;
#endif

  return 0;
}

int CSR_start(int argc, char const* argv[]){
  int N, NNZ;
 
  double *bvec,*xvec, *val;
  int *col, *ptr;
  int error;
  
  error=DisplaySolver();
  if(error!=0){
    Display_Err("error in start");
    return -1;
  }

  omp_set_num_threads(THREADS);
  printf("---- OpenMP set to %d ----\n", THREADS);
  
  error = UsageCheck(argc, argv);
  if(error!=0){
    Display_Err("error in start");
    return -1;
  }
  
  GetHead(col_path, ptr_path, bx_path, &N, &NNZ);
  
  bvec=Double1Malloc(N);
  xvec=Double1Malloc(N);

  val=Double1Malloc(NNZ);
  col=Intger1Malloc(NNZ);
  ptr=Intger1Malloc(N+1);

  GetData(col_path, ptr_path, bx_path, col, ptr, val, bvec, xvec, N, NNZ);


  error = SolverSelecter(val, col, ptr, bvec, xvec, N, EPS, I_MAX, KSKIP, FIX);

  if(error!=0){
    Display_Err("error in start");
    return(-1);
  }

  Double1Free(bvec);
  Double1Free(xvec);
  Double1Free(val);
  Intger1Free(col);
  Intger1Free(ptr);

  return 0;
}
