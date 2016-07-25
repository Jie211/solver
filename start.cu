#include "start.h"

char *c_matrix=NULL;
char *c_solver_outer=NULL;
char *c_solver_inner=NULL;
char *c_loop_outer=NULL;
char *c_loop_inner=NULL;
char *c_eps_outer=NULL;
char *c_eps_inner=NULL;
char *c_restart_outer=NULL;
char *c_restart_inner=NULL;
char *c_kskip_outer=NULL;
char *c_kskip_inner=NULL;
char *c_fix_inner=NULL;
char *c_fix_outer=NULL;
char *c_openmp_thread=NULL;
char *c_verbose=NULL;
char *c_cuda=NULL;


bool f_matrix=false;
bool f_solver_outer=false;
bool f_solver_inner=false;
bool f_loop_outer=false;
bool f_loop_inner=false;
bool f_eps_outer=false;
bool f_eps_inner=false;
bool f_restart_outer=false;
bool f_restart_inner=false;
bool f_kskip_outer=false;
bool f_kskip_inner=false;
bool f_fix_inner=false;
bool f_fix_outer=false;
bool f_openmp_thread=false;
bool f_verbose=false;
bool f_cuda=false;

bool S_CG=false;
bool S_CR=false;
bool S_GCR=false;
bool S_GMRES=false;
bool K_CG=false;
bool K_CR=false;
bool VP_CG=false;
bool VP_CR=false;
bool VP_GCR=false;
bool VP_GMRES=false;
bool IS_CG=false;
bool IS_CR=false;
bool IS_GCR=false;
bool IK_CG=false;
bool IK_CR=false;

bool INNER=false;

char *matrix=NULL;
char *solver_outer=NULL;
char *solver_inner=NULL;
int loop_outer=L_OUT;
int loop_inner=L_IN;
double eps_outer=E_OUT;
double eps_inner=E_IN;
int restart_outer=R_OUT;
int restart_inner=R_IN;
int kskip_outer=K_OUT;
int kskip_inner=K_IN;
int fix_outer=F_OUT;
int fix_inner=F_IN;
int openmp_thread=THREAD;
bool verbose=VERBOSE;
bool cuda=CUDA;

char bx_path[512];
char ptr_path[512];
char col_path[512];

double *d_val=NULL;
int *d_col=NULL, *d_ptr=NULL;

/* __device__ unsigned int count ; */
/* __shared__ double cache[16]; */


int CSR_start(int argc, char *argv[]){
  int N, NNZ;
 
  char setThreads[10];
  double *bvec,*xvec, *val;
  int *col, *ptr;
  int error;

  double ost1, oet1, ot1;
  double ost2, oet2, ot2;
  
  error=getCMD(argc, argv);
  if(error!=0){
    Display_Err("error in start");
    return -1;
  }
  /* sprintf(setThreads, "%d", THREADS); */
  sprintf(setThreads, "%d", openmp_thread);
  /* omp_set_num_threads(THREADS); */
  error=setenv("OMP_NUM_THREADS",setThreads,1);
  omp_set_num_threads(atoi(setThreads));
  if(error!=0){
    Display_Err("set omp error in start");
    return -1;
  }
  printf("---- OpenMP set to %d ----\n", openmp_thread);
  
  error = UsageCheck(matrix);
  if(error!=0){
    Display_Err("error in start");
    return -1;
  }
  
  ost1=gettimeofday_sec();

  GetHead(col_path, ptr_path, bx_path, &N, &NNZ);
  
  /* bvec=Double1Malloc(N); */
  /* xvec=Double1Malloc(N); */
  /*  */
  /* val=Double1Malloc(NNZ); */
  /* col=Intger1Malloc(NNZ); */
  /* ptr=Intger1Malloc(N+1); */

  checkCudaErrors(cudaMallocHost((void **)&bvec, sizeof(double)*N));
  checkCudaErrors(cudaMallocHost((void **)&xvec, sizeof(double)*N));

  checkCudaErrors(cudaMallocHost((void **)&val, sizeof(double)*NNZ));
  checkCudaErrors(cudaMallocHost((void **)&col, sizeof(int)*NNZ));
  checkCudaErrors(cudaMallocHost((void **)&ptr, sizeof(int)*(N+1)));



  GetData(col_path, ptr_path, bx_path, col, ptr, val, bvec, xvec, N, NNZ);

  ost2=gettimeofday_sec();
 
  checkCudaErrors( cudaMalloc((void **)&d_val, sizeof(double)*NNZ) );
  checkCudaErrors( cudaMalloc((void **)&d_col, sizeof(int)*NNZ) );
  checkCudaErrors( cudaMalloc((void **)&d_ptr, sizeof(int)*(N+1)) );
  checkCudaErrors( cudaMemcpy(d_val, val, sizeof(double)*NNZ, cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpy(d_col, col, sizeof(int)*NNZ, cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpy(d_ptr, ptr, sizeof(int)*(N+1), cudaMemcpyHostToDevice) );
  
  oet2=gettimeofday_sec();

  error = SolverSelecter(val, col, ptr, bvec, xvec, N, NNZ, eps_outer, loop_outer, kskip_outer, fix_outer);


  checkCudaErrors( cudaFree(d_val) );
  checkCudaErrors( cudaFree(d_col) );
  checkCudaErrors( cudaFree(d_ptr) );

  oet1=gettimeofday_sec();

  ot1=oet1-ost1;
  ot2=oet2-ost2;
  if(cuda){
    printf("@@@@@@@ all=%.3f, mainCopy=%.3f\n", ot1, ot2);
  }else{
    printf("@@@@@@@ all=%.3f\n", ot1);
  }
  if(error!=0){
    Display_Err("error in start");
    return(-1);
  }

  if(verbose){
    DisplayCMD();
  }

  /* Double1Free(bvec); */
  /* Double1Free(xvec); */
  /* Double1Free(val); */
  /* Intger1Free(col); */
  /* Intger1Free(ptr); */

  checkCudaErrors(cudaFreeHost(bvec));
  checkCudaErrors(cudaFreeHost(xvec));

  checkCudaErrors(cudaFreeHost(val));
  checkCudaErrors(cudaFreeHost(col));
  checkCudaErrors(cudaFreeHost(ptr));

  return 0;
}

void InputCMD(void){
  if(solver_outer==NULL){
  }else if( strcmp(solver_outer, "cg") == 0 ){
    S_CG=true;
  }else if(strcmp(solver_outer, "cr") == 0){
    S_CR=true;
  }else if(strcmp(solver_outer, "gcr") == 0){
    S_GCR=true;
  }else if(strcmp(solver_outer, "gmres") == 0){
    S_GMRES=true;
  }else if(strcmp(solver_outer, "kcg") == 0){
    K_CG=true;
  }else if(strcmp(solver_outer, "kcr") == 0){
    K_CR=true;
  }else if(strcmp(solver_outer, "vpcg") == 0){
    VP_CG=true;
  }else if(strcmp(solver_outer, "vpcr") == 0){
    VP_CR=true;
  }else if(strcmp(solver_outer, "vpgcr") == 0){
    VP_GCR=true;
  }else if(strcmp(solver_outer, "vpgmres") == 0){
    VP_GMRES=true;
  }
  if(solver_inner==NULL){
  }else if(strcmp(solver_inner, "cg") == 0){
    IS_CG=true;
  }else if(strcmp(solver_inner, "cr") == 0){
    IS_CR=true;
  }else if(strcmp(solver_inner, "gcr") == 0){
    IS_GCR=true;
  }else if(strcmp(solver_inner, "kcg") == 0){
    IK_CG=true;
  }else if(strcmp(solver_inner, "kcr") == 0){
    IK_CR=true;
  }
}
int CheckCMD(void){
  if(!f_matrix){
    Display_Mes("Must set Matrix name");
    return -1;
  }else if(!f_solver_outer){
    Display_Mes("Must set a Solver");
    return -1;
  }else if(VP_CG || VP_CR || VP_GCR || VP_GMRES){
    INNER=true;
  }else if( (VP_CG || VP_CR || VP_GCR || VP_GMRES)  && ((!IS_CG) && (!IS_CR) && (!IS_GCR) && (!IK_CG) && (!IK_CR) )){
    Display_Mes("If VP method is selected, Please select inner method. [-OuterSolver=]");
    return -1;
  }
  return 0;
}
void DisplayCMD(void){
  printf("*******************************************\n");
  printf(" Matrix: %s\n", matrix);
  printf(" OpenMPThreads: %d\n", openmp_thread);
  printf(" Verbose: %d\n", verbose);
  printf(" CUDA: %d\n", cuda);
  printf("*******************************************\n");
  printf(" OuterSolver: %s\n", solver_outer);
  printf(" OuterLoop: %d\n", loop_outer);
  printf(" OuterEPS: %.12e\n", eps_outer);
  printf(" OuterRestart: %d\n", restart_outer);
  printf(" OuterKskip: %d\n", kskip_outer);
  printf(" OuterFix: %d\n", fix_outer);
  if(f_solver_inner){
    printf("*******************************************\n");
    printf(" InnerSolver: %s\n", solver_inner);
    printf(" InnerLoop: %d\n", loop_inner);
    printf(" InnerEPS: %.12e\n", eps_inner);
    printf(" InnerRestart: %d\n", restart_inner);
    printf(" InnerKskip: %d\n", kskip_inner);
    printf(" InnerFix: %d\n", fix_inner);
  }
  printf("*******************************************\n");
}
int getCMD(int argc, char *argv[])
{
  /* int i; */

  if(argc==1){
<<<<<<< HEAD
    printf("Option: \n"
        "\tRequire:\n"
        "\t\t-M/--Matrix [matrix name]-> Matrix to solve\n"
        "\t\t\t-S/--OuterSolver [solver name]-> Select Outsider solver\n"
        "\t\t\t-s/--InnerSolver [solver name]-> Select Insider solver\n"
        "\tOther:\n"
        "\t\t-L/--OuterLoop [int]-> maxloop for Outer solver\n"
        "\t\t-l/--InnerLoop [int]-> maxloop for Inner solver\n"
        "\t\t-E/--OuterEPS [double]-> EPS for Outer solver\n"
        "\t\t-e/--InnerEPS [double]-> EPS for Inner solver\n"
        "\t\t-R/--OuterRestart [int]-> Restart counter for Outer solver\n"
        "\t\t-r/--InnerRestart [int]-> Restart counter for Outer solver\n"
        "\t\t-K/--OuterKskip [int]-> skip num for K-skip Outer solver\n"
        "\t\t-k/--InnerKskip [int]-> skip num for K-skip Inner solver\n"
        "\t\t-F/--OuterFix [1,2]-> BugFix for K-skip Outer solver(DEBUG)\n"
        "\t\t-f/--InnerFix [1,2]-> BugFix for K-skip Inner solver(DEBUG)\n"
        "\t\t-V/--Verbose [0,1]-> verbose mode\n"
        "\t\t-T/--Thread [num]-> Thread for OpenMP\n"
        "\t\t-C/--Cuda [0,1]-> Cuda mode\n");
    /* printf("Option: Matrix, OuterSolver, InnerSolver, OuterLoop, InnerLoop, OuterEPS, InnerEPS, OuterRestart, InnerRestart, OuterKskip, InnerKskip, OuterFix, InnerFix, Verbose, Cuda\n"); */
=======
    printf("Option: -Matrix=[matrix], -Verbose=[=0,1], -Cuda=[=0,1]\n");
    printf("-OuterSolver=[method], -OuterLoop=[loops], -OuterEPS=[eps], -OuterRestart[RestartTimes], -OuterKskip=[k], -OuterFix=[debug]\n");
    printf("-InnerSolver=[method], -InnerLoop=[loops], -InnerESP=[eps], -InnerRestart[RestartTimes], -InnerKskip=[k], -InnerFix=[debug]\n");
    printf("Method: cg, cr, gcr, gmres, kcg, kcr, vpcg, vpcr, vpgmres\n");
>>>>>>> 7d5c6de0d237c69d3aad1f15f115ac970c45ebdc
    return -1;
  }

  struct option longopts[] = {
    {"Matrix", required_argument, NULL, 'M'},
    {"OuterSolver", required_argument, NULL, 'S'},
    {"InnerSolver", optional_argument, NULL, 's'},
    {"OuterLoop", optional_argument, NULL, 'L'},
    {"InnerLoop", optional_argument, NULL, 'l'},
    {"OuterEPS", optional_argument, NULL, 'E'},
    {"InnerEPS", optional_argument, NULL, 'e'},
    {"OuterRestart", optional_argument, NULL, 'R'},
    {"InnerRestart", optional_argument, NULL, 'r'},
    {"OuterKskip", optional_argument, NULL, 'K'},
    {"InnerKskip", optional_argument, NULL, 'k'},
    {"OuterFix", optional_argument, NULL, 'F'},
    {"InnerFix", optional_argument, NULL, 'f'},
    {"Thread", optional_argument, NULL, 'T'},
    {"Verbose", optional_argument, NULL, 'V'},
    {"Cuda", optional_argument, NULL, 'C'},
    { 0,        0,                 0,     0  },
  };

  int opt;
  int longindex;
  while((opt=getopt_long_only(argc, argv, "M:S:s:L:l:E:e:R:r:K:k:F:f:T:V:C::", longopts, &longindex)) != -1){
    switch(opt){
      case 'M':
        f_matrix=true;
        c_matrix=optarg;
        matrix=optarg;
        break;
      case 'S':
        f_solver_outer=true;
        c_solver_outer=optarg;
        solver_outer=optarg;
        break;
      case 's':
        f_solver_inner=true;
        c_solver_inner=optarg;
        solver_inner=optarg;
        break;
      case 'L':
        f_loop_outer=true;
        c_loop_outer=optarg;
        loop_outer=atoi(c_loop_outer);
        break;
      case 'l':
        f_loop_inner=true;
        c_loop_inner=optarg;
        loop_inner=atoi(c_loop_inner);
        break;
      case 'E':
        f_eps_outer=true;
        c_eps_outer=optarg;
        eps_outer=atof(c_eps_outer);
        break;
      case 'e':
        f_eps_inner=true;
        c_eps_inner=optarg;
        eps_inner=atof(c_eps_inner);
        break;
      case 'R':
        f_restart_outer=true;
        c_restart_outer=optarg;
        restart_outer=atoi(c_restart_outer);
        break;
      case 'r':
        f_restart_inner=true;
        c_restart_inner=optarg;
        restart_inner=atoi(c_restart_inner);
        break;
      case 'K':
        f_kskip_outer=true;
        c_kskip_outer=optarg;
        kskip_outer=atoi(c_kskip_outer);
        break;
      case 'k':
        f_kskip_inner=true;
        c_kskip_inner=optarg;
        kskip_inner=atoi(c_kskip_inner);
        break;
      case 'F':
        f_fix_outer=true;
        c_fix_outer=optarg;
        fix_outer=atoi(c_fix_outer);
        break;
      case 'f':
        f_fix_inner=true;
        c_fix_inner=optarg;
        fix_inner=atoi(c_fix_inner);
        break;
      case 'T':
        f_openmp_thread=true;
        c_openmp_thread=optarg;
        openmp_thread=atoi(c_openmp_thread);
        break;
      case 'V':
        f_verbose=true;
        c_verbose=optarg;
        verbose=atoi(c_verbose);
        /* if(verbose!=true || verbose!=false){ */
        /*   verbose=0; */
        /*   Display_Err("verbose must set to 0 or 1"); */
        /* } */
        break;
      case 'C':
        f_cuda=true;
        c_cuda=optarg;
        cuda=atoi(c_cuda);
        break;

      default:
        printf("error \'%c\' \'%c\'\n", opt, optopt);
        return 1;
    }
  }

  InputCMD();
  if(CheckCMD()!=0)
    return -1;
  DisplayCMD();

  return 0;
}

