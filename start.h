#ifndef START_H_INCLUDED__
#define START_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include "./functions/blas.h"
#include "./functions/io.h"
#include "solvers.h"

// #define THREADS 8

#define S_IN "NO"
#define L_OUT 10000
#define L_IN 100
#define E_OUT 1e-8
#define E_IN 1e-1
#define R_OUT 1000
#define R_IN 10
#define K_OUT 2
#define K_IN 2
#define F_OUT 2
#define F_IN 2
#define THREAD 8

extern char *c_matrix;
extern char *c_solver_outer;
extern char *c_solver_inner;
extern char *c_loop_outer;
extern char *c_loop_inner;
extern char *c_eps_outer;
extern char *c_eps_inner;
extern char *c_restart_outer;
extern char *c_restart_inner;
extern char *c_kskip_outer;
extern char *c_kskip_inner;
extern char *c_fix_inner;
extern char *c_fix_outer;
extern char *c_openmp_thread;


extern bool f_matrix;
extern bool f_solver_outer;
extern bool f_solver_inner;
extern bool f_loop_outer;
extern bool f_loop_inner;
extern bool f_eps_outer;
extern bool f_eps_inner;
extern bool f_restart_outer;
extern bool f_restart_inner;
extern bool f_kskip_outer;
extern bool f_kskip_inner;
extern bool f_fix_inner;
extern bool f_fix_outer;
extern bool f_openmp_thread;

extern bool S_CG;
extern bool S_CR;
extern bool S_GCR;
extern bool K_CG;
extern bool K_CR;
extern bool VP_CG;
extern bool VP_CR;
extern bool VP_GCR;
extern bool IS_CG;
extern bool IS_CR;
extern bool IS_GCR;
extern bool IK_CG;
extern bool IK_CR;

extern bool INNER;

extern char *matrix;
extern char *solver_outer;
extern char *solver_inner;
extern int loop_outer;
extern int loop_inner;
extern double eps_outer;
extern double eps_inner;
extern int restart_outer;
extern int restart_inner;
extern int kskip_outer;
extern int kskip_inner;
extern int fix_outer;
extern int fix_inner;
extern int openmp_thread;

extern char bx_path[512];
extern char ptr_path[512];
extern char col_path[512];

extern int 
CSR_start(int argc, 
    char *argv[]);

extern int 
DisplaySolver(void);

extern int 
getCMD(int argc, 
    char *argv[]);

extern void 
InputCMD(void);

extern int 
CheckCMD(void);

extern void 
DisplayCMD(void);

extern int 
getCMD(int argc, 
    char *argv[]);


#endif //START_H_INCLUDED__

