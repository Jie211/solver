#ifndef IO_H_INCLUDED__
#define IO_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>
#include "start.h"

extern int 
UsageCheck(int argc, 
    char const*argv[]);

extern int 
FileFound(int argc, 
    char const* argv[]);

extern void
GetHead(const char *bx, 
    const char *col, 
    const char *ptr, 
    int *n, 
    int *nnz);
extern void 
GetData(const char *file1, 
    const char *file2, 
    const char *file3, 
    int *col, 
    int *ptr, 
    double *val, 
    double *b, 
    double *x, 
    int N, 
    int NZ);

extern FILE* 
FileInit(char *name, 
    char *mode);

extern void 
FileClose(FILE *fp);

extern void 
FileOutPutVec(FILE *fp, 
    double *vec, 
    int ndata);
#endif //IO_H_INCLUDED__

