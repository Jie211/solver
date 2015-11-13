#include "io.h"
#include <sys/stat.h>

int UsageCheck(int argc, char const* argv[]){
  int error = FileFound(argc, argv);
  if(error!=0){
    printf("filefound error\n");
    return -1;
  }

  return 0;
}
int FileFound(int argc, char const* argv[]){
  DIR *dir;
  struct dirent *dp;
  struct stat st;
  char path[512]="";
  char fullpath[512+512]="";
  char searchname[512]="";
  /* int result; */
  bool dir_found=false;
  bool bx=false;
  bool col=false;
  bool ptr=false;
  bool file_found=false;

  if(argc!=2){
    printf("./a.out [matrix file name]\n");
    return -1;
  }else{
    strcpy(searchname, argv[1]);
    strcpy(path, "../Matrix/CSR/");
  }

  if((dir=opendir(path))==NULL){
    perror("opendir");
    return -1;
  }
  for(dp=readdir(dir);dp!=NULL;dp=readdir(dir)){
    /* result = stat(dp->d_name,&st); */
    stat(dp->d_name,&st);
    /* if((st.st_mode & S_IFMT) == S_IFDIR){ */
    if(S_ISDIR(st.st_mode)){
      if(strcmp(dp->d_name, searchname) == 0){
        dir_found=true;
        break;
      }
    }
  }
  if(!dir_found){
    printf("Matrix directory not found\n");
  }else{
    printf("dirfound\n");
  }
  
  strcpy(fullpath, path);
  strcat(fullpath, searchname);
  strcat(fullpath, "/");
  if((dir=opendir(fullpath))==NULL){
    perror("opendir");
    return -1;
  }
  for(dp=readdir(dir);dp!=NULL;dp=readdir(dir)){
    /* result = stat(dp->d_name,&st); */
    stat(dp->d_name,&st);
    /* if((st.st_mode & S_IFMT) == S_IFDIR){ */
    if(S_ISDIR(st.st_mode)){
      if(strcmp(dp->d_name, "bx.txt") == 0){
        bx=true;
      }else if(strcmp(dp->d_name, "ColVal.txt") == 0){
        col=true;
      }else if(strcmp(dp->d_name, "Ptr.txt") == 0){
        ptr=true;
      }
      if(bx && col && ptr){
        file_found=true;
        break;
      }
    }
  }
  printf("%s\n",fullpath);
  if(!file_found){
    printf("Matrix file not found\n");
  }else{
    printf("filefound\n");
  }
  closedir(dir);
  return 0;
}

void GetHead(const char *bx, const char *col, const char *ptr, int *n, int *nnz)
{
  FILE *in1, *in2, *in3;

  if((in1 = fopen(bx, "r")) == NULL)
  {
    printf("head %s file open error\n", bx);
    exit(-1);
  }

  if((in2 = fopen(col, "r")) == NULL)
  {
    printf("head %s file open error\n", col);
    exit(-1);
  }

  if((in3 = fopen(ptr, "r")) == NULL)
  {
    printf("head %s file open error\n", ptr);
    exit(-1);
  }
  int N11, N12, N21, N22, N31, N32;
  int NZ1, NZ2, NZ3;

  fscanf(in1, "%d %d %d\n", &N11, &N12, &NZ1);
  fscanf(in2, "%d %d %d\n", &N21, &N22, &NZ2);
  fscanf(in3, "%d %d %d\n", &N31, &N32, &NZ3);

  if(N11!=N12)
  {
    printf("in %s N!=M \n", bx);
    exit(-1);
  }
  if(N21!=N22)
  {
    printf("in %s N!=M \n", col);
    exit(-1);
  }
  if(N31!=N32)
  {
    printf("in %s N!=M \n", ptr);
    exit(-1);
  }

  if(N11 != N21 || N21!=N31 || N31!=N11)
  {
    printf("N was not same in 3files\n");
    exit(-1);
  }

  if(NZ1 != NZ2 || NZ2!=NZ3 || NZ3!=NZ1)
  {
    printf("NNZ was not same in 3files\n");
    exit(-1);
  }
  *n = N11;
  *nnz = NZ1;

  fclose(in1);
  fclose(in2);
  fclose(in3);
}
void GetData(const char *file1, const char *file2, const char *file3, int *col, int *ptr, double *val, double *b, double *x, int N, int NZ)
{
  FILE *in1,*in2,*in3;
  if((in1 = fopen(file1, "r")) == NULL)
  {
    printf("%s file open error", file1);
    exit(0);
  }

  if((in2 = fopen(file2, "r")) == NULL)
  {
    printf("%s file open error", file2);
    exit(0);
  }

  if((in3 = fopen(file3, "r")) == NULL)
  {
    printf("%s file open error", file3);
    exit(0);
  }
  int getint;
  double getdouble, getdouble2;
  int skip1, skip2, skip3;

  fscanf(in1, "%d %d %d\n", &skip1, &skip2, &skip3);
  fscanf(in2, "%d %d %d\n", &skip1, &skip2, &skip3);
  fscanf(in3, "%d %d %d\n", &skip1, &skip2, &skip3);
  for(int i=0;i<NZ;i++)
  {
    fscanf(in1,"%d %le\n",&getint,&getdouble);
    col[i] = getint;
    val[i] = getdouble;
  }

  for(int i=0;i<N+1;i++)
  {
    fscanf(in2,"%d\n",&getint);
    ptr[i] = getint;
  }

  for(int i=0;i<N;i++)
  {
    fscanf(in3,"%le %le\n",&getdouble,&getdouble2);
    b[i] = getdouble;
    x[i] = getdouble2;
  }


  fclose(in1);
  fclose(in2);
  fclose(in3);
}
FILE* FileInit(char *name, char *mode){
  FILE *tmp;
  if((tmp = fopen(name, mode))==NULL){
    perror("File init error\n");
    exit(-1);
  }
  return (tmp);
}
void FileClose(FILE *fp){
  fclose(fp);
}
void FileOutPutVec(FILE *fp, double *vec, int ndata){
  int i;
  for(i=0;i<ndata;i++){
    fprintf(fp,"%d %.12e", i, vec[i]);
  }
}
