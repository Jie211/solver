#include "start.h"
int main(int argc, char const* argv[])
{
  int error=CSR_start(argc, argv);
  if(error!=0){
    printf("error in main\n");
    return -1;
  }
  return 0;
}
