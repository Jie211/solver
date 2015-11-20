#include "start.h"
int main(int argc, char *argv[])
{
  int error=CSR_start(argc, argv);
  if(error!=0){
    Display_Err("error in main");
    return -1;
  }
  return 0;
}
