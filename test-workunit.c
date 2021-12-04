
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "pso.h"

int main(int argc, char *argv[]) {

  size_t sz;

  sz = sizeof(workunit);

  printf("Size of workunit %lu\n", sz);
  
  return 0;

}
