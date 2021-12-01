
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "vector.h"

int main(int argc, char *argv[]) {

  size_t sz;

  sz = sizeof(vec3d);

  printf("Size of vec3d %ld\n", sz);
  
  return 0;

}
