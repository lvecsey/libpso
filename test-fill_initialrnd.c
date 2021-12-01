
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "pso.h"

int main(int argc, char *argv[]) {

  int rnd_fd;

  rnd_fd = open("/dev/urandom", O_RDONLY);

  {

    particlepack pp;

    long int num_particles;
    
    int retval;

    num_particles = 1000;
    
    retval = alloc_ppack(&pp, num_particles);
    if (retval == -1) {
      printf("FAIL");
      return -1;
    }

    retval = fill_initialrnd(&pp, num_particles, rnd_fd, -5.0, 5.0);
    if (retval == -1) {
      printf("FAIL");
      return -1;
    }

  }

  printf("SUCCESS");

  return 0;

}
