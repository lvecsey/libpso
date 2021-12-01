
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <math.h>

#include "mini_gxkit.h"

#include "pso_util.h"

int show_point(point3d_t *pnta, char *desc) {

  fprintf(stderr, "%s: %g %g %g\n", desc, pnta->x, pnta->y, pnta->z);

  return 0;

}

double calc_single(point3d_t *cmp_pnt1, point3d_t *cmp_pnt2) {

  double diff[3];

  double err;

  diff[0] = (cmp_pnt1->x - cmp_pnt2->x);
  diff[1] = (cmp_pnt1->y - cmp_pnt2->y);
  diff[2] = (cmp_pnt1->z - cmp_pnt2->z);      

  err = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

  return err;
  
}

double calc_totalerror(point3d_t *particles, long int num_particles, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra) {

  double sum;
  
  double total_err;

  double err;
  
  long int particleno;
  
  sum = 0.0;

  for (particleno = 0; particleno < num_particles; particleno++) {

    err = fitness_func(particles + particleno, ff_extra);

    sum += err;
    
  }

  total_err = (sum / num_particles);
  
  return total_err;

}

double fetch_rnd(int rnd_fd) {

  uint64_t rndval;

  int retval;
  
  retval = read(rnd_fd, &rndval, sizeof(uint64_t));
  if (retval != sizeof(uint64_t)) {
    perror("fetch_rnd read");
    return -1;
  }

  return -1.0 + 2.0 * rndval / 18446744073709551615.0;

}

int within_range(point3d_t *pnta, double blo, double bup) {

  if (pnta->x > bup || pnta->y > bup || pnta->z > bup) {
    return 0;
  }

  if (pnta->x < blo || pnta->y < blo || pnta->z < blo) {
    return 0;
  }
  
  return 1;

}
