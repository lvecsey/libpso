
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "bestpack.h"

#include "update_range.h"

#include "pso_util.h"

int update_rangepso(particlepack *pp, pthread_mutex_t *globalbest_mutex, double percent, int rnd_fd, long int particle_start, long int particle_end, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra, double blo, double bup, bestpack *tp_gbest) {

  point3d_t pnta, pntb;
  
  vec3d *vcur;
  point3d_t *xcur;

  long int particleno;

  update_param *up;

  double inertia_weight;

  double cur_err;

  point3d_t xnext;

  vec3d vnext;

  uint64_t rnds[6];
  
  double rands[6];

  ssize_t bytes_read;

  long int rndno;

  vcur = pp->vcur;
  xcur = pp->xcur;

  up = &(pp->up);

  inertia_weight = (1.0 - percent);

  for (particleno = particle_start; particleno < particle_end; particleno++) {
    
    pnta = (point3d_t) { .x = pp->pbest[particleno].x - xcur[particleno].x, .y = pp->pbest[particleno].y - xcur[particleno].y, .z = pp->pbest[particleno].z - xcur[particleno].z };
    pntb = (point3d_t) { .x = tp_gbest->best.x - xcur[particleno].x, .y = tp_gbest->best.y - xcur[particleno].y, .z = tp_gbest->best.z - xcur[particleno].z };

    bytes_read = read(rnd_fd, rnds, sizeof(uint64_t) * 6);
    if (bytes_read != sizeof(uint64_t) * 6) {
      perror("read");
      return -1;
    }

    for (rndno = 0; rndno < 6; rndno++) {
      rands[rndno] = rnds[rndno] / 18446744073709551615.0;
    }
    
    vnext[0] = inertia_weight * vcur[particleno][0] + up->c1 * up->accel_c * rands[0] * (pnta.x) + up->c2 * up->accel_c * rands[1] * (pntb.x);
    vnext[1] = inertia_weight * vcur[particleno][1] + up->c1 * up->accel_c * rands[2] * (pnta.y) + up->c2 * up->accel_c * rands[3] * (pntb.y);
    vnext[2] = inertia_weight * vcur[particleno][2] + up->c1 * up->accel_c * rands[4] * (pnta.z) + up->c2 * up->accel_c * rands[5] * (pntb.z);

    xnext = (point3d_t) { .x = xcur[particleno].x + vnext[0], .y = xcur[particleno].y + vnext[1], .z = xcur[particleno].z + vnext[2] };

    vcur[particleno][0] = vnext[0];
    vcur[particleno][1] = vnext[1];
    vcur[particleno][2] = vnext[2];    

    xcur[particleno] = xnext;

    cur_err = fitness_func(&xnext, ff_extra);

    if (cur_err < pp->fitness[particleno]) {

      pp->fitness[particleno] = cur_err;
      pp->pbest[particleno] = xnext;
	
      if (cur_err < tp_gbest->besterr) {
	tp_gbest->best = xnext;
	tp_gbest->besterr = cur_err;	
      }
	
    }
	
  }

  return 0;

}
