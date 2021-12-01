
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "update_range.h"

#include "pso_util.h"

int update_rangepso(particlepack *pp, pthread_mutex_t *bestparticle_mutex, double percent, int rnd_fd, long int particle_start, long int particle_end, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra, double blo, double bup) {

  point3d_t pnta, pntb;
  
  vec3d *vprev;
  point3d_t *xprev;
  
  vec3d *vcur;
  point3d_t *xcur;

  long int particleno;

  update_param *up;

  double inertia_weight;

  double c1, c2;

  double cur_err;

  point3d_t xnext;

  vec3d vnext;

  long int posno;
  
  vcur = pp->vcur;
  xcur = pp->xcur;

  up = &(pp->up);

  inertia_weight = (1.0 - percent);

  c1 = 2.0;
  c2 = 2.0;

  up->accel_c = 2.5;
  
  for (particleno = particle_start; particleno < particle_end; particleno++) {

    pnta = (point3d_t) { .x = pp->pbest[particleno].x - xcur[particleno].x, .y = pp->pbest[particleno].y - xcur[particleno].y, .z = pp->pbest[particleno].z - xcur[particleno].z };

    pthread_mutex_lock(bestparticle_mutex);
    pntb = (point3d_t) { .x = pp->gbest.x - xcur[particleno].x, .y = pp->gbest.y - xcur[particleno].y, .z = pp->gbest.z - xcur[particleno].z };
    pthread_mutex_unlock(bestparticle_mutex);

    pp->up.rand1 = fetch_rnd(rnd_fd);
    pp->up.rand2 = fetch_rnd(rnd_fd);
    vnext[0] = inertia_weight * vcur[particleno][0] + c1 * up->accel_c * up->rand1 * (pnta.x) + c2 * up->accel_c * up->rand2 * (pntb.x);

    pp->up.rand1 = fetch_rnd(rnd_fd);
    pp->up.rand2 = fetch_rnd(rnd_fd);
    vnext[1] = inertia_weight * vcur[particleno][1] + c1 * up->accel_c * up->rand1 * (pnta.y) + c2 * up->accel_c * up->rand2 * (pntb.y);

    pp->up.rand1 = fetch_rnd(rnd_fd);
    pp->up.rand2 = fetch_rnd(rnd_fd);
    vnext[2] = inertia_weight * vcur[particleno][2] + c1 * up->accel_c * up->rand1 * (pnta.z) + c2 * up->accel_c * up->rand2 * (pntb.z);

    xnext = (point3d_t) { .x = xcur[particleno].x + vnext[0], .y = xcur[particleno].y + vnext[1], .z = xcur[particleno].z + vnext[2] };
    
    if (within_range(&xnext, blo, bup)) {

      vcur[particleno][0] = vnext[0];
      vcur[particleno][1] = vnext[1];
      vcur[particleno][2] = vnext[2];    
      
      xcur[particleno] = xnext;
      
      cur_err = fitness_func(&xnext, ff_extra);

      if (cur_err < pp->fitness[particleno]) {

	pp->pbest[particleno] = xnext;
	pp->fitness[particleno] = cur_err;

	pthread_mutex_lock(bestparticle_mutex);	
	if (cur_err < pp->gbesterr) {
	  pp->gbest = xnext;
	  pp->gbesterr = cur_err;	
	}
	pthread_mutex_unlock(bestparticle_mutex);
	
      }
      
    }
	
  }
    
  return 0;

}
