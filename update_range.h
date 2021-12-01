#ifndef UPDATE_RANGE_H
#define UPDATE_RANGE_H

#include <pthread.h>

#include "mini_gxkit.h"

#include "particlepack.h"

int update_rangepso(particlepack *pp, pthread_mutex_t *bestparticle_mutex, double percent, int rnd_fd, long int particle_start, long int particle_end, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra, double blo, double bup);

#endif
