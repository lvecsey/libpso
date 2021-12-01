#ifndef PSO_UTIL_H
#define PSO_UTIL_H

#include "mini_gxkit.h"

double calc_single(point3d_t *cmp_pnt1, point3d_t *cmp_pnt2);

double calc_totalerror(point3d_t *particles, long int num_particles, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra);

double fetch_rnd(int rnd_fd);

int within_range(point3d_t *pnta, double blo, double bup);

#endif
