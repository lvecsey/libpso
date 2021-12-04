#ifndef FITNESSPACK_H
#define FITNESSPACK_H

#include "mini_gxkit.h"

typedef struct {

  double (*fitness_func)(point3d_t *pnta, void *ff_extra);
  
  void *ff_extra;
  
} fitnesspack;

#endif
