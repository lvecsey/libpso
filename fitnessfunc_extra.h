#ifndef FITNESSFUNC_EXTRA_H
#define FITNESSFUNC_EXTRA_H

typedef struct {

  double *inv_matrix;

  double *matrix;

  double blo, bup;
  
  double sf;
  
} fitnessfunc_extra;

#endif
