#ifndef EXTRAPARAMS_H
#define EXTRAPARAMS_H

#include <time.h>

typedef struct {

  double totalerr_avg;
  
  double totalerr_prev;

  double totalerr;

  long int progressdisp_freq;
  
  struct timespec now, then;
  
} extraparams;

#endif
