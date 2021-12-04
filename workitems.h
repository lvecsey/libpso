#ifndef WORKITEMS_H
#define WORKITEMS_H

#include "pso.h"

typedef struct {

  workunit *rptr_work;
  workunit *wptr_work;
  workunit *eptr_work;
  
} workptrs;

long int available_workitems(workunit *works, long int num_works, workunit *start_ptr, workunit *end_ptr);

#endif
