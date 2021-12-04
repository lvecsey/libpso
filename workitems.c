
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "pso.h"

#include "workitems.h"

long int available_workitems(workunit *works, long int num_works, workunit *start_ptr, workunit *end_ptr) {

  long int num_entries;

  if (start_ptr == end_ptr) {
    return 0;
  }
  
  if (start_ptr < end_ptr) {
      
    num_entries = (end_ptr - start_ptr);
	
  }

  else {

    num_entries = (num_works - (start_ptr - works));

    num_entries += (end_ptr - works);
	
  }

  return num_entries;

}  

