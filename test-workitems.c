
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "workitems.h"

int main(int argc, char *argv[]) {

  long int num_items;

  workunit *works;

  long int num_avail;

  long int num_collect;
  
  workunit *rptr_work;
  workunit *wptr_work;
  
  num_items = 84;

  num_collect = 28;
  
  works = malloc(sizeof(workunit) * num_items);

  rptr_work = works + 0;
  wptr_work = works + 0;
  
  num_avail = available_workitems(works, num_items, rptr_work, wptr_work);

  printf("%s: num_avail %ld\n", __FUNCTION__, num_avail);

  if (num_avail != 0) {
    printf("FAIL");
    return -1;
  }
  
  rptr_work = works + 5;
  wptr_work = works + 10;
  
  num_avail = available_workitems(works, num_items, rptr_work, wptr_work);

  printf("%s: num_avail %ld\n", __FUNCTION__, num_avail);

  if (num_avail != 5) {
    printf("FAIL");
    return -1;
  }

  rptr_work = works + 80;
  wptr_work = works + 5;
  
  num_avail = available_workitems(works, num_items, rptr_work, wptr_work);

  printf("%s: num_avail %ld\n", __FUNCTION__, num_avail);

  if (num_avail != 9) {
    printf("FAIL");
    return -1;
  }

  printf("SUCCESS");
  
  return 0;

}
