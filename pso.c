/*
    Use the PSO algorithm for a swarm of particles that matches a function
    Copyright (C) 2021  Lester Vecsey

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <stdint.h>

#include <pthread.h>

#include <math.h>

#include <time.h>

#include <errno.h>

#include <sys/queue.h>

#include <errno.h>

#include <string.h>

#include "mini_gxkit.h"

#include "extraparams.h"

#include "pso.h"

#include "pso_util.h"

#include "update_range.h"

#define debugprint fprintf

#define emptyprint(fp, str, ...) ;

#define dfprintf emptyprint

const long int debug_level = 1;

const long int def_numpart = 250000;

const long int def_numgenerations = 4000;

TAILQ_HEAD(tailhead, entry) head = TAILQ_HEAD_INITIALIZER(head);

int generate_work(pso *ps, long int num_workpacks, long int generationno);
workunit gen_workcalc(psorw *prw);
int show_threads(uint64_t *thread_states, long int num_threads);

int show_work(workunit *work);

double elap_cur(struct timespec *start, extraparams *ep);

int clear_wpstatus(workpack_markstatus *wp, long int num_packs);
int show_wpstatus(workpack_markstatus *wp, long int num_packs);

int init_pso(pso *ps, long int num_threads) {

  int retval;

  long int threadno;

  long int mutexno;
  
  ps->headp = &head;
  
  TAILQ_INIT(ps->headp);

  ps->mutex = malloc(sizeof(pthread_mutex_t) * num_mutex);
  if (ps->mutex == NULL) {
    perror("malloc");
    return -1;
  }
  
  for (mutexno = 0; mutexno < num_mutex; mutexno++) {
  
    retval = pthread_mutex_init(ps->mutex + mutexno, NULL);
    if (retval == -1) {
      perror("pthread_mutex_init");
      return -1;
    }

  }
      
  retval = pthread_cond_init(&(ps->worker_calccond), NULL);
  if (retval == -1) {
    perror("pthread_cond_init");
    return -1;
  }

  ps->worker_calc = malloc(sizeof(pthread_mutex_t) * num_threads);
  if (ps->worker_calc == NULL) {
    perror("malloc");
    return -1;
  }

  for (threadno = 0; threadno < num_threads; threadno++) {  

    retval = pthread_mutex_init(ps->worker_calc + threadno, NULL);
    if (retval == -1) {
      perror("pthread_mutex_init");
      return -1;
    }

  }
  
  retval = pthread_cond_init(&(ps->mainproc_cond), NULL);
  if (retval == -1) {
    perror("pthread_cond_init");
    return -1;
  }

  retval = pthread_mutex_init(&(ps->mainproc_mutex), NULL);
  if (retval == -1) {
    perror("pthread_mutex_init");
    return -1;
  }

  retval = pthread_cond_init(&(ps->worker_advancecond), NULL);
  if (retval == -1) {
    perror("pthread_cond_init");
    return -1;
  }
  
  ps->worker_advance = malloc(sizeof(pthread_mutex_t) * num_threads);
  if (ps->worker_advance == NULL) {
    perror("malloc");
    return -1;
  }

  for (threadno = 0; threadno < num_threads; threadno++) {  

    retval = pthread_mutex_init(ps->worker_advance + threadno, NULL);
    if (retval == -1) {
      perror("pthread_mutex_init");
      return -1;
    }

  }

  ps->workpacks_pergeneration = (pso_workpacks * num_threads);
  
  ps->prw.generationno = 0;
  
  ps->prw.workid = 0;

  ps->prw.thread_states = malloc(sizeof(uint64_t) * num_threads);
  if (ps->prw.thread_states == NULL) {
    perror("malloc");
    return -1;
  }

  ps->state = PSO_PROCESSING;
  
  return 0;
  
}

int alloc_ppack(particlepack *pp, long int num_particles) {

  pp->vcur = malloc(num_particles * sizeof(vec3d));
  if (pp->vcur == NULL) {
    perror("malloc");
    return -1;
  }

  pp->xcur = malloc(num_particles * sizeof(point3d_t));  
  if (pp->xcur == NULL) {
    perror("malloc");
    return -1;
  }

  pp->fitness = malloc(num_particles * sizeof(double));  
  if (pp->fitness == NULL) {
    perror("malloc");
    return -1;
  }

  pp->pbest = malloc(num_particles * sizeof(point3d_t));  
  if (pp->pbest == NULL) {
    perror("malloc");
    return -1;
  }

  return 0;

}

int alloc_threadprogress(generation_results *gr, long int num_generations) {

  gr->threadprogress_totalerrs = malloc(sizeof(double) * num_generations);
  if (gr->threadprogress_totalerrs == NULL) {
    perror("malloc");
    return -1;
  }
  
  return 0;
  
}

workpack_markstatus alloc_wpstatus(long int num_packs) {

  workpack_markstatus wp;

  wp.wpstatus = malloc(sizeof(long int) * num_packs);
  if (wp.wpstatus == NULL) {
    perror("malloc");
    exit(EXIT_FAILURE);
  }
  
  return wp;

}

int clear_wpstatus(workpack_markstatus *wp, long int num_packs) {
  
  memset(wp->wpstatus, 0, sizeof(long int) * num_packs);

  return 0;
  
}

int show_wpstatus(workpack_markstatus *wp, long int num_packs) {
  
  long int packno;

  fprintf(stderr, "%s: ", __FUNCTION__);
  
  for (packno = 0; packno < num_packs; packno++) {
    fprintf(stderr, "%ld ", wp->wpstatus[packno]);
  }

  fprintf(stderr, "\n");
  
  return 0;
  
}

point3d_t assign_pnt(point3d_t *a, point3d_t *b, point3d_t *c) {

  point3d_t pnta;

  pnta.x = (a->x + b->y + c->z) / 3.0;
  pnta.y = (a->x + b->y + c->z) / 3.0;
  pnta.z = (a->x + b->y + c->z) / 3.0;

  return pnta;
  
}

int output_asc(pso *ps, char *filename) {
  
  int asc_fd;
  mode_t mode;

  char strbuf[240];

  ssize_t bytes_written;

  size_t len;

  long int particleno;

  particlepack *pp;

  point3d_t *particles;

  int retval;

  point3d_t pnta;
  
  pp = &(ps->pp);
  
  fprintf(stderr, "%s: Output final vertex points %s\n", __FUNCTION__, filename);
    
  mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  asc_fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, mode);

  particles = pp->pbest;
  
  for (particleno = 0; particleno < ps->num_particles; particleno++) {

    pnta  = particles[particleno];
    
    retval = sprintf(strbuf, "%g %g %g\n", pnta.x, pnta.y, pnta.z);

    len = strlen(strbuf);
      
    bytes_written = write(asc_fd, strbuf, len);
    if (bytes_written != len) {
      perror("write");
      return -1;
    }
      
  }

  retval = close(asc_fd);
  if (retval == -1) {
    perror("close");
    return -1;
  }
  
  return 0;
  
}

int output_geomview(pso *ps, FILE *fp_out) {

  long int particleno;

  particlepack *pp;

  pixel_t color;
  
  fp_out = stdout;

  pp = &(ps->pp);
    
  printf("VECT\n");

  printf("%ld %ld %ld\n", ps->num_particles, ps->num_particles, ps->num_particles);

  for (particleno = 0; particleno < ps->num_particles; particleno++) {

    printf("%d ", 1);

  }

  putchar('\n');

  for (particleno = 0; particleno < ps->num_particles; particleno++) {

    printf("%d ", 1);

  }

  putchar('\n');

  for (particleno = 0; particleno < ps->num_particles; particleno++) {
    printf("%g %g %g\n", pp->pbest[particleno].x, pp->pbest[particleno].y, pp->pbest[particleno].z);
  }

  color.r = 32767.5 + 16383.75 * pp->xcur[particleno].y;
  color.g = 0;
  color.b = 0;
    
  for (particleno = 0; particleno < ps->num_particles; particleno++) {
    printf("%g %g %g %g\n", color.r / 65535.0, color.g / 65535.0, color.b / 65535.0, 1.0); 
  }

  return 0;

}

int fill_initialrnd(particlepack *pp, long int num_particles, int rnd_fd, double blo, double bup) {

  uint64_t rnds[6];

  int retval;

  long int particleno;

  double span;

  span = (bup - blo);
  
  for (particleno = 0; particleno < num_particles; particleno++) {

    retval = read(rnd_fd, &rnds, 6 * sizeof(uint64_t));
    if (retval != 6 * sizeof(uint64_t)) {
      perror("read");
      return -1;
    }
    
    pp->vcur[particleno][0] = -span + (2.0 * span * rnds[0]) / 18446744073709551615.0;
    pp->vcur[particleno][1] = -span + (2.0 * span * rnds[1]) / 18446744073709551615.0;
    pp->vcur[particleno][2] = -span + (2.0 * span * rnds[2]) / 18446744073709551615.0;    
    
    pp->xcur[particleno].x = blo + (span * rnds[3]) / 18446744073709551615.0;
    pp->xcur[particleno].y = blo + (span * rnds[4]) / 18446744073709551615.0;
    pp->xcur[particleno].z = blo + (span * rnds[5]) / 18446744073709551615.0;    
    
  }
  
  return 0;

}

#define RUNNING 0x1
#define HAVEWORK 0x2
#define EMPTYQUEUE 0x4
#define WAITCOND_NEXTGEN 0x8

typedef struct {

  uint64_t *state;

  pso *ps;

  particlepack *pp;

  long int threadno;

  int rnd_fd;

  double (*fitness_func)(point3d_t *pnta, void *extra);

  void *ff_extra;
  
  struct tailhead *headp;

  workunit *recalc;
  
  extraparams *ep;

  uint64_t *thread_states;

  long int *generationno;

} threadpack;

double elap_cur(struct timespec *start, extraparams *ep) {

  double elapsed;

  clock_gettime(CLOCK_REALTIME, &(ep->now));

  elapsed = (ep->now.tv_sec - start->tv_sec);

  return elapsed;
  
}

int psogen_csv(generation_results *gr, long int num_generations, char *outcsv_fn)  {

  FILE *csv_fp;

  char csv_delim;
  
  long int generationno;

  int retval;
  
  csv_fp = fopen(outcsv_fn, "w");
  if (csv_fp == NULL) {
    perror("fopen");
    return -1;
  }

  csv_delim = ',';

  fprintf(csv_fp, "Generation #%c Total Error\n", csv_delim);
    
  for (generationno = 0; generationno < num_generations; generationno++) {
      
    fprintf(csv_fp, "%ld%c %g\n", generationno, csv_delim, gr->threadprogress_totalerrs[generationno]);
      
  }
  
  retval = fclose(csv_fp);
  if (retval == -1) {
    perror("fclose");
    return -1;
  }

  return 0;

}

int postwork_progress(pso *ps, double totalerr, long int generationno, struct timespec *then, struct timespec *now, long int *progressdisp_freq) {

  {
    
    double percent;

    double remaining_sec;

    double estimated_totalsec;

    memcpy(then, now, sizeof(struct timespec));

    clock_gettime(CLOCK_REALTIME, now);

    if (now->tv_sec - then->tv_sec > 5) {

      progressdisp_freq[0] = 5;
	  
    }
	
    percent = generationno; percent /= ps->num_generations;

    {

      double elapsed_sec;

      elapsed_sec = (now->tv_sec - ps->start.tv_sec);

      estimated_totalsec = (ps->num_generations * elapsed_sec) / generationno;
	  
      remaining_sec = (estimated_totalsec - elapsed_sec);

    }
       
    ps->progress_func(percent, remaining_sec, generationno, totalerr);

  }
            
  return 0;

}

int advance_generation(pso *ps, particlepack *pp, long int workpacks_pergeneration) {

  int retval;

  pthread_mutex_lock(ps->mutex + MGENERATION);
  ps->prw.generationno++;
  pthread_mutex_unlock(ps->mutex + MGENERATION);

  pthread_mutex_lock(ps->mutex + MGENERATION);  
  retval = generate_work(ps, workpacks_pergeneration, ps->prw.generationno);
  if (retval == -1) {
    printf("%s: Trouble with call to generate_work.\n", __FUNCTION__);
    return -1;
  }
  pthread_mutex_unlock(ps->mutex + MGENERATION);
  
  return 0;

}

int shutdown_threads(uint64_t *thread_states, long int num_threads) {

  long int threadno;

  for (threadno = 0; threadno < num_threads; threadno++) {
    thread_states[threadno] = 0;
  }

  return 0;

}
  
int process_regwork(pso *ps, threadpack *tp, workunit *work) {

  double percent;

  particlepack *pp;

  pp = &(ps->pp);

  dfprintf(stderr, "%s: Processing work item %ld for generationno %ld\n", __FUNCTION__, work->packno, work->generationno);
    
  {
    
    percent = work->generationno; percent /= ps->num_generations;

    update_rangepso(pp, ps->mutex + MBESTPARTICLE, percent, tp->rnd_fd, work->particle_start, work->particle_end, tp->fitness_func, tp->ff_extra, ps->blo, ps->bup);

    pthread_mutex_lock(ps->mutex + MPACK);
    ps->wp.wpstatus[ work->packno ] |= 1;
    pthread_mutex_unlock(ps->mutex + MPACK);

    if (debug_level > 1) {
      dfprintf(stderr, "%s: Working on range %ld through %ld\n", __FUNCTION__, work->particle_start, work->particle_end);
    }
      
  }

  return 0;

}

int process_calcwork(pso *ps, threadpack *tp, workunit *work) {

  double totalerr;
  
  particlepack *pp;

  extraparams *ep;
  
  int retval;
  
  pp = &(ps->pp);

  ep = tp->ep;

  dfprintf(stderr, "%s: Processing work item %ld for generationno %ld\n", __FUNCTION__, work->packno, work->generationno);
  
  {

    pthread_mutex_lock(ps->mutex + MPARTICLE);
    totalerr = calc_totalerror(pp->pbest, ps->num_particles, tp->fitness_func, tp->ff_extra);
    pthread_mutex_unlock(ps->mutex + MPARTICLE);    

    if (isnan(totalerr)) {

      fprintf(stderr, "%s: Aborting because we have an error result of nan in one of the particles.\n", __FUNCTION__);

      exit(EXIT_FAILURE);

    }

    if (totalerr <= 0.0) {

      fprintf(stderr, "%s: Total error is quite low, early completion of PSO algorithm.\n", __FUNCTION__);

      return 1;
	
    }
    
    pthread_mutex_lock(ps->mutex + MEXTRAPAR);
    ep->totalerr_prev = ep->totalerr;
    ep->totalerr = totalerr;
    ep->totalerr_avg += totalerr;
    ep->totalerr_avg *= 0.5;
    pthread_mutex_unlock(ps->mutex + MEXTRAPAR);        
    
    pthread_mutex_lock(ps->mutex + MERRCALC);
    ps->gr.threadprogress_totalerrs[work->generationno] = totalerr;
    pthread_mutex_unlock(ps->mutex + MERRCALC);

    dfprintf(stderr, "%s: Set totalerr value into threadprogress_totalerrs.\n", __FUNCTION__);
    
    if (ep->totalerr > ep->totalerr_avg) {

      fprintf(stderr, "%s: Total error is heading in the wrong direction.\n", __FUNCTION__);
      fprintf(stderr, "%s: totalerr %g totalerr_prev %g totalerr_avg %g\n", __FUNCTION__, ep->totalerr, ep->totalerr_avg, ep->totalerr_prev);

      return -1;
	
    }

    pthread_mutex_lock(ps->mutex + MERRCALC);
    if ( (!(work->generationno % ep->progressdisp_freq)) && (ps->progress_func != NULL)) {
      retval = postwork_progress(ps, totalerr, work->generationno, &(ep->then), &(ep->now), &(ep->progressdisp_freq));
    }
    pthread_mutex_unlock(ps->mutex + MERRCALC);
    
    if (retval == -1) {
      dfprintf(stderr, "%s[%ld]: Leaving due to failure with postwork_progress.\n", __FUNCTION__, tp->threadno);
      return -1;
    }

    if (retval == 1) {
      printf("%s[%ld]: Early completion.\n", __FUNCTION__, tp->threadno);
      return 0;
      
    }
    
  }
  
  return 0;

}
  
void *pso_routine(void *extra) {

  void *ret;

  pso *ps;

  threadpack *tp;

  struct entry *n2;

  workunit work;
  
  int retval;

  uint64_t counter;

  workunit *works;

  long int num_works;

  workunit *rwork_ptr;

  workunit *ework_ptr;

  long int previous_generationno;
  
  uint64_t calcwait_counter;

  long int workno;
  
  tp = (threadpack*) extra;

  ps = tp->ps;

  ret = NULL;

  counter = 0;

  num_works = 10;

  n2 = NULL;
  
  works = malloc(sizeof(workunit) * num_works);
  if (works == NULL) {
    perror("malloc");
    return ret;
  }

  rwork_ptr = works + 0;
  ework_ptr = works + 1;

  previous_generationno = -1;

  calcwait_counter = 0;
  
  memset(&work, 0, sizeof(workunit));

  while (tp->state[0] & RUNNING) {

    dfprintf(stderr, "%s[%ld](counter %lu) Thread work.\n", __FUNCTION__, tp->threadno, counter);

    for (workno = 0; workno < num_works; workno++) {
    
      if (tp->state[0] & EMPTYQUEUE) {
	break;
      }

      dfprintf(stderr, "%s[%ld]: !EMPTYQUEUE\n", __FUNCTION__, tp->threadno);

      workno = 0;
      
      pthread_mutex_lock(ps->mutex + MQUEUE);
      if (!TAILQ_EMPTY(ps->headp)) {
	n2 = TAILQ_FIRST(ps->headp);
	work = n2->work;
	TAILQ_REMOVE(ps->headp, n2, entries);
	pthread_mutex_lock(ps->mutex + MTHREADSTATE);
	tp->state[0] |= HAVEWORK;
	pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      } else {
	pthread_mutex_lock(ps->mutex + MTHREADSTATE);
	tp->state[0] |= EMPTYQUEUE;
	pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      }
      pthread_mutex_unlock(ps->mutex + MQUEUE);

      if ( (tp->state[0] & HAVEWORK) && (tp->state[0] & RUNNING) ) {

	dfprintf(stderr, "%s[%ld](%gs): HAVEWORK(%ld) %ld %ld for generation %ld\n", __FUNCTION__, tp->threadno, elap_cur(&(ps->start), tp->ep), work.packno, work.particle_start, work.particle_end, work.generationno);
	
	free(n2);      

	retval = process_regwork(ps, tp, &work);
	if (retval == -1) {

	  fprintf(stderr, "%s[%ld]: Failure with call to process_regwork.\n", __FUNCTION__, tp->threadno);

	  pthread_mutex_lock(ps->mutex + MTHREADSTATE);
	  shutdown_threads(ps->prw.thread_states, ps->num_threads);
	  pthread_mutex_unlock(ps->mutex + MTHREADSTATE);

	  pthread_mutex_lock(ps->mutex + MTHREADSTATE);
	  show_threads(ps->prw.thread_states, ps->num_threads);
	  pthread_mutex_unlock(ps->mutex + MTHREADSTATE);      
	  
	  break;
	  
	}

	pthread_mutex_lock(ps->mutex + MTHREADSTATE);
	tp->state[0] &= (~HAVEWORK);
	pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
	
      }

    }

    if ( (tp->state[0] & EMPTYQUEUE) && (tp->state[0] & RUNNING) ) {

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);
      tp->state[0] |= WAITCOND_NEXTGEN;
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      
    }

    if (tp->state[0] & WAITCOND_NEXTGEN && (tp->state[0] & RUNNING)) {

      pthread_mutex_lock(ps->mutex + MGENERATION);
      previous_generationno = ps->prw.generationno;
      dfprintf(stderr, "%s[%ld]: WAITCOND_NEXTGEN for generationno %ld previous_generationno %ld\n", __FUNCTION__, tp->threadno, ps->prw.generationno, previous_generationno);
      pthread_mutex_unlock(ps->mutex + MGENERATION);

      if (previous_generationno >= (ps->num_generations - 1)) {

	dfprintf(stderr, "%s[%ld](%gs): Thread has reached last generation, no longer advancing.\n", __FUNCTION__, tp->threadno, elap_cur(&(ps->start), tp->ep));	
	
	break;
	
      }
      
      pthread_mutex_lock(ps->worker_advance + tp->threadno);

      while ( (tp->state[0] & WAITCOND_NEXTGEN) && (tp->state[0] & RUNNING)) {

	retval = pthread_cond_wait(&(ps->worker_advancecond), ps->worker_advance + tp->threadno);
	if (retval == -1) {
	  perror("pthread_cond_wait");
	  free(works);
	  return ret;
	}

	pthread_mutex_lock(ps->mutex + MTHREADSTATE);
	tp->state[0] &= ~(WAITCOND_NEXTGEN);
	pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
	
      }

      pthread_mutex_unlock(ps->worker_advance + tp->threadno);
      
      dfprintf(stderr, "%s[%ld](%gs): Received broadcast wakeup for next generation.\n", __FUNCTION__, tp->threadno, elap_cur(&(ps->start), tp->ep));

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);
      tp->state[0] &= ~(EMPTYQUEUE);
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);      
      
    }
    
    counter++;
    
  }
  
  free(works);
  
  dfprintf(stderr, "%s[%ld]: Leaving.\n", __FUNCTION__, tp->threadno);
  
  return ret;

}

int generate_work(pso *ps, long int num_workpacks, long int generationno) {

  long int varstart;

  long int particles_perpack;

  long int packno;
  
  struct entry *n1;  

  particles_perpack = (ps->num_particles / num_workpacks);
  
  for (packno = 0; packno < num_workpacks; packno++) {
      
    n1 = malloc(sizeof(struct entry));
    if (n1 == NULL) {
      perror("malloc");
      return -1;
    }

    varstart = (packno * particles_perpack);
    
    n1->work = (workunit) { .workid = ps->prw.workid, .generationno = generationno, .packno = packno, .particle_start = varstart, .particle_end = (varstart + particles_perpack), .wtype = WPROCESS_RANGE };
      
    pthread_mutex_lock(ps->mutex + MQUEUE);
    TAILQ_INSERT_TAIL(ps->headp, n1, entries);
    pthread_mutex_unlock(ps->mutex + MQUEUE);

    /*
    dfprintf(stderr, "%s: Generated workpack %ld with range %ld %ld for generation %ld\n", __FUNCTION__, n1->work.packno, n1->work.particle_start, n1->work.particle_end, n1->work.generationno);
    */
    
    ps->prw.workid++;

  }

  return 0;

}

workunit gen_workcalc(psorw *prw) {

  workunit work;

  work = (workunit) { .workid = -1, .generationno = prw->generationno, .packno = -1, .particle_start = 0, .particle_end = 0, .wtype = WCALCULATE_TOTALERR };

  return work;
  
}

int show_pso(pso *ps) {

  printf("Particles: %ld\n", ps->num_particles);  
  printf("Generations: %ld\n", ps->num_generations);
  printf("Threads: %ld\n", ps->num_threads);  

  return 0;

}

long int check_threadactivity(pso *ps, threadpack *tps) {

  long int allidle;

  long int threadno;
  
  for (threadno = 0; threadno < ps->num_threads; threadno++) {

    if (!(tps[threadno].state[0] & WAITCOND_NEXTGEN)) {

      break;

    }

  }

  allidle = (1 - (threadno < ps->num_threads));
    
  return allidle;
  
}

int show_work(workunit *work) {

  fprintf(stderr, "%s: workid %ld generationno %ld packno %ld particle_start %ld particle_end %ld wtype %ld\n", __FUNCTION__, work->workid, work->generationno, work->packno, work->particle_start, work->particle_end, work->wtype);

  return 0;
  
}
  
int show_threads(uint64_t *thread_states, long int num_threads) {

  long int threadno;

  dfprintf(stderr, "Thread summary: ");

  for (threadno = 0; threadno < num_threads; threadno++) {
    dfprintf(stderr, "%lu ", thread_states[threadno]);
  }

  dfprintf(stderr, "\n");
  
  return 0;
  
}

int set_initial(particlepack *pp, long int num_particles) {

  long int particleno;

  for (particleno = 0; particleno < num_particles; particleno++) {

    pp->fitness[particleno] = pp->fitness_func(pp->xcur + particleno, pp->ff_extra);
    
    pp->pbest[particleno] = pp->xcur[particleno];
      
  }
  
  pp->gbest = pp->xcur[0];
  pp->gbesterr = pp->fitness_func(&(pp->gbest), pp->ff_extra);
  
  for (particleno = 0; particleno < num_particles; particleno++) {

    if (pp->fitness[particleno] < pp->gbesterr) {

      pp->gbest = pp->xcur[particleno];
      pp->gbesterr = pp->fitness[particleno];
      
    }
      
  }

  return 0;

}

int waitnextgen_total(uint64_t *thread_states, long int num_threads) {

  long int threadno;

  long int total;

  total = 0;
  
  for (threadno = 0; threadno < num_threads; threadno++) {
    if (thread_states[threadno] & WAITCOND_NEXTGEN) {
      total++;
    }
  }

  return total;

}

enum { CALCCHECK = 1, PACKFILLED, ENSUREWAIT, GENWORK, CALCULATE, ADVANCE };

enum { MP_NONE, MP_CONTINUE, MP_SHUTDOWN };

int mp_monitor(pso *ps, particlepack *pp, threadpack *tp, extraparams *ep, long int *mainproc) {

  long int threadno;

  int retval;

  static uint64_t packscan_counter = 0;
  
  if (mainproc[0] == CALCCHECK) {

    dfprintf(stderr, "%s(%gs): CALCCHECK\n", __FUNCTION__, elap_cur(&(ps->start), ep));      

    pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
    for (threadno = 0; threadno < (ps->num_threads - 1); threadno++) {
      if (!(ps->prw.thread_states[threadno] & WAITCOND_NEXTGEN)) {
	break;
      }
    }
    pthread_mutex_unlock(ps->mutex + MTHREADSTATE);      
      
    if (threadno < (ps->num_threads - 1)) {

      dfprintf(stderr, "%s(%gs): CALCCHECK not all threads are in WAITCALC_COND.\n", __FUNCTION__, elap_cur(&(ps->start), ep));      

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
      show_threads(ps->prw.thread_states, ps->num_threads);
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      
      return MP_CONTINUE;
	
    }
      
    mainproc[0] = PACKFILLED;
      
  }
    
  if (mainproc[0] == PACKFILLED) {

    dfprintf(stderr, "%s(%gs): PACKFILLED\n", __FUNCTION__, elap_cur(&(ps->start), ep));

    pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
    show_threads(ps->prw.thread_states, ps->num_threads);
    pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      
    {

      long int packno;

      pthread_mutex_lock(ps->mutex + MPACK);
      for (packno = 0; packno < ps->workpacks_pergeneration; packno++) {
	if (!(ps->wp.wpstatus[packno])) {
	  packscan_counter++;
	  break;
	}
      }
      pthread_mutex_unlock(ps->mutex + MPACK);

      if (packno < ps->workpacks_pergeneration) {

	dfprintf(stderr, "%s(%gs): Still don't have pack filled for generationno %ld\n", __FUNCTION__, elap_cur(&(ps->start), ep), ps->prw.generationno);      

	for (threadno = 0; threadno < (ps->num_threads - 1); threadno++) {

	  if (ps->prw.thread_states[threadno] == (RUNNING | HAVEWORK)) {
	    break;
	  }
	  
	  if (ps->prw.thread_states[threadno] & WAITCOND_NEXTGEN) {
	    break;
	  }

	}

	dfprintf(stderr, "%s: threadno %ld (ps->num_threads - 1) %ld\n", __FUNCTION__, threadno, ps->num_threads - 1);
	  
	if (threadno == (ps->num_threads - 1)) {

	  dfprintf(stderr, "%s: No thread has gone into WAITCOND_NEXTGEN state yet!\n", __FUNCTION__);

	  pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
	  show_threads(ps->prw.thread_states, ps->num_threads);
	  pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      
	  return MP_CONTINUE;
	    
	}
	  
	dfprintf(stderr, "%s: At least one thread has gone into WAITCOND_NEXTGEN state.\n", __FUNCTION__);

	pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
	show_threads(ps->prw.thread_states, ps->num_threads);
	pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      
      }
	  
    }
      
    dfprintf(stderr, "%s: packscan_counter %lu\n", __FUNCTION__, packscan_counter); 

    mainproc[0] = ENSUREWAIT;
      
  }

  if (mainproc[0] == ENSUREWAIT) {

    long int incomplete;
    
    dfprintf(stderr, "%s(%gs): ENSUREWAIT\n", __FUNCTION__, elap_cur(&(ps->start), ep));

    incomplete = 0;
    
    for (threadno = 0; threadno < (ps->num_threads - 1); threadno++) {

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
      if (!(ps->prw.thread_states[threadno] & WAITCOND_NEXTGEN)) {
	incomplete |= 1;
      }
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);      

    }

    if (incomplete) {

      dfprintf(stderr, "%s(%gs): ENSUREWAIT Still waiting for all threads to become WAITCOND_NEXTGEN.\n", __FUNCTION__, elap_cur(&(ps->start), ep));

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
      show_threads(ps->prw.thread_states, ps->num_threads);
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      
      return MP_CONTINUE;
	
    }
      
    mainproc[0] = GENWORK;
      
  }
    
  if (mainproc[0] == GENWORK) {

    dfprintf(stderr, "%s(%gs): GENWORK\n", __FUNCTION__, elap_cur(&(ps->start), ep));
      
    pthread_mutex_lock(ps->mutex + MGENERATION);
    ps->recalc = gen_workcalc(&(ps->prw));
    pthread_mutex_unlock(ps->mutex + MGENERATION);
      
    mainproc[0] = CALCULATE;
	
  }

  if (mainproc[0] == CALCULATE) {

    dfprintf(stderr, "%s(%gs): CALCULATE\n", __FUNCTION__, elap_cur(&(ps->start), ep));

    pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
    show_threads(ps->prw.thread_states, ps->num_threads);
    pthread_mutex_unlock(ps->mutex + MTHREADSTATE);      
      
    retval = process_calcwork(ps, tp, &(ps->recalc));
    if (retval == -1) {

      fprintf(stderr, "%s[%ld]: Failure with call to process_calcwork.\n", __FUNCTION__, tp->threadno);

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);
      shutdown_threads(tp->thread_states, ps->num_threads);
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
	
      pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
      show_threads(tp->thread_states, ps->num_threads);
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);      
	
      return MP_SHUTDOWN;

    }
      
    mainproc[0] = ADVANCE;

  }
    
  if (mainproc[0] == ADVANCE) {

    long int wct;
      
    dfprintf(stderr, "%s(%gs): ADVANCE\n", __FUNCTION__, elap_cur(&(ps->start), ep));

    dfprintf(stderr, "%s: Advancing to next generation and signalling broadcast.\n", __FUNCTION__);

    for (threadno = 0; threadno < (ps->num_threads - 1); threadno++) {

      if (!(ps->prw.thread_states[threadno] & EMPTYQUEUE)) {

	break;

      }	  

    }

    if (threadno < (ps->num_threads - 1)) {

      dfprintf(stderr, "%s: Not all of the threads have made it to the EMPTYQUEUE waiting position.\n", __FUNCTION__);

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
      show_threads(ps->prw.thread_states, ps->num_threads);
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
      
      return MP_CONTINUE;
	
    }
      
    for (threadno = 0; threadno < (ps->num_threads - 1); threadno++) {

      if (!(ps->prw.thread_states[threadno] & WAITCOND_NEXTGEN)) {

	break;

      }	  

    }

    wct = waitnextgen_total(ps->prw.thread_states, ps->num_threads);

    dfprintf(stderr, "%s: ADVANCE waitcalc_total %ld\n", __FUNCTION__, wct);
      
    if (threadno < (ps->num_threads - 1) && wct == 0) {

      fprintf(stderr, "%s: Not all of the threads have made it to nextgen waiting condition.\n", __FUNCTION__);

      pthread_mutex_lock(ps->mutex + MTHREADSTATE);      
      show_threads(ps->prw.thread_states, ps->num_threads);
      pthread_mutex_unlock(ps->mutex + MTHREADSTATE);      	

      return MP_CONTINUE;

    }
      
    {
	
      retval = advance_generation(ps, pp, ps->workpacks_pergeneration);
      if (retval == -1) {
	printf("%s[%ld]: Trouble with call to advance_generation.\n", __FUNCTION__, tp->threadno);
	return -1;
      }

      pthread_mutex_lock(ps->mutex + MPACK);
      clear_wpstatus(&(ps->wp), ps->workpacks_pergeneration);
      pthread_mutex_unlock(ps->mutex + MPACK);
	
    }

    dfprintf(stderr, "%s(%gs): Broadcast signalling advancecond.\n", __FUNCTION__, elap_cur(&(ps->start), ep));
	       
    retval = pthread_cond_broadcast(&(ps->worker_advancecond));
    if (retval == -1) {
      perror("pthread_cond_broadcast");
      return -1;
    }

    mainproc[0] = PACKFILLED;
      
  }

  return MP_NONE;

}
  
int process_pso(pso *ps, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra) {

  long int threadno;

  pthread_t *threads;

  threadpack *tps;
    
  int retval;

  int rnd_fd;

  particlepack *pp;

  extraparams ep;

  uint64_t wcounter;

  uint64_t mainstate;

  long int mainproc;

  threadpack *tp;

  uint64_t counter;

  workunit work;

  long int num_works;

  long int workno;

  struct entry *n2;
  
  ep.totalerr_prev = INFINITY;
  ep.totalerr = INFINITY;
  ep.totalerr_avg = INFINITY;
  ep.progressdisp_freq = 20;
  
  clock_gettime(CLOCK_REALTIME, &(ps->start));
  
  pp = &(ps->pp);
  
  threads = malloc(sizeof(pthread_t) * ps->num_threads);
  if (threads == NULL) {
    perror("malloc");
    return -1;
  }

  tps = malloc(sizeof(threadpack) * ps->num_threads);
  if (tps == NULL) {
    perror("malloc");
    free(threads);
    return -1;
  }

  {
    long int totalerrno;
    for (totalerrno = 0; totalerrno < ps->num_generations; totalerrno++) {
      ps->gr.threadprogress_totalerrs[totalerrno] = NAN;
    }
  }

  rnd_fd = open("/dev/urandom", O_RDONLY);
  if (rnd_fd == -1) {
    perror("open");
    free(threads);
    free(tps);
    return -1;
  }

  fprintf(stderr, "%s: Filling with initial random values.\n", __FUNCTION__);
  
  retval = fill_initialrnd(pp, ps->num_particles, rnd_fd, ps->blo, ps->bup);
  if (retval == -1) {
    printf("%s: Trouble with call to fill_initialrnd.\n", __FUNCTION__);
    return -1;
  }

  retval = set_initial(pp, ps->num_particles);
    
  {
  
    fprintf(stderr, "%s: Generating work units for first generation.\n", __FUNCTION__);

    ps->prw.generationno = 0;

    retval = generate_work(ps, ps->workpacks_pergeneration, ps->prw.generationno);
    if (retval == -1) {
      printf("%s: Trouble with call to generate_work.\n", __FUNCTION__);
      return -1;
    }

  }
  
  for (threadno = 0; threadno < ps->num_threads; threadno++) {

    tps[threadno].state = ps->prw.thread_states + threadno;
    
    tps[threadno].state[0] = RUNNING;
    tps[threadno].ps = ps;
    tps[threadno].pp = &(ps->pp);
    tps[threadno].threadno = threadno;
    tps[threadno].rnd_fd = rnd_fd;
    tps[threadno].fitness_func = fitness_func;
    tps[threadno].ff_extra = ff_extra;
    tps[threadno].headp = ps->headp;
    tps[threadno].recalc = &(ps->recalc);
    
    tps[threadno].ep = &ep;

    tps[threadno].thread_states = ps->prw.thread_states;

    tps[threadno].generationno = &(ps->prw.generationno);

  }

  clear_wpstatus(&(ps->wp), pso_workpacks);
  
  for (threadno = 0; threadno < (ps->num_threads - 1); threadno++) {

    retval = pthread_create(threads + threadno, NULL, pso_routine, tps + threadno);
    
  }

  wcounter = 0;

  enum { MAINPROC_WORK = 0x1 };

  mainstate = MAINPROC_WORK;
  
  mainproc = PACKFILLED;

  tp = tps + (ps->num_threads - 1);
  
  if (mainstate & MAINPROC_WORK) {

    tp->state[0] = RUNNING;
    
  }

  counter = 0;

  memset(&work, 0, sizeof(workunit));

  num_works = 10;

  n2 = NULL;
  
  for (;;) {

    dfprintf(stderr, "%s[%ld]: generationno %ld tp->state %lu mainproc %ld\n", __FUNCTION__, counter, ps->prw.generationno, tp->state[0], mainproc);
    
    counter++;
    
    if (mainstate & MAINPROC_WORK) {
      if (!(tp->state[0] & RUNNING)) {
	fprintf(stderr, "%s: We have shutdown.\n", __FUNCTION__);
	break;
      }
    }
    
    while ( (tp->state[0] & RUNNING) && (mainstate & MAINPROC_WORK) ) {

      for (workno = 0; workno < num_works; workno++) {
      
	if (!(tp->state[0] & EMPTYQUEUE) && (tp->state[0] & RUNNING) ) {
      
	  pthread_mutex_lock(ps->mutex + MQUEUE);
	  if (!TAILQ_EMPTY(ps->headp)) {
	    n2 = TAILQ_FIRST(ps->headp);
	    work = n2->work;
	    TAILQ_REMOVE(ps->headp, n2, entries);
	    pthread_mutex_lock(ps->mutex + MTHREADSTATE);	
	    tp->state[0] |= HAVEWORK;
	    pthread_mutex_unlock(ps->mutex + MTHREADSTATE);	
	  } else {
	    pthread_mutex_lock(ps->mutex + MTHREADSTATE);	
	    tp->state[0] |= EMPTYQUEUE;
	    pthread_mutex_unlock(ps->mutex + MTHREADSTATE);	
	  }
	  pthread_mutex_unlock(ps->mutex + MQUEUE);

	}
      
	if ( (tp->state[0] & HAVEWORK) && (tp->state[0] & RUNNING) ) {

	  free(n2);
	
	  retval = process_regwork(ps, tp, &work);
	  if (retval == -1) {
	    printf("%s: Failure with call to process_regwork.\n", __FUNCTION__);
	    break;
	  }

	  pthread_mutex_lock(ps->mutex + MTHREADSTATE);	
	  tp->state[0] &= (~HAVEWORK);
	  pthread_mutex_unlock(ps->mutex + MTHREADSTATE);	
	
	  wcounter++;

	  break;
	  
	}

      }
	
      if ( (tp->state[0] & EMPTYQUEUE) && (tp->state[0] & RUNNING) ) {

	pthread_mutex_lock(ps->mutex + MTHREADSTATE);	
	tp->state[0] |= WAITCOND_NEXTGEN;
	pthread_mutex_unlock(ps->mutex + MTHREADSTATE);	
	
	break;
	
      }
      
      if ( (tp->state[0] & WAITCOND_NEXTGEN) && (tp->state[0] & RUNNING) ) {

	pthread_mutex_lock(ps->mutex + MTHREADSTATE);
	tp->state[0] &= ~(EMPTYQUEUE);
	pthread_mutex_unlock(ps->mutex + MTHREADSTATE);	
	
	break;
	
      }
            
    }

    pthread_mutex_lock(ps->mutex + MGENERATION);
    if (ps->prw.generationno >= (ps->num_generations - 1)) {
      dfprintf(stderr, "\n%s: Completed generations.\n", __FUNCTION__);
      ps->state = PSO_COMPLETE;
    }
    pthread_mutex_unlock(ps->mutex + MGENERATION);

    if (ps->state != PSO_PROCESSING) {
      dfprintf(stderr, "%s: Main processing loop complete.\n", __FUNCTION__);
      break;
    }

    retval = mp_monitor(ps, pp, tp, &ep, &mainproc);
    switch(retval) {
    case MP_NONE: break;
    case MP_CONTINUE:
      dfprintf(stderr, "%s: Received continue from mp_monitor.\n", __FUNCTION__);
      break;
    case -1:
      dfprintf(stderr, "%s: Trouble with call to mp_followup to monitor thread activity.\n", __FUNCTION__);
      return -1;
    case MP_SHUTDOWN:
      dfprintf(stderr, "%s: mp_monitor has shut down main processing loop.\n", __FUNCTION__);      
      
      ps->state = PSO_COMPLETE;
      break;
    }
    
    pthread_mutex_lock(ps->mutex + MTHREADSTATE);
    dfprintf(stderr, "%s: ps->prw.generationno %ld ps->num_generations %ld\n", __FUNCTION__, ps->prw.generationno, ps->num_generations);
    pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
    
  }

  {

    dfprintf(stderr, "%s: Final total error calculation.\n", __FUNCTION__);
    
    pthread_mutex_lock(ps->mutex + MGENERATION);
    ps->recalc = gen_workcalc(&(ps->prw));
    pthread_mutex_unlock(ps->mutex + MGENERATION);
  
    retval = process_calcwork(ps, tp, &(ps->recalc));
    if (retval == -1) {

      fprintf(stderr, "%s[%ld]: Failure with call to process_calcwork.\n", __FUNCTION__, tp->threadno);

      return -1;

    }

  }
    
  dfprintf(stderr, "%s: Closing down threads.\n", __FUNCTION__);
  
  pthread_mutex_lock(ps->mutex + MTHREADSTATE);
  shutdown_threads(ps->prw.thread_states, ps->num_threads);
  pthread_mutex_unlock(ps->mutex + MTHREADSTATE);
  
  dfprintf(stderr, "%s: Joining threads.\n", __FUNCTION__);

  for (threadno = 0; threadno < (ps->num_threads - 1); threadno++) {
    retval = pthread_join(threads[threadno], NULL);
    if (retval == -1) {
      perror("pthread_join");
      return -1;
    }
  }
      
  dfprintf(stderr, "\n%s: PSO algorithm complete.\n", __FUNCTION__);   

  {

    double elapsed;

    elapsed = (ep.now.tv_sec - ps->start.tv_sec);

    printf("Elapsed runtime %gs\n", elapsed);

  }
  
  free(threads);
  free(tps);
  
  retval = close(rnd_fd);
  if (retval == -1) {
    perror("close");
    return -1;
  }
  
  return 0;
  
}

int close_pso(pso *ps) {

  int retval;

  long int threadno;

  long int mutexno;
  
  free(ps->pp.vcur);
  free(ps->pp.xcur);

  free(ps->pp.fitness);
  free(ps->pp.pbest);

  free(ps->gr.threadprogress_totalerrs);

  for (mutexno = 0; mutexno < num_mutex; mutexno++) {
    retval = pthread_mutex_destroy(ps->mutex + mutexno);
  }

  free(ps->mutex);
  
  retval = pthread_cond_destroy(&(ps->worker_calccond));

  for (threadno = 0; threadno < ps->num_threads; threadno++) {

    retval = pthread_mutex_destroy(ps->worker_calc + threadno);  

  }
  
  retval = pthread_cond_destroy(&(ps->mainproc_cond));
  retval = pthread_mutex_destroy(&(ps->mainproc_mutex));  
  if (retval == -1) {
    perror("pthread_mutex_destroy");
    return -1;
  }

  retval = pthread_cond_destroy(&(ps->worker_advancecond));
  
  for (threadno = 0; threadno < ps->num_threads; threadno++) {

    retval = pthread_mutex_destroy(ps->worker_advance + threadno);  

  }

  free(ps->worker_calc);
  
  free(ps->worker_advance);
  
  free(ps->prw.thread_states);
  
  return 0;
    
}
  
int debug_frame(particlepack *pp, long int num_particles) {

  long int particleno;
  
  for (particleno = 0; particleno < num_particles; particleno += 1000) {

    fprintf(stderr, "%s: Fitness[%ld] %g\n", __FUNCTION__, particleno, pp->fitness[particleno]);
    fprintf(stderr, "%s: x[%ld] %g %g %g\n", __FUNCTION__, particleno, pp->xcur[particleno].x, pp->xcur[particleno].y, pp->xcur[particleno].z);      
      
  }

  return 0;
    
}
