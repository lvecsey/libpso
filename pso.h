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

#ifndef PSO_H
#define PSO_H

#include <sys/queue.h>

#include "particlepack.h"

#define pso_workpacks 100

enum { WPROCESS_RANGE = 1, WCALCULATE_TOTALERR };

typedef struct {

  long int workid;

  long int generationno;

  long int packno;

  long int particle_start, particle_end;

  long int wtype;
  
} workunit;

struct entry {

  workunit work;
  
  TAILQ_ENTRY(entry) entries;

};

typedef struct {

  double *threadprogress_totalerrs;
  
} generation_results;

typedef struct {

  long int *wpstatus;
  
} workpack_markstatus;

typedef struct {

  long int generationno;

  long int workid;
  
  uint64_t *thread_states;

} psorw;

enum { MQUEUE, MERRCALC, MPARTICLE, MPACK, MEXTRAPAR, MBESTPARTICLE, MGENERATION, MTHREADSTATE };

#define num_mutex 8

enum { PSO_INIT, PSO_PROCESSING, PSO_COMPLETE };

typedef struct {

  uint64_t state;
  
  long int num_particles;

  long int num_generations;
  
  long int num_threads;

  particlepack pp;

  double blo, bup;
  
  struct tailhead *headp;

  workunit recalc;

  pthread_mutex_t *mutex;
    
  pthread_cond_t worker_calccond;
  pthread_mutex_t *worker_calc;

  pthread_cond_t mainproc_cond;
  pthread_mutex_t mainproc_mutex;

  pthread_cond_t worker_advancecond;
  pthread_mutex_t *worker_advance;

  generation_results gr;

  workpack_markstatus wp;

  long int workpacks_pergeneration;
  
  int (*progress_func)(double percent, double eta, long int generationno, double combined_err);

  struct timespec start;

  psorw prw;
  
} pso;

#include "mini_gxkit.h"

extern const long int def_numpart;

extern const long int def_numgenerations;

int init_pso(pso *ps, long int num_threads);

int alloc_ppack(particlepack *pp, long int num_particles);

int alloc_threadprogress(generation_results *gr, long int num_generations);

workpack_markstatus alloc_wpstatus(long int num_packs);

int process_pso(pso *ps, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra);

int show_pso(pso *ps);

int output_asc(pso *ps, char *filename);

int output_geomview(pso *ps, FILE *fp_out);

int psogen_csv(generation_results *gr, long int num_generations, char *outcsv_fn);

int close_pso(pso *ps);

double calc_totalerror(point3d_t *particles, long int num_particles, double (*fitness_func)(point3d_t *pnta, void *extra), void *ff_extra);

double fetch_rnd(int rnd_fd);
int fill_initialrnd(particlepack *pp, long int num_particles, int rnd_fd, double blo, double bup);
int debug_frame(particlepack *pp, long int num_particles);

#endif
