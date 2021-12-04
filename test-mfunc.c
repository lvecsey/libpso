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

#include <string.h>

#include <math.h>

#include "mini_gxkit.h"

#include "mfuncs.h"

#include "pso.h"

#include "pso_util.h"

#include "fitnesspack.h"

#include "fitnessfunc_extra.h"

double pyramid_fitnessfunc(point3d_t *pnta, void *extra) {

  point3d_t expected_pntc;

  double err;

  fitnessfunc_extra *ff_extra;

  ff_extra = (fitnessfunc_extra*) extra;

  expected_pntc = math3d_pyramid(pnta->x, pnta->y);

  expected_pntc.y *= ff_extra->sf;

  if (1 == within_range(pnta, ff_extra->blo, ff_extra->bup)) {
  
    err = calc_single(pnta, &expected_pntc);

  }

  else {

    err = 1e6;

  }
  
  return err;
  
}

double bumps_fitnessfunc(point3d_t *pnta, void *extra) {

  point3d_t expected_pntc;

  double err;

  fitnessfunc_extra *ff_extra;

  ff_extra = (fitnessfunc_extra*) extra;

  expected_pntc = math3d_bumps(pnta->x, pnta->y);

  expected_pntc.y *= ff_extra->sf;

  if (1 == within_range(pnta, ff_extra->blo, ff_extra->bup)) {
  
    err = calc_single(pnta, &expected_pntc);

  }

  else {

    err = 1e6;

  }
    
  return err;
  
}

double sinc_fitnessfunc(point3d_t *pnta, void *extra) {

  point3d_t expected_pntc;

  double err;

  fitnessfunc_extra *ff_extra;

  ff_extra = (fitnessfunc_extra*) extra;

  expected_pntc = math3d_sinc(pnta->x, pnta->y);

  expected_pntc.y *= ff_extra->sf;
  
  if (1 == within_range(pnta, ff_extra->blo, ff_extra->bup)) {
  
    err = calc_single(pnta, &expected_pntc);

  }

  else {

    err = 1e6;

  }
  
  return err;
  
}

double torus_fitnessfunc(point3d_t *pnta, void *extra) {

  point3d_t expected_pntc;

  double err;

  fitnessfunc_extra *ff_extra;

  ff_extra = (fitnessfunc_extra*) extra;

  expected_pntc = math3d_torus(pnta->x, pnta->y);

  expected_pntc.z *= ff_extra->sf;
  
  if (isnan(expected_pntc.z)) {
    expected_pntc.z = 0.0;
  }
  
  if (1 == within_range(pnta, ff_extra->blo, ff_extra->bup)) {
  
    err = calc_single(pnta, &expected_pntc);

  }

  else {

    err = 1e6;

  }
  
  return err;
  
}

double tube_fitnessfunc(point3d_t *pnta, void *extra) {

  point3d_t expected_pntc;

  double err;

  fitnessfunc_extra *ff_extra;

  ff_extra = (fitnessfunc_extra*) extra;

  expected_pntc = math3d_tube(pnta->x, pnta->y);

  expected_pntc.z *= ff_extra->sf;
  
  if (isnan(expected_pntc.z)) {
    expected_pntc.z = 0.0;
  }

  if (1 == within_range(pnta, ff_extra->blo, ff_extra->bup)) {
  
    err = calc_single(pnta, &expected_pntc);

  }

  else {

    err = 1e6;

  }
  
  return err;
  
}

double stairs_fitnessfunc(point3d_t *pnta, void *extra) {

  point3d_t expected_pntc;

  double err;

  fitnessfunc_extra *ff_extra;

  ff_extra = (fitnessfunc_extra*) extra;

  expected_pntc = math3d_stairs(pnta->x, pnta->y);

  expected_pntc.z *= ff_extra->sf;
  
  if (isnan(expected_pntc.z)) {
    expected_pntc.z = 0.0;
  }

  if (1 == within_range(pnta, ff_extra->blo, ff_extra->bup)) {
  
    err = calc_single(pnta, &expected_pntc);

  }

  else {

    err = 1e6;

  }
  
  return err;
  
}

int progress(double percent, double eta, long int generationno, double combined_err) {

  fprintf(stderr, "\r%s(%0.04g%%) eta %.8gs Generation %ld total error %g       ", __FUNCTION__, 100.0 * percent, eta, generationno, combined_err);
  
  return 0;

}

#define def_numparticles (2000 * pso_workpacks)
#define def_numgenerations 250
#define def_numthreads 5

#define def_ascfn "mfunc.asc"

int main(int argc, char *argv[]) {

  pso ps;

  int retval;

  char *mfunc_str;
  
  char *outasc_fn;
  
  char *outcsv_fn;

  char *env_DEBUG;

  memset(&ps, 0, sizeof(pso));
  
  ps.state = PSO_INIT;
  
  ps.num_particles = argc>1 ? strtol(argv[1],NULL,10) : def_numparticles;
  
  ps.num_generations = argc>2 ? strtol(argv[2],NULL,10) : def_numgenerations;

  ps.num_threads = argc>3 ? strtol(argv[3],NULL,10) : def_numthreads;

  outasc_fn = argc>4 ? argv[4] : def_ascfn;

  mfunc_str = argc>5 ? argv[5] : "pyramid";
  
  outcsv_fn = "pso_generations.csv";
  
  ps.blo = -3.0;

  ps.bup = 3.0;

  ps.progress_func = progress;

  {
  
    ps.num_particles /= pso_workpacks;
    ps.num_particles *= pso_workpacks;

    show_pso(&ps);
    
  }

  {

    env_DEBUG = getenv("DEBUG");

    if (env_DEBUG != NULL) {

      ps.debug_level = strtol(env_DEBUG, NULL, 10);

      if (ps.debug_level < 0 || ps.debug_level > 10) {
	ps.debug_level = 0;
      }

      fprintf(stderr, "%s: Using debug level %ld\n", __FUNCTION__, ps.debug_level);
      
    }
    
  }
  
  fprintf(stderr, "%s: Initializing pso including worker thread items.\n", __FUNCTION__);
  
  retval = init_pso(&ps, ps.num_threads);
  if (retval == -1) {
    printf("%s: Trouble with call to init_pso.\n", __FUNCTION__);
    printf("FAIL");
    return -1;
  }

  fprintf(stderr, "%s: Allocating memory for particles.\n", __FUNCTION__);
  
  retval = alloc_ppack(&(ps.pp), ps.num_particles);
  if (retval == -1) {
    printf("%s: Trouble with call to alloc_ppack.\n", __FUNCTION__);
    printf("FAIL");
    return -1;
  }

  retval = alloc_threadprogress(&(ps.gr), ps.num_generations);
  if (retval == -1) {
    printf("%s: Trouble with call to alloc_threadprogress.\n", __FUNCTION__);
    printf("FAIL");
    return -1;
  }

  retval = alloc_wpstatus(&(ps.wp), ps.workpacks_pergeneration);
  if (retval == -1) {
    printf("%s: Trouble with call to alloc_wpstatus.\n", __FUNCTION__);
    printf("FAIL");
    return -1;
  }

  fprintf(stderr, "%s: Processing PSO algorithm.\n", __FUNCTION__);

  {

    fitnessfunc_extra ff_extra;

    fitnesspack *fpack;
    
    size_t len;

    fpack = &(ps.pp.fpack);
    
    ff_extra.matrix = NULL;
    ff_extra.inv_matrix = NULL;
    ff_extra.blo = ps.blo;
    ff_extra.bup = ps.bup;
    ff_extra.sf = 1.0;

    len = strlen(mfunc_str);

    fpack->ff_extra = &ff_extra;
    
    if (!strncmp(mfunc_str, "pyramid", len)) {
    
      fpack->fitness_func = pyramid_fitnessfunc;

    }

    if (!strncmp(mfunc_str, "bumps", len)) {
    
      fpack->fitness_func = bumps_fitnessfunc;

    }

    if (!strncmp(mfunc_str, "sinc", len)) {
    
      fpack->fitness_func = sinc_fitnessfunc;

    }

    if (!strncmp(mfunc_str, "torus", len)) {
    
      fpack->fitness_func = torus_fitnessfunc;

    }

    if (!strncmp(mfunc_str, "tube", len)) {
    
      fpack->fitness_func = tube_fitnessfunc;

    }

    if (!strncmp(mfunc_str, "stairs", len)) {
    
      fpack->fitness_func = stairs_fitnessfunc;

    }
    
    ps.pp.up = tuning_param(0.875, 2.0, 2.0);

    retval = setfull_computation(ps.packs_region, ps.workpacks_pergeneration);
    if (retval == -1) {
      fprintf(stderr, "%s: Region out of range?\n", __FUNCTION__);
      return -1;
    }
        
    retval = process_pso(&ps, fpack->fitness_func, fpack->ff_extra);
    if (retval == -1) {
      printf("%s: Trouble with processing of pso.\n", __FUNCTION__);
      printf("FAIL");
      return -1;
    }

  }
    
  printf("%s: Writing output values for csv.\n", __FUNCTION__);

  retval = psogen_csv(&(ps.gr), ps.num_generations, outcsv_fn);  
  if (retval == -1) {
    printf("%s: Trouble with call to psogen_csv.\n", __FUNCTION__);
    printf("FAIL");
    return -1;
  }
  
  {

    fprintf(stderr, "%s: Outputting asc result.\n", __FUNCTION__);
    
    retval = output_asc(&(ps.pp), ps.num_particles, outasc_fn);
    if (retval == -1) {
      printf("%s: Trouble with outputting of asc file.\n", __FUNCTION__);
      printf("FAIL");
      return -1;
    }

  }
    
  retval = close_pso(&ps);
  if (retval == -1) {
    printf("%s: Trouble with call to close_pso.\n", __FUNCTION__);
    printf("FAIL");
    return -1;
  }
  
  printf("SUCCESS");  
  
  return 0;

}
