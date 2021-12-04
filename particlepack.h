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

#ifndef PARTICLEPACK_H
#define PARTICLEPACK_H

#include "vector.h"

#include "mini_gxkit.h"

#include "fitnesspack.h"

typedef struct {

  double accel_c;
  double c1, c2;
  
} update_param;

typedef struct {

  vec3d *vcur;
  point3d_t *xcur;

  point3d_t *pbest;
  double *fitness;

  point3d_t gbest;
  double gbesterr;
  
  update_param up;

  fitnesspack fpack;
  
} particlepack;

int alloc_ppack(particlepack *pp, long int num_particles);

#endif
