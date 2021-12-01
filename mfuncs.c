
#include <math.h>

#include "mini_gxkit.h"

#include "mfuncs.h"

#include "sincfunc.h"

point3d_t math3d_sinc(double x, double y) {

  point3d_t pnta;

  double hyp;

  double sf;
  
  sf = 15.0;

  hyp = sf * sqrt(x * x + y * y);
  
  pnta = (point3d_t) { .x = x, .y = y, .z = sinc(hyp) };
  
  return pnta;

}

point3d_t math3d_tube(double x, double y) {

  point3d_t pnta;

  pnta = (point3d_t) { .x = x, .y = y, .z = 1.0 / (15.0 * (x * x + y * y)) };
  
  return pnta;

}

point3d_t math3d_torus(double x, double y) {

  point3d_t pnta;

  double value1;

  value1 = (0.6 - pow(x * x + y * y, 0.5));

  pnta = (point3d_t) { .x = x, .y = y, .z = pow(0.4 * 0.4 - (value1 * value1), 0.5)  };
  
  return pnta;

}

point3d_t math3d_bumps(double x, double y) {

  point3d_t pnta;

  pnta = (point3d_t) { .x = x, .y = y, .z = sin(5.0 * x) * cos(5.0 * y) / 5.0 };
  
  return pnta;

}

point3d_t math3d_cone(double x, double y) {

  point3d_t pnta;

  pnta = (point3d_t) { .x = x, .y = y, .z = pow(x * x + y * y, 0.5) };
  
  return pnta;
  
}

point3d_t math3d_pyramid(double x, double y) {

  point3d_t pnta;

  pnta = (point3d_t) { .x = x, .y = y, .z = 1.0 - fabs(x + y) - fabs(y - x) };
  
  return pnta;

}

point3d_t math3d_ripple(double x, double y) {

  point3d_t pnta;

  pnta = (point3d_t) { .x = x, .y = y, .z = sin(10.0 * (x * x + y * y)) / 10.0 };
  
  return pnta;

}
