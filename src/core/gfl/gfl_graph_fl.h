/****************************************************************************
 *  Taken from the GFL package by Tansey and Scott                          *
 *  Modified by Yannick Spill on 05/30/2017                                 *
 ****************************************************************************/

#ifndef GRAPH_FL_H
#define GRAPH_FL_H

#ifdef __cplusplus
extern "C" {
#endif
  
#include <math.h>
#include "gfl_tf.h"
#include "gfl_utils.h"
#include <omp.h>

#define VARYING_PENALTY_DELAY 50


int graph_fused_lasso (int n, double *y,
                        int ntrails, unsigned *trails, unsigned *breakpoints,
                        double lam, double *alpha, double inflate,
                        int maxsteps, double converge,
                        double *beta);

int graph_fused_lasso_weight (int n, double *y, double *w,
                        int ntrails, unsigned *trails, unsigned *breakpoints,
                        double lam, double *alpha, double inflate,
                        int maxsteps, double converge,
                        double *beta);

int graph_fused_lasso_warm (int n, double *y,
                        int ntrails, unsigned *trails, unsigned *breakpoints,
                        double lam, double *alpha, double inflate,
                        int maxsteps, double converge,
                        double *beta, double *z, double *u);

int graph_fused_lasso_weight_warm (int n, double *y, double *w,
                        int ntrails, unsigned *trails, unsigned *breakpoints,
                        double lam, double *alpha, double inflate,
                        int maxsteps, double converge,
                        double *beta, double *z, double *u);


void update_beta(int n, double *y, double *z, double *u, int *nzmap, int *zmap, double alpha, double *beta);
void update_beta_weight(int n, double *y, double *w, double *z, double *u, int *nzmap, int *zmap, double alpha, double *beta);
void update_z(int ntrails, unsigned *trails, unsigned *breakpoints, double *beta, double *u, double lam, double *ybuf, double *wbuf, double *tf_dp_buf, double *z);
void update_u(int n, double *beta, double *z, int *zmap, int *nzmap, double *u);
double primal_resnorm(int n, double *beta, double *z, int *nzmap, int *zmap);
double dual_resnorm(int nz, double *z, double *zold, double alpha);

#ifdef __cplusplus
}
#endif
  
#endif




