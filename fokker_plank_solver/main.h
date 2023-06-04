#ifndef __MAIN_H__
#define __MAIN_H__

#include "solver.h"

extern "C" void get_grids(double *E_out, double *mu_out,  double *z_out, double *B_loop,double* n_loop);
extern "C" void solve(double* time, double* S_in, double *J_out, double *f_out, double *n_fast_out);
extern "C" void create_solver(long const *general_param, double const *grid_param, double const *inf_param, double const *loop_param, double const * turb_param);
extern "C" void delete_solver();
#endif // __MAIN_H__
