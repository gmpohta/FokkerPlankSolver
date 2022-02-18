#ifndef __MAIN_H__
#define __MAIN_H__

#include <windows.h>
#include "solver.h"

#ifdef BUILD_DLL
    #define DLL_EXPORT __declspec(dllexport)
#else
    #define DLL_EXPORT __declspec(dllimport)
#endif

Solver* p_sol=NULL;

extern "C" void DLL_EXPORT get_grids(double *E_out, double *mu_out,  double *z_out, double *B_loop,double* n_loop);
extern "C" void DLL_EXPORT solve(double* time, double* S_in, double *J_out, double *f_out, double *n_fast_out);
extern "C" void DLL_EXPORT create_solver(long const *general_param, double const *grid_param, double const *inf_param, double const *loop_param, double const * turb_param);
extern "C" void DLL_EXPORT delete_solver();
#endif // __MAIN_H__
