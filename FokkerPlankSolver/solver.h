#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED
#include<math.h>
#include "MyTypes.h"
#include "const.h"
#include <vector>

class Solver{
public:
    unsigned long nE,nm,nz,nt;
    bool switch_turb; //on/off turbulence
    bool switch_reverse_current; //on/off reverse current
    Solver(const long *general_param, const double  *grid_param, const double  *inf_param,
           const double *loop_param, const double *turb_param);
    //Pointers to array of double
    double get_time()  {return (*time)[nt]*tmax;}; //s
    void get_E(double* E_out){for (size_t ii=0; ii<nE+1; ii+=1){E_out[ii]=(*E_)(ii)*me*c*c/qe*1e-6;}}; //output in MeV
    void get_z(double* z_out){for (size_t kk=0; kk<2*nz+1; kk+=1){z_out[kk]=((*z_)(kk)-LEFT_END_LOOP)*zmax;}}; //output in m
    void get_mu(double* mu_out){for (size_t jj=0; jj<nm+1; jj+=1){mu_out[jj]=(*mu_)(jj);}};
    void calc_n_fast();
    void get_dlnB(double* dlnB_out) {for (size_t kk=0;kk<2*nz+1; kk+=1){dlnB_out[kk]=(*dlnB)(kk)/zmax;}}; //output in m-1
    void get_n_background(double* n_loop_out) {for (size_t kk=0;kk<2*nz+1; kk+=1){n_loop_out[kk]=(*n_background)(kk);}}; //output in m-3
    //Time step for f and J
    void f_step();
    void J_step();
    bool change_step_time_or_break();
    //
    void define_initial_f(double * S_in, double* f_in, double* J_in,double* n_fast_in);
    //
    void clear_step();
    ~Solver();

private:
    //An arrays of pointers to solutions
    Array3d* external_f=NULL;
    Array3d* source=NULL;
    Array1d* external_J=NULL;
    double const RIGHT_END_LOOP = 0.5; //dimensionless coordinate
    double const LEFT_END_LOOP = -0.5; //dimensionless coordinate
    double n_loop_min; //min density to which we normalize
    //functions defining grids
    void define_mu_grid();
    void define_E_grid();
    void define_z_grid();
    //
    double dmu,dkE,dkz;
    //
    void find_bondary_index();
    //boundary index, index witch correspond start-1 and end+1 loop coordinates
    long n_z_end,n_z_start;
    // number of points along the loop for input params
    unsigned long nz_loop;
    //Turbulence
    void define_turbulent_diff(double const *turb_param);
    Array1d* lamb_turb=NULL;
    //
    double tmax;
    double dt_old,dt,dt_max;
    const double k_t_change=1.5;
    double gamma;// For variable time step
    //Diffusion
    double lambda_e(unsigned kk){return lambda0_e/(*n_background)(kk)/(*coulomb_log)(kk);};
    //for current
    double integral_for_reverse_current(unsigned long step,unsigned long kk,unsigned ct);
    void reccurrent_calc(double* Ui,long const recc_number, long const N);
    //Density of non thermal electrons
    Array1d* n_fast=NULL;
    double integral_for_n_fast(long kk);
    //To define parameters distributions along the loop
    double n_inf,kT_el_inf, conductivity_inf,log_inf;
    void define_z_loop(const long shift_ind, double const *loop_param);
    void define_n_background(const long shift_ind, const double  *loop_param);
    void define_coulomb_log(const long shift_ind, double const *loop_param);
    void define_conductivity(const long shift_ind, double const *loop_param);
    void define_kT_el(const long shift_ind, double const *loop_param);
    void define_dlnB(const long shift_ind, double const *loop_param);
    //
    double zmax, Emin;
    //For logarifmical grids
    double r_gridz,r_gridE;
    //Grid for t (dimensionless var)
    std::vector<double>* time;
    double t_tmp;
    //For mu
    Array1d* mu_=NULL;
    //For E
    Array1d* dEdk_=NULL;
    Array1d* E_=NULL;
    Array1d* kE_=NULL;
    //For s
    Array1d* z_=NULL;
    Array1d* dzdk_=NULL;
    //Temporary array for the layer in time
    Array3d* ftmp=NULL;
    Array1d* Jtmp=NULL;
    bool add_layer();
    //Distribution magnetic field, density thermal plasma in calculations area
    Array1d* kT_el=NULL;
    Array1d* conductivity =NULL;
    Array1d* coulomb_log=NULL;
    Array1d* dlnB=NULL;
    Array1d* n_background =NULL;
    //Distribution input magnetic field and density
    Array1d* B_loop=NULL;
    Array1d* z_loop=NULL;
    //Dimensionless speed of electrons
    double speed(unsigned ii) {return pow((*E_)(ii)*((*E_)(ii)+2.0),0.5)/((*E_)(ii)+1.0);};
    //
    double k1(unsigned ii,unsigned jj){return -c*tmax/zmax*(*mu_)(jj)*speed(ii);};
    double k2(unsigned ii,unsigned jj,unsigned kk){double b =speed(ii);
        return -2.0*c*tmax*(*mu_)(jj)/lambda_e(kk)/pow(b,3)/pow((*E_)(ii)+1.0,2)+c*tmax/zmax*b/2.0*(1.0-(*mu_)(jj)*(*mu_)(jj))*(*dlnB)(kk)-switch_turb*c*tmax/(*lamb_turb)(ii)*b*(*mu_)(jj);};
    double k3(unsigned ii,unsigned kk){return c*tmax/speed(ii)/lambda_e(kk);};
    double k4(unsigned ii,unsigned jj,unsigned kk){double b =speed(ii);
        return c*tmax/lambda_e(kk)/pow(b,3)/pow((*E_)(ii)+1.0,2)*(1.0-(*mu_)(jj)*(*mu_)(jj))+switch_turb*c*tmax/(*lamb_turb)(ii)*b/2.0*(1.0-(*mu_)(jj)*(*mu_)(jj));};
    double k5(unsigned ii,unsigned kk){return -c*tmax/lambda_e(kk)/pow((*E_)(ii)*(*E_)(ii)+2.0*(*E_)(ii),1.5);};
    double kJ2(unsigned ii, unsigned jj,unsigned kk){
        return -qe*qe*n_loop_min*tmax/me/(*conductivity)(kk)*(1.0-(*mu_)(jj)*(*mu_)(jj))/speed(ii);};
    double kJ3(unsigned ii, unsigned jj,unsigned kk){double b =speed(ii);
        return -qe*qe*n_loop_min*tmax/me/(*conductivity)(kk)*(*mu_)(jj)*b;};
};

class Spline{
//Class spline interpolate
public:
    Spline(const double* const xi, const double* const yi,long n);
    double* calcspln(Array1d* const xin,long const n);
    double calcspln(double const x_in);
    ~Spline();
private:
    long n_;
    double* a_;
    double* b_;
    double* c_;
    double* d_;
    double* x_;
};
#endif



