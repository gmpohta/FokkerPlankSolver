
#include "MyTypes.h"
#include "solver.h"
#include "const.h"
#include <vector>
#include <iostream>
double* tridiag_alg(double const * a,double const * b, double const * c, double const * r,unsigned long beg_ind,unsigned long end_ind);

/******Class Solver*****/
Solver::~Solver(){
    delete time;
    delete mu_;
    //Log grid
    //for E
    delete dEdk_;
    delete kE_;
    delete E_;
    //for z
    delete z_;
    delete dzdk_;
    //
    delete dlnB;
    delete B_loop;
    delete z_loop;
    delete n_background;
    delete conductivity;
    delete coulomb_log;
    //
    delete lamb_turb;
}

Solver::Solver(const long *general_param, const double *grid_param, const double *inf_param, const double *loop_param, const double * turb_param){
    ///general_param = [NE, Nmu, Nz_output, Nz_loop_param, switch_reverse_current, switch_turb]
    ///grid_param=[Emin,zmax,tmax,dt0,dt_max,gamma, r_gridE, r_gridz]
    ///inf_param=[n_inf, log_inf, conductivity_inf]
    ///loop_param=[*z_loop,*B_loop,*n_loop,*coulomb_log,*conductivity]
    ///turb_param=[*E_turb, *lamb_turb]

    nE=general_param[0]; //number of intervals
    nm=general_param[1];
    nz=general_param[2];
    nz_loop=general_param[3]; // number of point
    switch_reverse_current=general_param[4];
    switch_turb=general_param[5];

    define_turbulent_diff(turb_param);

    ///Solution area
    Emin=grid_param[0]*qe/me/c/c*1e6; ///input in MeV, convert to dimensionless
    zmax=grid_param[1];/// Length of the loop, input in m
    tmax=grid_param[2];
    dt=grid_param[3]/tmax;
    dt_max=grid_param[4]/tmax;
    dt_old=0.0;
    gamma=grid_param[5];
    ///Coefficient for the logarithmic grid in energy and z
    r_gridE=grid_param[6]*qe/me/c/c*1e6; ///input in MeV, convert to dimensionless
    r_gridz=grid_param[7]/zmax; //input in m, convert to dimensionless

    //arrays grid (dimensionless variables)
    time=new std::vector<double>;
    define_mu_grid();
    define_z_grid();
    define_E_grid();
    //
    n_inf=inf_param[0];///input in m-3
    log_inf=inf_param[1];
    conductivity_inf=inf_param[2];///input in Omh-1m-1
    //
    define_z_loop(0,loop_param);
    find_bondary_index();
    define_dlnB(1*nz_loop,loop_param);
    define_n_background(2*nz_loop,loop_param);
    define_coulomb_log(3*nz_loop,loop_param);
    define_conductivity(4*nz_loop,loop_param);

    //Set first moment in time
    nt=0;
    time->push_back(0.0);
    t_tmp=(*time)[nt]+dt;
}

void Solver::clear_step(){
    delete ftmp;
    delete Jtmp;
    ftmp=NULL;
    Jtmp=NULL;
}

void Solver::define_initial_f(double * S_in, double * f_in, double* J_in, double* n_fast_in){
    external_f=new Array3d(nE+1,nm+1,2*nz+1,f_in);
    external_J=new Array1d(2*nz+1,J_in);
    n_fast=new Array1d(2*nz+1,n_fast_in);
    ftmp=new Array3d(nE+1,nm+1,2*nz+1);
    Jtmp=new Array1d(2*nz+1);
    source=new Array3d(nE+1,nm+1,2*nz+1,S_in);

    for (long ii=0; ii<nE+1;ii++){
        for (long jj=0;jj<nm+1;jj++){
            for(long kk=0;kk<2*nz+1;kk++){
                (*external_f)(ii,jj,kk)*=me*c*c/qe/n_loop_min*1e-6; //input in m-3MeV-1, convert to dimensionless
            }
        }
    }
    for (long kk=0; kk<2*nz+1; kk+=1){(*external_J)(kk)*=1/n_loop_min/qe/c;} //input in A/m2, convert to dimensionless
}

void Solver::define_turbulent_diff(double const *turb_param){
    lamb_turb=new Array1d(nE+1);
    for (long ii=0; ii<nE+1;ii++){
        (*lamb_turb)(ii)=turb_param[ii];///input in m
    }
}

void Solver::define_z_loop(const long shift_ind, double const *loop_param){
    z_loop = new Array1d(nz_loop);
    double z_start=loop_param[0+shift_ind];
    double z_end= loop_param[nz_loop-1+shift_ind];

    for (long itt=0; itt<nz_loop; itt+=1){
        (*z_loop)(itt)=RIGHT_END_LOOP-(z_end-loop_param[itt+shift_ind])/(z_end-z_start);  //transform z=[0,zmax] in m to z=[-0.5,0.5](dimensionless)
    }
}

void Solver::define_n_background(const long shift_ind, const double* loop_param){
    ///input in m-3
    n_background  = new Array1d(int(2*nz+1));

    double *tmp= new double[nz_loop];
    for (long itt=0; itt<nz_loop; itt+=1){
        tmp[itt]=loop_param[itt+shift_ind];
    }
    Spline* sp_n = new Spline(z_loop->getp_data(), tmp, nz_loop);
    double* interpolated_n = sp_n->calcspln(z_,2*nz+1);

    double z_r=(*z_)(n_z_end-1);
    double n_r=interpolated_n[n_z_end-1];
    double z_L=(*z_)(n_z_start+1);
    double n_L=interpolated_n[n_z_start+1];
    double k_exp=700.0;

    for (unsigned long kk=0; kk<2*nz+1;kk++){
        if ((*z_)(kk)<=RIGHT_END_LOOP and (*z_)(kk)>=LEFT_END_LOOP){
            (*n_background)(kk)= interpolated_n[kk];
        }else if ((*z_)(kk)>RIGHT_END_LOOP){
            (*n_background)(kk)= n_inf-(n_inf-n_r)*exp(-k_exp*pow((*z_)(kk)-z_r,2));
        }else{
            (*n_background)(kk)= n_inf-(n_inf-n_L)*exp(-k_exp*pow((*z_)(kk)-z_L,2));
        }
    }

    delete tmp;
    delete sp_n;
    delete interpolated_n;

    //Find min loop density
    n_loop_min=(*n_background)(0);
    for (long kk=0; kk<2*nz+1; kk+=1){
       if ((*n_background)(kk)<n_loop_min) n_loop_min=(*n_background)(kk);
    }
}
void Solver::find_bondary_index(){
    //Finding indexes which correspond first points without the loop
    n_z_end=0;
    n_z_start=0;
    for (long kk=0; kk<2*nz+1; kk+=1){
        if ((*z_)(kk)<=RIGHT_END_LOOP and (*z_)(kk)>=LEFT_END_LOOP){
            if (n_z_end==0)
                n_z_start=kk;
            n_z_end+=1;
        }
    }
    n_z_end+=n_z_start;
    n_z_start-=1;
}
void Solver::define_coulomb_log(const long shift_ind, double const *loop_param){
    coulomb_log=  new Array1d(2*nz+1);

    double *tmp= new double[nz_loop];
    for (long itt=0; itt<nz_loop; itt+=1){
        tmp[itt]=loop_param[itt+shift_ind];
    }
    Spline* sp_log = new Spline(z_loop->getp_data(),tmp,nz_loop);
    double* interpolated_log = sp_log->calcspln(z_,2*nz+1);
    for (unsigned long kk=0; kk<2*nz+1;kk++){
        if ((*z_)(kk)<=RIGHT_END_LOOP and (*z_)(kk)>=LEFT_END_LOOP)
            (*coulomb_log)(kk)= interpolated_log[kk];
        else
            (*coulomb_log)(kk)= log_inf;
    }

    delete tmp;
    delete sp_log;
    delete interpolated_log;
}

void Solver::define_conductivity(const long shift_ind, double const *loop_param){
    ///input in Omh-1m-1
    conductivity = new Array1d(2*nz+1);

    double *tmp= new double[nz_loop];
    for (long itt=0; itt<nz_loop; itt+=1){
        tmp[itt]=loop_param[itt+shift_ind];
    }
    Spline* sp = new Spline(z_loop->getp_data(),tmp,nz_loop);
    double* interpolated = sp->calcspln(z_,2*nz+1);
    for (unsigned long kk=0; kk<2*nz+1;kk++){
        if ((*z_)(kk)<=RIGHT_END_LOOP and (*z_)(kk)>=LEFT_END_LOOP)
            (*conductivity)(kk)= interpolated[kk];
        else
            (*conductivity)(kk)= conductivity_inf;
    }

    delete tmp;
    delete sp;
    delete interpolated;
}
void Solver::define_kT_el(const long shift_ind, double const *loop_param){
    kT_el= new Array1d(2*nz+1);

    double *tmp= new double[nz_loop];
    for (long itt=0; itt<nz_loop; itt+=1){
        tmp[itt]=loop_param[itt+shift_ind];
    }
    Spline* sp = new Spline(z_loop->getp_data(),tmp,nz_loop);
    double* interpolated = sp->calcspln(z_,2*nz+1);
    for (unsigned long kk=0; kk<2*nz+1;kk++){
        if ((*z_)(kk)<=RIGHT_END_LOOP and (*z_)(kk)>=LEFT_END_LOOP)
            (*kT_el)(kk)= interpolated[kk];
        else
            (*kT_el)(kk)= kT_el_inf;
    }

    delete tmp;
    delete sp;
    delete interpolated;
}
void Solver::define_dlnB(const long shift_ind, const double *loop_param){
    dlnB = new Array1d(2*nz+1);
    B_loop=new Array1d(nz_loop);
    for (long itt=0; itt<nz_loop; itt+=1){
        (*B_loop)(itt)=loop_param[itt+shift_ind];
    }
    //Calculating dlnB/dz  on a non-uniform grid using the Newton interpolate polynomial
    Array1d* tmp_dlnB = new Array1d(nz_loop);
    double u01=0;
    double u12=0;
    double u012=0;
    for (long kk=0; kk<nz_loop-2; kk+=1){
        u01=((*B_loop)(kk)-(*B_loop)(kk+1))/((*z_loop)(kk)-(*z_loop)(kk+1));
        u12=((*B_loop)(kk+1)-(*B_loop)(kk+2))/((*z_loop)(kk+1)-(*z_loop)(kk+2));
        u012=(u01-u12)/((*z_loop)(kk)-(*z_loop)(kk+2));
        (*tmp_dlnB)(kk)=(u01+((*z_loop)(kk)-(*z_loop)(kk+1))*u012)/(*B_loop)(kk);
    }
    u01=((*B_loop)(nz_loop-1)-(*B_loop)(nz_loop-2))/((*z_loop)(nz_loop-1)-(*z_loop)(nz_loop-2));
    u12=((*B_loop)(nz_loop-2)-(*B_loop)(nz_loop-3))/((*z_loop)(nz_loop-2)-(*z_loop)(nz_loop-3));
    u012=(u01-u12)/((*z_loop)(nz_loop-1)-(*z_loop)(nz_loop-3));

    (*tmp_dlnB)(nz_loop-2)=(u01+((*z_loop)(nz_loop-2)-(*z_loop)(nz_loop-1))*u012)/(*B_loop)(nz_loop-2);
    (*tmp_dlnB)(nz_loop-1)=(u01+((*z_loop)(nz_loop-1)-(*z_loop)(nz_loop-2))*u012)/(*B_loop)(nz_loop-1);

    //interpolate dlnB/dz
    Spline* sp_dlnB = new Spline(z_loop->getp_data(),tmp_dlnB->getp_data(),nz_loop);
    double* interpolated_dlnB = sp_dlnB->calcspln(z_,2*nz+1);
    for (unsigned long kk=0; kk<2*nz+1;kk++){
        if ((*z_)(kk)<=RIGHT_END_LOOP and (*z_)(kk)>=LEFT_END_LOOP)
            (*dlnB)(kk)= interpolated_dlnB[kk];
        else
            (*dlnB)(kk)= 0;
    }

    ///linear extrapolate to points z(n_z_start) and z(n_z_end)
    ///We cannot just assume (*dlnB)(n_z_start)=calcspln((*z_)(n_z_start)) becose (*z_)(n_z_start)<-1 (outside the loop boundary)
    ///Therefore, we do a linear extrapolation by 2 points: by z = -0.5 and z = (*z_)(n_z_start+1)
    double first=LEFT_END_LOOP;
    double second=(*z_)(n_z_start+1);
    double y_second=(*dlnB)(n_z_start+1);
    double y_first = sp_dlnB->calcspln(first);
    (*dlnB)(n_z_start)= y_first+(y_second-y_first)/(second-first)*((*z_)(n_z_start)-first);
    ///The same for right end to loop
    first=RIGHT_END_LOOP;
    second=(*z_)(n_z_end-1);
    y_first = sp_dlnB->calcspln(first);
    y_second=(*dlnB)(n_z_end-1);
    (*dlnB)(n_z_end)= y_first+(y_second-y_first)/(second-first)*((*z_)(n_z_end)-first);

    delete sp_dlnB;
    delete interpolated_dlnB;
    delete tmp_dlnB;
}

void Solver::define_mu_grid(){
    mu_ = new Array1d(nm+1);

    double delt=2.0/nm;
    for (long j=0; j<nm+1;j++){(*mu_)(j)=j*delt-1;};
    dmu=2.0/nm;
}

void Solver::define_E_grid(){
    dEdk_ = new Array1d(nE+1);
    kE_ = new Array1d(nE+1);
    E_ = new Array1d(nE+1);

    double kminE=1.0-exp(-Emin/r_gridE);
    dkE=(1.0-kminE)/nE;

    for (long ii=0;ii<nE+1;ii++){
        (*kE_)(ii)=dkE*ii+kminE;
        (*E_)(ii)=-r_gridE*log(1.0-(*kE_)(ii));
        (*dEdk_)(ii)=r_gridE/(1.0-(*kE_)(ii));
    }
    (*E_)(nE)=INFINITY;
}

void Solver::define_z_grid(){
    dzdk_ = new Array1d(int(2*nz+1));
    z_ = new Array1d(int(2*nz+1));
    double kz;
    dkz=2.0/2.0/nz; //kz in range [-1,1]
    for (long kk=0;kk<nz+1;kk++) {
        kz=dkz*kk-1;
        (*z_)(kk)=r_gridz*log(1.0+kz);
        (*dzdk_)(kk)=r_gridz/(1.0+kz);
    }
    for (long kk=nz+1;kk<2*nz+1;kk++){
        kz=dkz*kk-1;
        (*z_)(kk)=-r_gridz*log(1.0-kz);
        (*dzdk_)(kk)=r_gridz/(1.0-kz);
    }
    (*z_)(2*nz)=INFINITY;
    (*z_)(0)=INFINITY;
}

void Solver::f_step(){
    //Note ct+1 is next time moment
    for (long i=0; i<nE+1;i++){
        for (long j=0;j<nm+1;j++){
            (*ftmp)(i,j,0)=0.0;
            (*ftmp)(i,j,2*nz)=0.0;
        }
    }
    for (long k=0; k<2*nz+1;k++){
        for (long j=0;j<nm+1;j++){
            (*ftmp)(nE,j,k)=0.0;
        }
    }
    Array3d& f_old=*external_f;
    Array1d& J_old=*external_J;
    double S_old;
    Array3d *f1=new Array3d(nE+1,nm+1,2*nz+1);
    Array3d *f2=new Array3d(nE+1,nm+1,2*nz+1);
    for (int jj=0; jj<nm+1; jj++)
        for (int kk=0; kk<2*nz+1;kk++)
            (*f1)(nE,jj,kk)=0.0;
    double a1,a2,aJ2,aJ3,a3,a4,a5;
    double op_mu2,op_mu1,op_E,op_z,sou;
    for (long ii=nE-1;ii>=0;ii--){ //reversed(range(nE))
        for (long kk=1; kk<2*nz;kk++){ //
            a3=k3(ii,kk);
            a5=k5(ii,kk);
            for (long jj=1;jj<nm;jj++){
                //We need to take it at the middle moment of time, but we dont know the end moment, so we take it at the initial moment
                S_old=(*source)(ii,jj,kk)*me*c*c/qe/n_loop_min*1e-6*tmax; //input in m-3MeV-1s-1 convert to dimensionless
                a1=k1(ii,jj);
                a2=k2(ii,jj,kk);
                aJ2=kJ2(ii,jj,kk)*J_old(kk);
                aJ3=kJ3(ii,jj,kk)*J_old(kk);
                a4=k4(ii,jj,kk);

                op_mu2=a4/dmu/dmu*(f_old(ii,jj+1,kk)-2.0*f_old(ii,jj,kk)+f_old(ii,jj-1,kk));
                op_E=(a3+aJ3)/dkE/(*dEdk_)(ii)*(f_old(ii+1,jj,kk)-f_old(ii,jj,kk));
                op_mu1=(a2+aJ2)/2.0/dmu*(f_old(ii,jj+1,kk)-f_old(ii,jj-1,kk));
                op_z=a1/2.0/dkz/(*dzdk_)(kk)*(f_old(ii,jj,kk+1)-f_old(ii,jj,kk-1));
                sou=S_old+op_mu1+op_mu2+op_E+op_z+a5*f_old(ii,jj,kk);
                (*f1)(ii,jj,kk)=(sou+(a3+aJ3)*dt/dkE/(*dEdk_)(ii)*(*f1)(ii+1,jj,kk))/(1.0+(a3+aJ3)*dt/dkE/(*dEdk_)(ii)-a5*dt);
            }
        }
    }
    double *ay=new double[nm+1];
    double *by=new double[nm+1];
    double *cy=new double[nm+1];
    double *ry=new double[nm+1];

    for (long ii=0;ii<nE;ii++){
        for (long kk=1; kk<2*nz;kk++){
            for (long jj=1;jj<nm;jj++){
                a2=k2(ii,jj,kk);
                aJ2=kJ2(ii,jj,kk)*J_old(kk);
                a4=k4(ii,jj,kk);

                ay[jj]=(a2+aJ2)*dt/2.0/dmu-a4*dt/dmu/dmu;
                cy[jj]=-(a2+aJ2)*dt/2.0/dmu-a4*dt/dmu/dmu;
                ry[jj]=(*f1)(ii,jj,kk);
                by[jj]=1.0+2.0*a4*dt/dmu/dmu;
            }
            by[1]=1.0+(k2(ii,1,kk)+kJ2(ii,1,kk)*J_old(kk))*dt/dmu;
            cy[1]=-(k2(ii,1,kk)+kJ2(ii,1,kk)*J_old(kk))*dt/dmu;
            by[nm-1]=1.0-(k2(ii,nm-1,kk)+kJ2(ii,nm-1,kk)*J_old(kk))*dt/dmu;
            ay[nm-1]=(k2(ii,nm-1,kk)+kJ2(ii,nm-1,kk)*J_old(kk))*dt/dmu;

            double* tmp=tridiag_alg(ay,by,cy,ry,1,nm);
            for (long jj=1;jj<nm;jj++){(*f2)(ii,jj,kk)=tmp[jj-1];};
            (*f2)(ii,0,kk)=2.0*(*f2)(ii,1,kk)-(*f2)(ii,2,kk);
            (*f2)(ii,nm,kk)=2.0*(*f2)(ii,nm-1,kk)-(*f2)(ii,nm-2,kk);
            delete tmp;
        }
    }

    double* az=new double[2*nz+1];
    double* bz=new double[2*nz+1];
    double* cz=new double[2*nz+1];
    double* rz=new double[2*nz+1];
    for (long ii=0;ii<nE;ii++){
        for (long jj=0; jj<nm+1;jj++){
            a1=k1(ii,jj);
            for (long kk=1; kk<2*nz;kk++){
                az[kk]=a1*dt/2.0/dkz/(*dzdk_)(kk);
                cz[kk]=-a1*dt/2.0/dkz/(*dzdk_)(kk);
                rz[kk]=(*f2)(ii,jj,kk)*dt-a1*dt/2.0/dkz/(*dzdk_)(kk)*(f_old(ii,jj,kk+1)-f_old(ii,jj,kk-1))+f_old(ii,jj,kk);
                bz[kk]=1.0;
            }
            double* tmp=tridiag_alg(az,bz,cz,rz,1,2*nz);
            for (size_t kk=1;kk<2*nz;kk++){(*ftmp)(ii,jj,kk)=tmp[kk-1];};
            delete tmp;
        }
    }

    delete ay;
    delete by;
    delete cy;
    delete ry;
    delete az;
    delete bz;
    delete cz;
    delete rz;
    delete f1;
    delete f2;
}

double Solver::integral_for_reverse_current(unsigned long step,unsigned long kk,unsigned ct){
    double out=0.0;
    double bet=0.0;
    double hx=(*kE_)(step)-(*kE_)(0);
    double hy=(*mu_)(step)-(*mu_)(0);
    Array3d& f_old=*external_f;

    for (long ii=step; ii<nE; ii+=step){
        bet=speed(ii);
        out+=hy/2.0*hx*((*ftmp)(ii,0,kk)+f_old(ii,0,kk))/2.0*bet*(*mu_)(0)*(*dEdk_)(ii);
        out+=hy/2.0*hx*((*ftmp)(ii,nm,kk)+f_old(ii,nm,kk))/2.0*bet*(*mu_)(nm)*(*dEdk_)(ii);
    }
    for (long jj=step;jj<nm; jj+=step){
        bet=speed(0);
        out+=hx/2.0*hy*((*ftmp)(0,jj,kk)+f_old(0,jj,kk))/2.0*bet*(*mu_)(jj)*(*dEdk_)(0);
        for (long ii=step;ii<nE;ii+=step){
            bet=speed(ii);
            out+=hy*hx*((*ftmp)(ii,jj,kk)+f_old(ii,jj,kk))/2.0*bet*(*mu_)(jj)*(*dEdk_)(ii);
        }
        //There should be one more line here for i=nE, but f(nE,jj,kk,ct+1) =0 so there is no line
    }
    return out;
}

void Solver::reccurrent_calc(double* Ui,long const recc_number,long const N){
    double R=0.0;
    for (long i=recc_number; i<N; i+=1){
        R=(Ui[i]-Ui[i-1])/(pow(2.0,2*recc_number)-1.0);
        if (R!=R) R=0.0;
        Ui[i]+=R;
    }
}

bool Solver::add_layer(){
    nt+=1;
    time->push_back(t_tmp);

     for (long ii=0; ii<nE+1;ii++){
        for (long jj=0;jj<nm+1;jj++){
            for(long kk=0;kk<2*nz+1;kk++){
                (*external_f)(ii,jj,kk)=(*ftmp)(ii,jj,kk);
            }
        }
    }
    if (switch_reverse_current){
        for (long kk=0; kk<2*nz+1; kk+=1){(*external_J)(kk)= (*Jtmp)(kk);}//copy from tmp to external
    }
    return 1;
}

bool Solver::change_step_time_or_break(){
    Array3d& f_old=*external_f;
    double fmax=0.0;
    double df_max=0.0;
    double df=0.0;
    bool out=0;
    /*
    If max(fold-fcurr)/fold<=gamma then go to the next time step.
    if the previous and current steps are not equal, then we do not decrease the step, if the previous and current steps are equal, then we increase the step by  a factor of 1.5
    If max(fold-fcurr)/fold>gamma then decrease the time step by the factor of 1.5 times and try to calculate again
    */

    for(int ii=0;ii<nE+1;ii++) {
        for(int jj=0;jj<nm+1;jj++) {
            for(int kk=0;kk<2*nz+1;kk++) {
                if ((*z_)(kk)>=LEFT_END_LOOP and (*z_)(kk)<=RIGHT_END_LOOP and std::abs(f_old(ii,jj,kk))>fmax){
                    fmax=std::abs(f_old(ii,jj,kk));
                }
            }
        }
    }

    for(int ii=0;ii<nE+1;ii++) {
        for(int jj=0;jj<nm+1;jj++) {
            for(int kk=0;kk<2*nz+1;kk++) {
                if ((*z_)(kk)>=LEFT_END_LOOP and (*z_)(kk)<=RIGHT_END_LOOP){
                    df=std::abs((f_old(ii,jj,kk)-(*ftmp)(ii,jj,kk))/fmax);
                    if (df==df and df!=1.0/0.0 and df>df_max){
                        df_max=df;
                    }
                }
            }
        }
    }

    if (df_max<=gamma){
        out=add_layer();
        if (dt!=dt_old){
            dt_old=dt;
            dt*=k_t_change;
        }
    }
    else{
        dt_old=dt;
        dt/=k_t_change;
    }
    if(dt>dt_max){dt=dt_max;}
    t_tmp=(*time)[nt]+dt;
    if (t_tmp>1.0){
        t_tmp=1.0;
        dt=dt_old;
    }
    return out;
}

void Solver::J_step(){
    int r=2;
    long N=0;// Number of grids (1 grid is spent for each refinement, and each grid is r times larger than the previous one)
    long N1=long(log(nE+1)/log(r))+1;
    long N2=long(log(nm+1)/log(r))+1;
    if (N1>N2)
        N=N2;
    else
        N=N1;

    for (long kk=0;kk<2*nz+1;kk+=1){
        double* Ui=new double[N];
        long step=1;
        //Compute integrals for each grid with steps: dm, dE; r*dm, r*dE; etc..
        for (long cc=N-1;cc>=0;cc-=1){
            Ui[cc]=integral_for_reverse_current(step,kk,nt);
            step*=r;
        }
        for (long cc=1;cc<N-2;cc+=1){//There is N-2 all possible clarifications
            reccurrent_calc(Ui,cc,N);
        }
        (*Jtmp)(kk)=-Ui[N-1];
        delete Ui;
    }
}

double Solver::integral_for_n_fast(long kk){
    double out=0.0;
    double hx=(*kE_)(1)-(*kE_)(0);
    double hy=(*mu_)(1)-(*mu_)(0);

    for (long ii=1; ii<nE; ii+=1){
        out+=hy/2.0*hx*(*external_f)(ii,0,kk)*(*dEdk_)(ii);
        out+=hy/2.0*hx*(*external_f)(ii,nm,kk)*(*dEdk_)(ii);
    }
    for (long jj=1;jj<nm; jj+=1){
        out+=hx/2.0*hy*(*external_f)(0,jj,kk)*(*dEdk_)(0);
        for (long ii=1;ii<nE;ii+=1){
            out+=hy*hx*(*external_f)(ii,jj,kk)*(*dEdk_)(ii);
        }
        //There should be one more line here for i=nE, but f(nE,jj,kk,ct+1) =0 so there is no line
    }
    return out;
}

void Solver::calc_n_fast(){
    for (long kk=0; kk<2*nz+1; kk+=1){
        (*n_fast)(kk)=integral_for_n_fast(kk)*n_loop_min; //m-3
    }
    for (long ii=0; ii<nE+1;ii++){
        for (long jj=0;jj<nm+1;jj++){
            for(long kk=0;kk<2*nz+1;kk++){
                (*external_f)(ii,jj,kk)*=n_loop_min/me/c/c*qe*1e6; //convert to m-3MeV-1
            }
        }
    }
    for (long kk=0; kk<2*nz+1; kk+=1){(*external_J)(kk)*=n_loop_min*qe*c;} //convert to A/m2
}

double* tridiag_alg(double const * a,double const * b, double const * c, double const * r,unsigned long beg_ind,unsigned long end_ind){
    //To solve three-diagonal system of linear equations
    unsigned long n=end_ind-beg_ind;
    double * alf = new double[n];
    double * beta=new double[n];
    double * out=new double[n];
    alf[0]=-c[beg_ind]/b[beg_ind];
    beta[0]=r[beg_ind]/b[beg_ind];
    for (unsigned long ii=1; ii<n-1; ii++){
        alf[ii]=-c[beg_ind+ii]/(b[beg_ind+ii]+alf[ii-1]*a[beg_ind+ii]);
        beta[ii]=(r[beg_ind+ii]-a[beg_ind+ii]*beta[ii-1])/(b[beg_ind+ii]+alf[ii-1]*a[beg_ind+ii]);
    }
    beta[n-1]=(r[n-1+beg_ind]-a[n-1+beg_ind]*beta[n-2])/(b[n-1+beg_ind]+alf[n-2]*a[n-1+beg_ind]);
    out[n-1]=beta[n-1];
    for (long ii=n-2; ii>=0; ii-=1){
        out[ii]=alf[ii]*out[ii+1]+beta[ii];
    }
    delete alf;
    delete beta;

    return out;
}

Spline::~Spline(){
    delete a_;
    delete b_;
    delete c_;
    delete d_;
    delete x_;
}

Spline::Spline(const double* const xi, const double* const yi,long n){
    n_ =n;
    a_ = new double[n_];
    b_ = new double[n_];
    c_ = new double[n_];
    d_ = new double[n_];
    x_ = new double[n_];

    for (long itt=0; itt<n_; itt+=1){
        a_[itt]=yi[itt];
        x_[itt]=xi[itt];
    }
    c_[0]=0.0;
    c_[n_-1]=0.0;

    double * A=new double[n_-1];
    Array1d * B=new Array1d(n_-1);
    double hi=0.0;
    double hi1=0.0;
    double d=0.0;
    double b=0.0;
    double z=0.0;

    (*B)(0)=c_[0];
    for (long itt=1; itt<n_-1; itt+=1){
        hi=xi[itt]-xi[itt-1];
        hi1=xi[itt+1]-xi[itt];
        b=2.0*(hi+hi1);
        d=3.0*((yi[itt+1]-yi[itt])/hi1-(yi[itt]-yi[itt-1])/hi);
        z=b+hi*A[itt-1];
        A[itt]=-hi/z;
        (*B)(itt)=(d-hi*(*B)(itt-1))/z;
    }
    for (long itt=n_-2; itt>=1; itt-=1){
        c_[itt]=(*B)(itt)+A[itt]*c_[itt+1];
    }
    for (long itt=0; itt<n_-1; itt+=1){
        hi1=xi[itt+1]-xi[itt];
        d_[itt]=(c_[itt+1]-c_[itt])/hi1/3.0;
        b_[itt]=(yi[itt+1]-yi[itt])/hi1-hi1/3.0*(2.0*c_[itt]+c_[itt+1]);
    }
    delete(A);
    delete(B);
}

double* Spline::calcspln(Array1d * const xin,long const nin){
    double* out= new double[nin];
    double a_tmp,b_tmp,c_tmp, d_tmp, x_tmp;

    long lbnd, ubnd, cmpval;
    for (int itt=0; itt<nin; itt+=1){
        double & xi=(*xin)(itt);
        if (xi!=xi){
            out[itt]=NAN;
        }
        else{
            if (xi<x_[0] or xi>x_[n_-1]){
                a_tmp=NAN;
                b_tmp=0;
                x_tmp=0;
                c_tmp=0;
                d_tmp=0;
            }
            else{
                lbnd=0;
                ubnd=n_;
                while (ubnd-lbnd>1){
                    cmpval=long((ubnd+lbnd)/2);
                    if (xi<=x_[cmpval])
                        ubnd=cmpval;
                    else
                        lbnd=cmpval;
                }
                a_tmp=a_[lbnd];
                b_tmp=b_[lbnd];
                x_tmp=x_[lbnd];
                c_tmp=c_[lbnd];
                d_tmp=d_[lbnd];
            }
            out[itt]=a_tmp+b_tmp*(xi-x_tmp)+(c_tmp+d_tmp*(xi-x_tmp))*(xi-x_tmp)*(xi-x_tmp);}
    }
    return out;
}

double Spline::calcspln(double const x_in){
    double out = 0;
    double a_tmp,b_tmp,c_tmp, d_tmp, x_tmp;

    long l_bnd, r_bnd, cmpval;
    if (x_in!=x_in){
        out=NAN;
    }
    else{
        if (x_in<x_[0] or x_in>x_[n_-1]){
            a_tmp=NAN;
            b_tmp=0;
            x_tmp=0;
            c_tmp=0;
            d_tmp=0;
        }
        else{
            l_bnd=0;
            r_bnd=n_;
            while (r_bnd-l_bnd>1){
                cmpval=long((r_bnd+l_bnd)/2);
                if (x_in<=x_[cmpval])
                    r_bnd=cmpval;
                else
                    l_bnd=cmpval;
            }
            a_tmp=a_[l_bnd];
            b_tmp=b_[l_bnd];
            x_tmp=x_[l_bnd];
            c_tmp=c_[l_bnd];
            d_tmp=d_[l_bnd];
        }
        out=a_tmp+b_tmp*(x_in-x_tmp)+(c_tmp+d_tmp*(x_in-x_tmp))*(x_in-x_tmp)*(x_in-x_tmp);
    }

    return out;
}

