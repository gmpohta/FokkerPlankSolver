#ifndef CONST_H_INCLUDED
#define CONST_H_INCLUDED
const double c=3e8;  // Speed of light, m/s
const double qe=1.6e-19;
const double me=9.1e-31;
const double eps0 = 8.85e-12;
const double PI = 3.141592653589793;
const double kbol=1.38e-23;  //J/K
const double r0_e=1.0/4.0/PI/eps0*qe*qe/me/c/c; ///Classical electron radius
const double lambda0_e=1.0/4.0/PI/r0_e/r0_e; /// use in coulonb diffusion for electrons
const double Rau=1.49597870700e11; /// Astronomical unit, m

#endif
