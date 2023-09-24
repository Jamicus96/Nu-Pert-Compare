#ifndef Parameters_H
#define Parameters_H

#include <complex>

// oscillation parameters
extern double t12, t13, t23, Dmsq21, Dmsq31, delta;
extern double Dmsq32, Dmsqee, eps;
extern double c12, s12, c12sq, s12sq, s212, c212;
extern double c13, s13, c13sq, s13sq, s213, c213;
extern double c23, s23, c23sq, s23sq, s223, c223;
extern double cd, sd;
extern double mo_sign;
extern std::complex<double> eid;
// extern double a0_vac, a1_vac, H_ee, Y_ee, R_H_em, I_H_em, R_Y_em, I_Y_em;

void Recalc_Parameters();
// void JP_Prob_Constants();

// unit conversions
extern double eVsqkm_to_GeV;
extern double YerhoE2a;

#endif
