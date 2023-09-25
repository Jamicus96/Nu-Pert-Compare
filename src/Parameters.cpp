#include <cmath>
#include <complex>

#include "Parameters.h"

// http://globalfit.astroparticles.es
double s12sq = 0.32;
double s13sq = 0.022;
double s23sq = 0.55;
double Dmsq21 = 7.5e-5; // eV^2
double Dmsq31 = 2.50e-3; // eV^2
double delta = -0.4 * M_PI; // chosen to avoid any accidental cancellations

double t12, t13, t23;
double Dmsq32, Dmsqee, eps;
double c12, s12, c12sq, s212, c212;
double c13, s13, c13sq, s213, c213;
double c23, s23, c23sq, s223, c223;
double cd, sd;
double mo_sign;
std::complex<double> eid;

double a0_vac, a1_vac, H_ee, Y_ee, R_H_em, I_H_em, R_Y_em, I_Y_em;

void Recalc_Parameters()
{
	c12sq = 1 - s12sq;
	c13sq = 1 - s13sq;
	c23sq = 1 - s23sq;

	s12 = sqrt(s12sq);
	s13 = sqrt(s13sq);
	s23 = sqrt(s23sq);

	c12 = sqrt(c12sq);
	c13 = sqrt(c13sq);
	c23 = sqrt(c23sq);

	t12 = acos(c12);
	t13 = acos(c13);
	t23 = acos(c23);

	s212 = sin(2 * t12);
	s213 = sin(2 * t13);
	s223 = sin(2 * t23);

	c212 = cos(2 * t12);
	c213 = cos(2 * t13);
	c223 = cos(2 * t23);

	Dmsqee = Dmsq31 - s12sq * Dmsq21;
	Dmsq32 = Dmsq31 - Dmsq21;
	eps = Dmsq21 / Dmsqee;

	cd = cos(delta);
	sd = sin(delta);
	eid = exp(std::complex<double>(0, delta));

	mo_sign = Dmsq31 > 0 ? 1 : -1;

	// Extra
	JP_Prob_Constants();
}

/**
 * @brief Computes constants needed for JP function.
 */
void JP_Prob_Constants() {

    // Compute vacuum constants common to all flavours
    a0_vac = (Dmsq21*Dmsq21*Dmsq21 + Dmsq31*Dmsq31*Dmsq31) / 27.0 - (Dmsq21*Dmsq21 * Dmsq31 + Dmsq21 * Dmsq31*Dmsq31) / 18.0;
    a1_vac = (Dmsq21*Dmsq21 + Dmsq31*Dmsq31 - Dmsq21 * Dmsq31) / 9.0;
    H_ee = Dmsq21 * (s12*s12 * c13*c13 - 1.0 / 3.0) + Dmsq31 * (s13*s13 - 1.0 / 3.0);
    Y_ee = (Dmsq21*Dmsq21 * (s12*s12 * c13*c13 - 1.0 / 3.0) + Dmsq31*Dmsq31 * (s13*s13 - 1.0 / 3.0) + 2.0 * Dmsq21 * Dmsq31 * (c12*c12 * c13*c13 - 1.0 / 3.0)) / 3.0;

    // Compute extra mu->e constants
    R_H_em = Dmsq21 * s12 * c13 * (c12 * c23 - s12 * s23 * s13 * cd) + Dmsq31 * s13 * s23 * c13 * cd;
    I_H_em = Dmsq21 * s12*s12 * s13 * s23 * c13 * sd - Dmsq31 * s13 * s23 * c13 * sd;

    R_Y_em = (s12 * c13 * (c12 * c23 - s12 * s13 * s23 * cd) * Dmsq21*Dmsq21
			+ s13 * s23 * c13 * cd * Dmsq31*Dmsq31
			- 2.0 * c12 * c13 * (s12 * c23 + s13 * s23 * c12 * cd) * Dmsq21 * Dmsq31) / 3.0;

    I_Y_em = (s12*s12 * s13 * s23 * c13 * sd * Dmsq21*Dmsq21
			- s13 * s23 * c13 * sd * Dmsq31*Dmsq31
			+ 2.0 * s13 * s23 * c12*c12 * c13 * sd * Dmsq21 * Dmsq31) / 3.0;
}

double eVsqkm_to_GeV = 1e-9 / 1.97327e-7 * 1e3;
double YerhoE2a = 1.52e-4;

