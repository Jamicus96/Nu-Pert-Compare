#include <cmath>
#include <vector>

#include "JP.h"
#include "Parameters.h"

#define square(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

namespace JP {

/**
 * @brief Computes constants needed for JP function.
 */
void JP_Prob_Constants(double& a0_vac, double& a1_vac, double& H_ee, double& Y_ee, double& R_H_em, double& I_H_em, double& R_Y_em, double& I_Y_em, double& PHASE) {

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

    // Extra
    PHASE = 2.0 * M_PI / 3.0;
}

double Pmue(double a, double L, double E)
{

    // Compute vacuum values first, independent of a, L or E
    double a0_vac, a1_vac, H_ee, Y_ee, R_H_em, I_H_em, R_Y_em, I_Y_em, PHASE;
    JP_Prob_Constants(a0_vac, a1_vac, H_ee, Y_ee, R_H_em, I_H_em, R_Y_em, I_Y_em, PHASE);

    // Compute all the matter-corrected quantities and then the transition probability
	double a0, a1, sqrt_a1, eigen[3], R_X[3], I_X[3], arcCos, L4E, Theta_10, Theta_20, Theta_21, denom;

    // Make  matter corrections (a0 and a1 have nothing to do with a)
    a /= 3.0;
    a0 = a0_vac + 1.5 * (Y_ee * a + H_ee * square(a)) + cube(a); // a0
    a1 = a1_vac + H_ee * a + square(a); // a1

    // Get eigenvalues of H, and constants X
    sqrt_a1 = sqrt(a1);
    arcCos = acos(a0 / (sqrt_a1 * a1)) / 3.0;
    sqrt_a1 *= 2.0;

    for (unsigned int i = 0; i < 3; ++i) {
        eigen[i] = sqrt_a1 * cos(arcCos - PHASE * i);
        denom = square(eigen[i]) - a1;
        R_X[i] = ((eigen[i] + a) * R_H_em + R_Y_em) / denom;
        I_X[i] = ((eigen[i] + a) * I_H_em + I_Y_em) / denom;
    }

    // Compute Theta constants (not mixing angles)
    L4E = eVsqkm_to_GeV * L / (4 * E);
    Theta_10 = (eigen[1] - eigen[0]) * L4E;
    Theta_20 = (eigen[2] - eigen[0]) * L4E;
    Theta_21 = (eigen[2] - eigen[1]) * L4E;

    // Compute probability
    return 2.0 * (- 2.0 * ((R_X[1]*R_X[0] + I_X[1]*I_X[0]) * square(sin(Theta_10))
                            + (R_X[2]*R_X[0] + I_X[2]*I_X[0]) * square(sin(Theta_20))
                            + (R_X[2]*R_X[1] + I_X[2]*I_X[1]) * square(sin(Theta_21)))
                            + ((I_X[1]*R_X[0] - R_X[1]*I_X[0]) * sin(2.0 * Theta_10)
                            + (I_X[2]*R_X[0] - R_X[2]*I_X[0]) * sin(2.0 * Theta_20)
                            + (I_X[2]*R_X[1] - R_X[2]*I_X[1]) * sin(2.0 * Theta_21))) / 9.0;
}
}

