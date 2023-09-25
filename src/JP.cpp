#include <cmath>
#include <vector>

#include "JP.h"
#include "Parameters.h"

#define square(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

namespace JP {
double Pmue(double a, double L, double E)
{
    // Compute all the matter-corrected quantities and then the transition probability
	double a0, a1, sqrt_a1, eigen[3], R_X[3], I_X[3], arcCos, L4E, Theta_10, Theta_20, Theta_21, denom, PHASE;

    // Make  matter corrections (a0 and a1 have nothing to do with a)
    a /= 3.0;
    a0 = a0_vac + 1.5 * (Y_ee * a + H_ee * square(a)) + cube(a); // a0
    a1 = a1_vac + H_ee * a + square(a); // a1

    // Get eigenvalues of H, and constants X
    sqrt_a1 = sqrt(a1);
    arcCos = acos(a0 / (sqrt_a1 * a1)) / 3.0;
    sqrt_a1 *= 2.0;
    PHASE = 2.0 * M_PI / 3.0;

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

