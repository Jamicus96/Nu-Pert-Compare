#ifndef JP_H
#define JP_H

// From https://arxiv.org/abs/2309.06900

namespace JP {
void JP_Prob_Constants(double& a0_vac, double& a1_vac, double& H_ee, double& Y_ee, double& R_H_em, double& I_H_em, double& R_Y_em, double& I_Y_em, double& PHASE);
void Pmue(std::vector<double>& a, std::vector<double>& L, std::vector<double>& E, std::vector<double>& P);
double Pmue(double a, double L, double E);
}

#endif
