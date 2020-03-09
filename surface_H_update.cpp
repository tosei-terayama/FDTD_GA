#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void surface_H_update(double **newEr, double **newEth, double **newEphi,
                double **Htheta, double **Hphi,
                double z_real, double z_imag)
{
    double coeff_1 = (z_real/2.0) + (z_imag/Dt);
    double coeff_2 = (z_real/2.0) - (z_imag/Dt);
    double val_1, val_2, val_3, val_4;
    double sin_th;

    double ri_1 = dist(0);
    double ri_2 = dist(0.5);
    double ri_3 = dist(1);

    val_2 = Dt/MU0/ri_2/delta_r;

    for(int j = 0; j < Ntheta + 1; j++){
        sin_th = std::sin(th(j));
        val_1 = Dt/MU0/ri_2/sin_th/delta_phi;
        for(int k = 0; k < Nphi; k++){
            Htheta[j][k] = ((1.0 - val_2*ri_1*coeff_2)*Htheta[j][k] - val_1*(newEr[j][k + 1] - newEr[j][k])
                        + val_2*ri_3*newEphi[j][k])/(1.0 + val_2*ri_1*coeff_1);
        }
    }

    val_3 = Dt/MU0/ri_2/delta_r;
    val_4 = Dt/MU0/ri_2/delta_theta;

    for(int j = 0; j < Ntheta; j++){
        for(int k = 0; k < Nphi + 1; k++){
            Hphi[j][k] = ((1.0 - val_3*ri_1*coeff_2)*Hphi[j][k] + val_4*(newEr[j + 1][k] - newEr[j][k])
                          - val_3*ri_3*newEth[j][k])/(1.0 + val_3*ri_1*coeff_1);
        }
    }

    /*for(int j  = 0; j < Ntheta + 1; j++){
        sin_th1 = std::sin(th(j));
        val_1 = Dt/MU0/ri_2/sin_th1/delta_phi;
        for(int k = 0; k < Nphi; k++){
            Htheta[j][k] = ((1.0 - val_2*ri_1*coeff_2)*Htheta[j][k] - val_1*(newEr[j][k + 1] - newEr[j][k])
                + ri_3*val_2*newEphi[j][k])/(1.0 + val_2*ri_1*coeff_1);
        }
    }

    val_3 = Dt/MU0/ri_2/delta_r;
    val_4 = Dt/MU0/ri_2/delta_theta;

    for(int j = L; j < Ntheta - L; j++){
        for(int k = L + 1; k < Nphi - L; k++){
            Hphi[j][k] = ((1.0 - ri_1*val_3*coeff_2)*Hphi[j][k] + val_4*(newEr[j + 1][k] - newEr[j][k])
                - ri_3*val_3*newEth[j][k])/(1.0 + val_3*ri_1*coeff_1);
        }
    }*/
}

