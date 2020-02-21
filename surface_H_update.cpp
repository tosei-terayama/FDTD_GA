#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void surface_H(double **newEr, double **newEth, double **newEphi,
                double **Hth_r1, double **Hphi_r1,
                double z_real, double z_imag)
{
    double coeff_1 = (z_real/2.0) + (z_imag/Dt);
    double coeff_2 = (z_real/2.0) - (z_imag/Dt);
    double val_1, val_2, val_3, val_4;
    double sin_th1;

    double ri_1 = dist(0);
    double ri_2 = dist(0.5);
    double ri_3 = dist(1);

    val_2 = Dt/MU0/ri_2/delta_r;
    
    for(int j = L + 1; j <= Ntheta - L - 1; j++){
        sin_th1 = std::sin(th(j));
        val_1 = Dt/MU0/ri_2/sin_th1/delta_phi;
        for(int k = L; k <= Nphi - L - 1; k++){
            Hth_r1[j][k] = ((1.0 - ri_1*val_2*coeff_2)*Hth_r1[j][k] - val_1*(newEr[j][k + 1] - newEr[j][k])
             + ri_3*val_2*newEphi[j][k])/(1.0 + val_2*ri_1*coeff_1);
        }
    }

    val_3 = Dt/MU0/ri_2/delta_r;
    val_4 = Dt/MU0/ri_2/delta_theta;

    for(int j = L; j <= Ntheta - L - 1; j++){
        for(int k = L + 1; k <= Nphi - L - 1; k++){
            Hphi_r1[j][k] = ((1.0 - ri_1*val_3*coeff_2)*Hphi_r1[j][k] - ri_2*val_3*newEth[j][k]
             + val_4*(newEr[j + 1][k] - newEr[j][k]))/(1.0 + ri_1*val_3*coeff_1);
        }
    }

}
