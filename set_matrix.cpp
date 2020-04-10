#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "fdtd3d.h"

void set_matrix(
    std::complex <double> zj, double ***Cmat, double ***Fmat,
    double* Nh, double* Ny){
    double omg_c = E_Q*B_abs/E_M;

    for(int ir = Nr - ion_L; ir < Nr; ir++){
        int i = ir - (Nr - ion_L);
        double Alt = ir*delta_r;
        double omg_p = E_Q*std::sqrt(Ny[i]/E_M/EPS0);
        std::complex <double> omg = omega - zj * Nh[i];
        std::complex <double> diag_comp = omega/(omg_c*omg_c - omg*omg);
        std::complex <double> offd_comp = zj * omg_c / (omg_c*omg_c - omg*omg);
        std::complex <double> coef = zj * EPS0 * omg_p*omg_p;
        std::cout << i*delta_r*1e-3 << " " << omg_p/2.0/M_PI*1e-3 << std::endl;

        Eigen::Matrix3d Sigma = Eigen::Matrix3d::Zero(3,3);
        Sigma(0, 0) = real( coef*diag_comp );
        Sigma(1, 1) = real( coef*diag_comp );
        Sigma(0, 1) = real( coef*offd_comp );
        Sigma(1, 0) = real( -1.0*coef*offd_comp );
        Sigma(2, 2) = real( -1.0*coef/omg );

        Eigen::Matrix3d A = 
                EPS0/Dt*Eigen::Matrix3d::Identity(3,3) + 0.5*Sigma;
        Eigen::Matrix3d B =
                EPS0/Dt*Eigen::Matrix3d::Identity(3,3) - 0.5*Sigma;
        Eigen::Matrix3d C = A.inverse()*B;
        Eigen::Matrix3d F = 1.0/Dt*A.inverse();

        for(int m = 0; m < 3; m++){
            for(int n = 0; n < 3; n++){
                Cmat[i][m][n] = C(m, n);
                Fmat[i][m][n] = F(m, n);
            }
        }
    }

}