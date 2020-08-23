#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Core>

#include "fdtd3d.h"

void set_perturbation(perturbation P_info, double*** noise_Nh, double* Nh){

    Eigen::Vector3d R_c;
    Eigen::Vector3d R;
    Eigen::Vector3d R_d;

    double z_0{ dist(P_info.r0()) };
    double R_theta{ 0.0 };
    int lower_r = (int)Alt_lower_ionosphere/1.0e3;

    R_c(0) = z_0*std::sin(th(P_info.th0()))*std::cos(ph(P_info.phi0()));
    R_c(1) = z_0*std::sin(th(P_info.th0()))*std::sin(ph(P_info.phi0()));
    R_c(2) = z_0*std::cos(th(P_info.th0()));

    for(int i = 0; i <= ion_L; i++){
        double z{ dist(i + lower_r) };

        for(int j = 0; j <= Ntheta; j++){
            for(int k = 0; k <= Nphi; k++){
                R(0) = z*std::sin(th(j))*std::cos(ph(k));
                R(1) = z*std::sin(th(j))*std::sin(ph(k));
                R(2) = z*std::cos(th(j));

                R_d = (R/R.norm()) * R_c.norm();
                R_theta = std::acos( (R_c.dot(R_d))/(R_c.norm()*R_d.norm()) );

                // 情報落ち対策 //
                if( std::abs(((R_c.dot(R_d)/(R_c.norm()*R_d.norm())))) > 1.0 ) R_theta = std::acos(1.0);

                double d_h{ z_0*R_theta };

                double enhance = P_info.alpha()*std::exp(- (std::pow(d_h, 2.0)/2.0/std::pow(P_info.sig_h(), 2.0)))
                                            *std::exp(- (std::pow(z - z_0, 2.0)/2.0/std::pow(P_info.sig_r(), 2.0)));

                noise_Nh[i][j][k] = Nh[i] + Nh[i] * enhance;

                std::cout << j << " " << k << std::endl;
            }
        }
    }
    std::cout << "loop OK" << std::endl;

}