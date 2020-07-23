#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Core>

#include "fdtd3d.h"

void set_perturbation(perturbation P_info, double*** dis_Nh, double* Nh){
    
    double p_r0 = P_info.r0();
    double p_th0 = P_info.th0();
    double p_phi0 = P_info.phi0();
    double Alpha = P_info.alpha();
    double Sigma_r = P_info.sig_r();
    double Sigma_h = P_info.sig_h();

    int r0 = (int)p_r0 - (int)Alt_lower_ionosphere/1.0e3;
    int th0 = (int)p_th0;
    int phi0 = (int)p_phi0;

    Eigen::Vector3d ion_r0;
    Eigen::Vector3d unit_r0;
    Eigen::Vector3d ion_r;
    Eigen::Vector3d unit_r;

    ion_r0(0) = dist(p_r0)*std::sin(th(p_th0))*std::cos(ph(p_phi0));
    ion_r0(1) = dist(p_r0)*std::sin(th(p_th0))*std::sin(ph(p_phi0));
    ion_r0(2) = dist(p_r0)*std::cos(th(p_th0));

    unit_r0 = ion_r0/ion_r0.norm();
    double ion_index(0.0);
    double horiz_l(0.0);
    double unit_th(0.0);

    int i, j, k;
    for(i = r0 - P_info.range_r(); i <= r0 + P_info.range_r(); i++){
        for(j = th0 - P_info.range_th(); j <= th0 + P_info.range_th(); j++){
            for(k = phi0 - P_info.range_phi(); k <= phi0 + P_info.range_phi(); k++){

                ion_index = (double)i + Alt_lower_ionosphere/1.0e3;
                ion_r(0) = dist(ion_index)*std::sin(th(double(j)))*std::cos(ph(double(k)));
                ion_r(1) = dist(ion_index)*std::sin(th(double(j)))*std::sin(ph(double(k)));
                ion_r(2) = dist(ion_index)*std::cos(th(double(j)));

                unit_r = ion_r/ion_r.norm();

                unit_th = std::acos(unit_r0.dot(unit_r));

                //情報落ち対策//
                if(std::abs(unit_r0.dot(unit_r)) > 1) unit_th = std::acos(1.0);

                horiz_l = ion_r0(0)*unit_th;
                dis_Nh[i][j][k] = Alpha * std::exp(-std::pow(dist(p_r0) - dist(ion_index), 2.0)/2.0/std::pow(Sigma_r,2.0))
                                * std::exp(-std::pow(horiz_l, 2.0)/2.0/std::pow(Sigma_h, 2.0));
            }
        }
        std::cout << i << " " << dis_Nh[i][(int)p_th0][(int)p_phi0] << std::endl;
    }
    std::cout << "check dist" << std::endl;

    for(int i = 0; i < ion_L; i++){
        for(int j = 0; j < Ntheta; j++){
            for(int k = 0; k < Nphi; k++){

                if( ((i < r0 + P_info.range_r())&&(i > r0 - P_info.range_r())) &&
                ((j < th0 + P_info.range_th())&&(j > th0 - P_info.range_th())) &&
                ((k < phi0 + P_info.range_phi())&&(k > phi0 - P_info.range_phi())) ){
                    dis_Nh[i][j][k] = Nh[i]*dis_Nh[i][j][k];
                }
                else dis_Nh[i][j][k] = Nh[i];
            }
        }
    }

}