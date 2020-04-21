#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "fdtd3d.h"

void make_rot_mat(double **R2_1, double **invR1_2, double b_th, double b_phi)
{   
    Eigen::MatrixXd r1(3, 3);
    Eigen::MatrixXd r2(3, 3);
    Eigen::MatrixXd r21(3, 3);
    Eigen::MatrixXd inv_r12(3, 3);
    Eigen::MatrixXd inv_r1(3, 3);
    Eigen::MatrixXd inv_r2(3, 3);
    Eigen::MatrixXd unit(3, 3);

    r1(0, 0) = std::cos(b_th*M_PI/180.0);
    r1(0, 1) = 0.0;
    r1(0, 2) = std::sin(b_th*M_PI/180.0);
    r1(1, 0) = 0.0;
    r1(1, 1) = 1.0;
    r1(1, 2) = 0.0;
    r1(2, 0) = -std::sin(b_th*M_PI/180.0);
    r1(2, 1) = 0.0;
    r1(2, 2) = std::cos(b_th*M_PI/180.0);

    r2(0, 0) = std::cos(b_phi*M_PI/180.0);
    r2(0, 1) = -std::sin(b_phi*M_PI/180.0);
    r2(0, 2) = 0.0;
    r2(1, 0) = std::sin(b_phi*M_PI/180.0);
    r2(1, 1) = std::cos(b_phi*M_PI/180.0);
    r2(1, 2) = 0.0;
    r2(2, 0) = 0.0;
    r2(2, 1) = 0.0;
    r2(2, 2) = 1.0;

    //unit_matrix//
    unit = Eigen::MatrixXd::Identity(3, 3);

    inv_r1 = r1.fullPivLu().solve(unit);
    inv_r2 = r2.fullPivLu().solve(unit);

    r21 = r2*r1;
    inv_r12 = inv_r1*inv_r2;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            R2_1[i][j] = r21(i, j);
            invR1_2[i][j] = inv_r12(i, j);
        }
    }
    
}