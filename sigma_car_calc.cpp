#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "fdtd3d.h"

void sig_car_calc(double ***sigma_car, double ***sigma_real, double **R2_1, double **invR1_2)
{

    Eigen::MatrixXd r21(3, 3);
    Eigen::MatrixXd inv_r12(3, 3);
    Eigen::MatrixXd sigma_car_mat(3, 3);
    Eigen::MatrixXd sigma_real_mat(3, 3);
    Eigen::MatrixXd A(3, 3);

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            r21(i, j) = R2_1[i][j];
            inv_r12(i, j) = invR1_2[i][j];
        }
    }

    for(int i = 0; i < ion_L; i++){

        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                sigma_real_mat(j, k) = sigma_real[i][j][k];
            }
        }
        
        sigma_car_mat = r21*sigma_real_mat*inv_r12;

        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                sigma_car[i][j][k] = sigma_car_mat(j, k);
            }
        }

    }
}