#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "fdtd3d.h"

void coordinate_trans(double ****Cmat_r, double ****Fmat_r, double ****Cmat_th, double ****Fmat_th,
                      double ****Cmat_phi, double ****Fmat_phi, double ***sig_car, double ***sig_car_r)
{
    double **mat = memory_allocate2d(3, 3, 0.0);
    double **mat_r = memory_allocate2d(3, 3, 0.0);

    Eigen::MatrixXd R1(3, 3);
    Eigen::MatrixXd R2(3, 3);
    Eigen::MatrixXd A(3, 3);
    Eigen::MatrixXd invA(3, 3);
    Eigen::MatrixXd B(3, 3);
    Eigen::MatrixXd C(3, 3);
    Eigen::MatrixXd F(3, 3);
    Eigen::MatrixXd sig_car_mat(3, 3);
    Eigen::MatrixXd sig_car_mat_r(3, 3);
    Eigen::MatrixXd sig_sph_mat(3, 3);

    //Make Unit Matrix//
    Eigen::MatrixXd unit(3, 3);
    unit = Eigen::MatrixXd::Identity(3, 3);

    double theta(0.0), phi(0.0);

    for(int i = 0; i < ion_L; i++){

        /*for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                mat[j][k] = sig_car[i][j][k];
                mat_r[j][k] = sig_car_r[i][j][k];
            }
        }*/

        for(int m = 0; m < 3; m++){
            for(int n = 0; n < 3; n++){
                sig_car_mat(m, n) = sig_car[i][m][n];
                sig_car_mat_r(m, n) = sig_car_r[i][m][n];
            }
        }

        for(int j = 1; j < Ntheta; j++){
            theta = th(j);
            for(int k = 1; k < Nphi; k++){
                phi = k*delta_phi;

                R1(0, 0) = std::sin(theta)*std::cos(phi);
                R1(0, 1) = std::sin(theta)*std::sin(phi);
                R1(0, 2) = std::cos(theta);
                R1(1, 0) = std::cos(theta)*std::cos(phi);
                R1(1, 1) = std::cos(theta)*std::sin(phi);
                R1(1, 2) = -std::sin(theta);
                R1(2, 0) = -std::sin(phi);
                R1(2, 1) = std::cos(phi);
                R1(2, 2) = 0.0;

                R2 = R1.transpose();

                sig_sph_mat = R1*sig_car_mat_r*R2;

                A = (EPS0/Dt)*unit + sig_sph_mat/2.0;
                B = (EPS0/Dt)*unit - sig_sph_mat/2.0;

                invA = A.fullPivLu().solve(unit);

                C = invA*B;
                F = invA/Dt;

                for(int m = 0; m < 9; m++){
                    
                    Cmat_r[i][j][k][m] = C(m/3, m%3);
                    Fmat_r[i][j][k][m] = F(m/3, m%3);

                }
                
            }
        }

        for(int j = 0; j < Ntheta; j++){
            theta = th(j + 0.5);
            for(int k = 1; k < Nphi; k++){
                phi = k*delta_phi;

                R1(0, 0) = std::sin(theta)*std::cos(phi);
                R1(0, 1) = std::sin(theta)*std::sin(phi);
                R1(0, 2) = std::cos(theta);
                R1(1, 0) = std::cos(theta)*std::cos(phi);
                R1(1, 1) = std::cos(theta)*std::sin(phi);
                R1(1, 2) = -std::sin(theta);
                R1(2, 0) = -std::sin(phi);
                R1(2, 1) = std::cos(phi);
                R1(2, 2) = 0.0;

                R2 = R1.transpose();

                sig_sph_mat = R1*sig_car_mat_r*R2;

                A = (EPS0/Dt)*unit + sig_sph_mat/2.0;
                B = (EPS0/Dt)*unit - sig_sph_mat/2.0;

                invA = A.fullPivLu().solve(unit);

                C = invA*B;
                F = invA/Dt;

                for(int m = 0; m < 9; m++){
                    
                    Cmat_th[i][j][k][m] = C(m/3, m%3);
                    Fmat_th[i][j][k][m] = F(m/3, m%3);

                }

            }
        }

        for(int j = 1; j < Ntheta; j++){
            theta = th(j);
            for(int k = 0; k < Nphi; k++){
                phi = (k + 0.5)*delta_phi;

                R1(0, 0) = std::sin(theta)*std::cos(phi);
                R1(0, 1) = std::sin(theta)*std::sin(phi);
                R1(0, 2) = std::cos(theta);
                R1(1, 0) = std::cos(theta)*std::cos(phi);
                R1(1, 1) = std::cos(theta)*std::sin(phi);
                R1(1, 2) = -std::sin(theta);
                R1(2, 0) = -std::sin(phi);
                R1(2, 1) = std::cos(phi);
                R1(2, 2) = 0.0;

                R2 = R1.transpose();

                sig_sph_mat = R1*sig_car_mat_r*R2;

                A = (EPS0/Dt)*unit + sig_sph_mat/2.0;
                B = (EPS0/Dt)*unit - sig_sph_mat/2.0;

                invA = A.fullPivLu().solve(unit);

                C = invA*B;
                F = invA/Dt;

                for(int m = 0; m < 9; m++){
                    
                    Cmat_phi[i][j][k][m] = C(m/3, m%3);
                    Fmat_phi[i][j][k][m] = F(m/3, m%3);

                }

            }
        }
    }

}
