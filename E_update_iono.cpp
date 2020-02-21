#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "fdtd3d.h"

double E_update_iono(double **sigma_car, double Er, double Eth, double Eph,
                   double nDr, double nDth, double nDph,
                   double oDr, double oDth, double oDph,
                   int flag, double *Cmat, double *Fmat)
{
    double retE(0.0);

    Eigen::MatrixXd C(3, 3);
    Eigen::MatrixXd F(3, 3);
    Eigen::VectorXd o_E(3);
    Eigen::VectorXd n_E(3);
    Eigen::VectorXd o_D(3);
    Eigen::VectorXd n_D(3);

    //Make E vector//
    o_E(0) = Er;
    o_E(1) = Eth;
    o_E(2) = Eph;

    //Make D vetors//
    n_D(0) = nDr;
    n_D(1) = nDth;
    n_D(2) = nDph;

    o_D(0) = oDr;
    o_D(1) = oDth;
    o_D(2) = oDph;

    for(int m = 0; m < 9; m++){
        C(m/3, m%3) = Cmat[m];
        F(m/3, m%3) = Fmat[m];
    }

    n_E = C*o_E + F*(n_D - o_D);

    switch(flag){
        case 1:
        retE = n_E(0);
        break;

        case 2:
        retE = n_E(1);
        break;

        case 3:
        retE = n_E(2);
        break;

        default:
        std::cout << "flag = " << flag << std::endl << "Error!!" << std::endl;
        exit(0);
    }

    return retE;
}