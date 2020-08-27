#include <iostream>
#include "fdtd3d.h"


void delete_PML(
    double ****Dr_theta1, double ****Dr_theta2, double ****Dr_phi,
    double ****Dtheta_phi, double ****Dtheta_r,
    double ****Dphi_r, double ****Dphi_theta,
    double ****Hr_theta1, double ****Hr_theta2, double ****Hr_phi,
    double ****Htheta_phi, double ****Htheta_r,
    double ****Hphi_r, double ****Hphi_theta)
{

    for(int i = 0; i <= 1; i++){
        delete_3d(Dr_theta1[i], Nr, L);
        delete_3d(Dr_theta2[i], Nr, L);
        delete_3d(Dr_phi[i], Nr, L);
        delete_3d(Dtheta_phi[i], Nr, L);
        delete_3d(Dtheta_r[i], Nr, L);
        delete_3d(Dphi_r[i], Nr, L);
        delete_3d(Dphi_theta[i], Nr, L);

        delete_3d(Hr_theta1[i], Nr + 1, L);
        delete_3d(Hr_theta2[i], Nr + 1, L);
        delete_3d(Hr_phi[i], Nr + 1, L);
        delete_3d(Htheta_phi[i], Nr, L + 1);
        delete_3d(Htheta_r[i], Nr, L + 1);
        delete_3d(Hphi_r[i], Nr, L);
        delete_3d(Hphi_theta[i], Nr, L);
    }
    
    for(int i = 2; i <= 3; i++){
        delete_3d(Dr_theta1[i], Nr, Ntheta - 2*L - 1);
        delete_3d(Dr_theta2[i], Nr, Ntheta - 2*L - 1);
        delete_3d(Dr_phi[i], Nr, Ntheta - 2*L - 1);
        delete_3d(Dtheta_phi[i], Nr, Ntheta - 2*L);
        delete_3d(Dtheta_r[i], Nr, Ntheta - 2*L);
        delete_3d(Dphi_r[i], Nr, Ntheta - 2*L - 1);
        delete_3d(Dphi_theta[i], Nr, Ntheta - 2*L - 1);

        delete_3d(Hr_theta1[i], Nr + 1, Ntheta - 2*L);
        delete_3d(Hr_theta2[i], Nr + 1, Ntheta - 2*L);
        delete_3d(Hr_phi[i], Nr + 1, Ntheta - 2*L);
        delete_3d(Htheta_phi[i], Nr, Ntheta - 2*L - 1);
        delete_3d(Htheta_r[i], Nr, Ntheta - 2*L - 1);
        delete_3d(Hphi_r[i], Nr, Ntheta - 2*L);
        delete_3d(Hphi_theta[i], Nr, Ntheta - 2*L);
    }

    delete[] Dr_theta1;
    delete[] Dr_theta2;
    delete[] Dr_phi;
    delete[] Dtheta_phi;
    delete[] Dtheta_r;
    delete[] Dphi_r;
    delete[] Dphi_theta;

    delete[] Hr_theta1;
    delete[] Hr_theta2;
    delete[] Hr_phi;
    delete[] Htheta_phi;
    delete[] Htheta_r;
    delete[] Hphi_r;
    delete[] Hphi_theta;
}