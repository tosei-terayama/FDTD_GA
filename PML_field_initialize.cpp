#include <iostream>
#include "fdtd3d.h"

void PML_field_initialize(
    double ****Dr_theta1, double ****Dr_theta2, double ****Dr_phi,
    double ****Dtheta_phi, double ****Dtheta_r,
    double ****Dphi_r, double ****Dphi_theta,
    double ****Hr_theta1, double ****Hr_theta2, double ****Hr_phi,
    double ****Htheta_phi, double ****Htheta_r,
    double ****Hphi_r, double ****Hphi_theta)
{
    for(int i = 0; i <= 1; i++){
        //D components in PML(Theta direction)// 
        Dr_theta1[i] = memory_allocate3d(Nr, L, Nphi - 1, 0.0);
        Dr_theta2[i] = memory_allocate3d(Nr, L, Nphi - 1, 0.0);
        Dr_phi[i] = memory_allocate3d(Nr, L, Nphi - 1, 0.0);
        Dtheta_phi[i] = memory_allocate3d(Nr, L, Nphi - 1, 0.0);
        Dtheta_r[i] = memory_allocate3d(Nr, L, Nphi - 1, 0.0);
        Dphi_r[i] = memory_allocate3d(Nr, L, Nphi, 0.0);
        Dphi_theta[i] = memory_allocate3d(Nr, L, Nphi, 0.0);

        //H compornents in PML(Theta direction)//
        Hr_theta1[i] = memory_allocate3d(Nr + 1, L, Nphi, 0.0);
        Hr_theta2[i] = memory_allocate3d(Nr + 1, L, Nphi, 0.0);
        Hr_phi[i] = memory_allocate3d(Nr + 1, L, Nphi, 0.0);
        Htheta_phi[i] = memory_allocate3d(Nr, L + 1, Nphi, 0.0);
        Htheta_r[i] = memory_allocate3d(Nr, L + 1, Nphi, 0.0);
        Hphi_r[i] = memory_allocate3d(Nr, L, Nphi + 1, 0.0);
        Hphi_theta[i] = memory_allocate3d(Nr, L, Nphi + 1, 0.0);
    }

    //PML region (Phi direction)//
    for(int i = 2; i <= 3; i++){
        //D components in PML(Phi direction)//
        Dr_theta1[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L, 0.0);
        Dr_theta2[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L, 0.0);
        Dr_phi[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L, 0.0);
        Dtheta_phi[i] = memory_allocate3d(Nr, Ntheta - 2*L, Nphi - 1, 0.0);
        Dtheta_r[i] = memory_allocate3d(Nr, Ntheta - 2*L, Nphi - 1, 0.0);
        Dphi_r[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L, 0.0);
        Dphi_theta[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L, 0.0);

        //H components in PML(Phi direction)//
        Hr_theta1[i] = memory_allocate3d(Nr + 1, Ntheta - 2*L, L, 0.0);
        Hr_theta2[i] = memory_allocate3d(Nr + 1, Ntheta - 2*L, L, 0.0);
        Hr_phi[i] = memory_allocate3d(Nr + 1, Ntheta - 2*L, L, 0.0);
        Htheta_phi[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L, 0.0);
        Htheta_r[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L, 0.0);
        Hphi_r[i] = memory_allocate3d(Nr, Ntheta - 2*L, L + 1, 0.0);
        Hphi_theta[i] = memory_allocate3d(Nr, Ntheta - 2*L, L + 1, 0.0);
    }
        
}