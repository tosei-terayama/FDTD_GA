#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void E_update(
  double**** E_r, double**** E_theta, double**** E_phi,
  double**** D_r, double**** D_theta, double**** D_phi,
  int New, int Old, double***** Cmat, double***** Fmat)
{
  double ***newD_r = D_r[New], ***oldD_r = D_r[Old];
  double ***newD_th = D_theta[New], ***oldD_th = D_theta[Old];
  double ***newD_ph = D_phi[New], ***oldD_ph = D_phi[Old];
  double maxE(0.0);
  double interpol_Er(0.0), interpol_Eth(0.0), interpol_Eph(0.0);
  double interpol_nDr(0.0), interpol_nDth(0.0), interpol_nDph(0.0),
  interpol_oDr(0.0), interpol_oDth(0.0), interpol_oDph(0.0);
  constexpr int Ir { 0 }, Ith{ 1 }, Iph{ 2 };

  for(int i = 0; i < Nr - ion_L; i++){
    for(int j = 1; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
        E_r[New][i][j][k] = E_r[Old][i][j][k] + 
          (newD_r[i][j][k] - oldD_r[i][j][k])/EPS0;

        if(maxE < std::abs(E_r[New][i][j][k])){
          maxE = std::abs(E_r[New][i][j][k]);
        }

      }
    }
  }

  for(int i = Nr - ion_L; i < Nr; i++){
    int m = i - (Nr - ion_L);
    for(int j = 1; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
          
          interpol_Eth = (
            E_theta[Old][i  ][j][k] + E_theta[Old][i  ][j-1][k] + 
            E_theta[Old][i+1][j][k] + E_theta[Old][i+1][j-1][k])/4.0;
          interpol_nDth = (
            newD_th[i  ][j][k] + newD_th[i  ][j-1][k] + 
            newD_th[i+1][j][k] + newD_th[i+1][j-1][k])/4.0;
          interpol_oDth = (
            oldD_th[i  ][j][k] + oldD_th[i  ][j-1][k] + 
            oldD_th[i+1][j][k] + oldD_th[i+1][j-1][k])/4.0;

          interpol_Eph = (
            E_phi[Old][i  ][j][k] + E_phi[Old][i  ][j][k-1] + 
            E_phi[Old][i+1][j][k] + E_phi[Old][i+1][j][k-1])/4.0;
          interpol_nDph = (
            newD_ph[i  ][j][k] + newD_ph[i  ][j][k-1] + 
            newD_ph[i+1][j][k] + newD_ph[i+1][j][k-1])/4.0;
          interpol_oDph = (
            oldD_ph[i  ][j][k] + oldD_ph[i  ][j][k-1] + 
            oldD_ph[i+1][j][k] + oldD_ph[i+1][j][k-1])/4.0;
          
          E_r[New][i][j][k] = 
              Cmat[m][j][k][Ir][Ir] * E_r[Old][i][j][k] + 
              Cmat[m][j][k][Ir][Ith] * interpol_Eth + 
              Cmat[m][j][k][Ir][Iph] * interpol_Eph + 
              Fmat[m][j][k][Ir][Ir] * (newD_r[i][j][k] - oldD_r[i][j][k]) + 
              Fmat[m][j][k][Ir][Ith] * (interpol_nDth - interpol_oDth) + 
              Fmat[m][j][k][Ir][Iph] * (interpol_nDph - interpol_oDph);

          if(maxE < std::abs(E_r[New][i][j][k])){
            maxE = std::abs(E_r[New][i][j][k]);
          }

      }
    }
  }

  for(int i = 1; i < Nr - ion_L; i++){
    for(int j = 0; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
          E_theta[New][i][j][k] = E_theta[Old][i][j][k] + 
            (newD_th[i][j][k] - oldD_th[i][j][k])/EPS0;

          if(maxE < std::abs(E_theta[New][i][j][k])){
            maxE = std::abs(E_theta[New][i][j][k]);
          }
      }
    }
  }

  for(int i = Nr - ion_L; i < Nr; i++){
    int m = i - (Nr - ion_L);
    for(int j = 0; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){

          interpol_Er = (
            E_r[Old][i][j  ][k] + E_r[Old][i-1][j  ][k] + 
            E_r[Old][i][j+1][k] + E_r[Old][i-1][j+1][k])/4.0;
          interpol_nDr = (
            newD_r[i][j  ][k] + newD_r[i-1][j  ][k] +
            newD_r[i][j+1][k] + newD_r[i-1][j+1][k])/4.0;
          interpol_oDr = (
            oldD_r[i][j  ][k] + oldD_r[i-1][j  ][k] + 
            oldD_r[i][j+1][k] + oldD_r[i-1][j+1][k])/4.0;
          
          interpol_Eph = (
            E_phi[Old][i][j  ][k] + E_phi[Old][i][j  ][k-1] + 
            E_phi[Old][i][j+1][k] + E_phi[Old][i][j+1][k-1])/4.0;
          interpol_Eph = (
            newD_ph[i][j  ][k] + newD_ph[i][j  ][k-1] +
            newD_ph[i][j+1][k] + newD_ph[i][j+1][k-1])/4.0;
          interpol_Eth = (
            oldD_ph[i][j  ][k] + oldD_ph[i][j  ][k-1] + 
            oldD_ph[i][j+1][k] + oldD_ph[i][j+1][k-1])/4.0;

          E_theta[New][i][j][k] = 
              Cmat[m][j][k][Ith][Ir] * interpol_Er +
              Cmat[m][j][k][Ith][Ith] * E_theta[Old][i][j][k] +
              Cmat[m][j][k][Ith][Iph] * interpol_Eph +
              Fmat[m][j][k][Ith][Ir] * (interpol_nDr - interpol_oDr) +
              Fmat[m][j][k][Ith][Ith] * (newD_th[i][j][k] - oldD_th[i][j][k]) +
              Fmat[m][j][k][Ith][Iph] * (interpol_nDph - interpol_oDph);
          
          if(maxE < std::abs(E_theta[New][i][j][k])){
            maxE = std::abs(E_theta[New][i][j][k]);
          }

      }
    } 
  }

  for(int i = 1; i < Nr - ion_L; i++){
    for(int j = 1; j < Ntheta; j++){
      for(int k = 0; k < Nphi; k++){
        E_phi[New][i][j][k] = E_phi[Old][i][j][k] + 
              (newD_ph[i][j][k] - oldD_ph[i][j][k])/EPS0;

        if(maxE < std::abs(E_phi[New][i][j][k])){
          maxE = std::abs(E_phi[New][i][j][k]);
        }

      }
    }
  }

  for(int i = Nr - ion_L; i < Nr; i++){
    int m = i - (Nr - ion_L);
    for(int j = 1; j < Ntheta; j++){
      for(int k = 0; k < Nphi; k++){

          interpol_Er = (
            E_r[Old][i][j][k  ] + E_r[Old][i-1][j][k  ] + 
            E_r[Old][i][j][k+1] + E_r[Old][i-1][j][k+1])/4.0;
          interpol_nDr = (
            newD_r[i][j][k  ] + newD_r[i-1][j][k  ] + 
            newD_r[i][j][k+1] + newD_r[i-1][j][k+1])/4.0;
          interpol_oDr = (
            oldD_r[i][j][k  ] + oldD_r[i-1][j][k  ] + 
            oldD_r[i][j][k+1] + oldD_r[i-1][j][k+1])/4.0;

          interpol_Eth = (
            E_theta[Old][i][j][k  ] + E_theta[Old][i][j-1][k  ] + 
            E_theta[Old][i][j][k+1] + E_theta[Old][i][j-1][k+1])/4.0;
          interpol_nDth = (
            newD_th[i][j][k  ] + newD_th[i][j-1][k  ] + 
            newD_th[i][j][k+1] + newD_th[i][j-1][k+1])/4.0;
          interpol_oDth = (
            oldD_th[i][j][k  ] + oldD_th[i][j-1][k  ] + 
            oldD_th[i][j][k+1] + oldD_th[i][j-1][k+1])/4.0;

          E_phi[New][i][j][k] =
                Cmat[m][j][k][Iph][Ir] * interpol_Er + 
                Cmat[m][j][k][Iph][Ith] * interpol_Eth + 
                Cmat[m][j][k][Iph][Iph] * E_phi[Old][i][j][k] + 
                Fmat[m][j][k][Iph][Ir] * (interpol_nDr - interpol_oDr) + 
                Fmat[m][j][k][Iph][Ith] * (interpol_nDth - interpol_oDth) +
                Fmat[m][j][k][Iph][Iph] * (newD_ph[i][j][k] - oldD_ph[i][j][k]);

          if(maxE < std::abs(E_phi[New][i][j][k])){
            maxE = std::abs(E_phi[New][i][j][k]);
          }

      }
    }
  }
  
  if(maxE > 1.0e12) exit(0);
  
}
