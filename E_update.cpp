#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void E_update(
  double**** E_r, double**** E_theta, double**** E_phi,
  double**** D_r, double**** D_theta, double**** D_phi,
  int New, int Old, double*** Cmat, double*** Fmat)
{
  double ***newD_r = D_r[New], ***oldD_r = D_r[Old];
  double ***newD_th = D_theta[New], ***oldD_th = D_theta[Old];
  double ***newD_ph = D_phi[New], ***oldD_ph = D_phi[Old];
  double theta(0.0), phi(0.0);
  double maxE(0.0);
  double interpol_Er(0.0), interpol_Eth(0.0), interpol_Eph(0.0);
  double interpol_nDr(0.0), interpol_nDth(0.0), interpol_nDph(0.0),
  interpol_oDr(0.0), interpol_oDth(0.0), interpol_oDph(0.0);
  int flag(0);
  int iono_flag(0);
  int I, J, K;
  constexpr int Ir { 0 }, Ith{ 1 }, Iph{ 2 };

  for(int i = 0; i < Nr - ion_L; i++){
    for(int j = 1; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
        E_r[New][i][j][k] = E_r[Old][i][j][k] + 
          (newD_r[i][j][k] - oldD_r[i][j][k])/EPS0;
        
        if(maxE < std::abs(E_r[New][i][j][k])){
          flag = 1;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_r[New][i][j][k]);
        }

      }
    }
  }

   iono_flag = 1;
  for(int i = Nr - ion_L; i < Nr; i++){
    int m = i - (Nr - ion_L);
    for(int j = 1; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
          phi = k*delta_phi;
          
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
              Cmat[m][Ir][Ir] * E_r[Old][i][j][k] + 
              Cmat[m][Ir][Ith] * interpol_Eth + 
              Cmat[m][Ir][Iph] * interpol_Eph + 
              Fmat[m][Ir][Ir] * (newD_r[i][j][k] - oldD_r[i][j][k]) + 
              Fmat[m][Ir][Ith] * (interpol_nDth - interpol_oDth) + 
              Fmat[m][Ir][Iph] * (interpol_nDph - interpol_oDph);
           
         /*E_r[NEW][i][j][k] = E_update_iono(sigma_car_r[i - Nr + ion_L], E_r[OLD][i][j][k], interpol_Eth, interpol_Eph,
         D_r[NEW][i][j][k], interpol_nDth, interpol_nDph, D_r[OLD][i][j][k], interpol_oDth, interpol_oDph,
         iono_flag, Cmat_r[i - Nr + ion_L][j][k], Fmat_r[i - Nr + ion_L][j][k]);*/

          if(maxE < std::abs(E_r[New][i][j][k])){
          flag = 1;
          I = i;
          J = j;
          K = k;
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
          flag = 2;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_theta[New][i][j][k]);
        }

      }
    }
  }

   iono_flag = 2;
  for(int i = Nr - ion_L; i < Nr; i++){
    int m = i - (Nr - ion_L);
    for(int j = 0; j < Ntheta; j++){
      //theta = th(j + 0.5);
      for(int k = 1; k < Nphi; k++){
          //phi = k*delta_phi;

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

          /*E_theta[NEW][i][j][k] = E_update_iono(sigma_car[i - Nr + ion_L], interpol_Er, E_theta[OLD][i][j][k], interpol_Eph,
          interpol_nDr, D_theta[NEW][i][j][k], interpol_nDph, interpol_oDr, D_theta[OLD][i][j][k], interpol_oDph,
          iono_flag, Cmat_th[i - Nr + ion_L][j][k], Fmat_th[i - Nr + ion_L][j][k]);*/

          E_theta[New][i][j][k] = 
              Cmat[m][Ith][Ir] * interpol_Er +
              Cmat[m][Ith][Ith] * E_theta[Old][i][j][k] +
              Cmat[m][Ith][Iph] * interpol_Eph +
              Fmat[m][Ith][Ir] * (interpol_nDr - interpol_oDr) +
              Fmat[m][Ith][Ith] * (newD_th[i][j][k] - oldD_th[i][j][k]) +
              Fmat[m][Ith][Iph] * (interpol_nDph - interpol_oDph);

          if(maxE < std::abs(E_theta[New][i][j][k])){
          flag = 2;
          I = i;
          J = j;
          K = k;
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
          flag = 3;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_phi[New][i][j][k]);
      }

      }
    }
  }

  iono_flag = 3;
  for(int i = Nr - ion_L; i < Nr; i++){
    int m = i - (Nr - ion_L);
    for(int j = 1; j < Ntheta; j++){
      //theta = th(j);
      for(int k = 0; k < Nphi; k++){
          //phi = (k + 0.5)*delta_phi;

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

          /*E_phi[NEW][i][j][k] = E_update_iono(sigma_car[i - Nr + ion_L], interpol_Er, interpol_Eth, E_phi[OLD][i][j][k],
          interpol_nDr, interpol_nDth, D_phi[NEW][i][j][k], interpol_oDr, interpol_oDth, D_phi[OLD][i][j][k],
          iono_flag, Cmat_phi[i - Nr + ion_L][j][k], Fmat_phi[i - Nr + ion_L][j][k]);*/

          E_phi[New][i][j][k] =
                Cmat[m][Iph][Ir] * interpol_Er + 
                Cmat[m][Iph][Ith] * interpol_Eth + 
                Cmat[m][Iph][Iph] * E_phi[Old][i][j][k] + 
                Fmat[m][Iph][Ir] * (interpol_nDr - interpol_oDr) + 
                Fmat[m][Iph][Ith] * (interpol_nDth - interpol_oDth) +
                Fmat[m][Iph][Iph] * (newD_ph[i][j][k] - oldD_ph[i][j][k]);

          if(maxE < std::abs(E_phi[New][i][j][k])){
          flag = 3;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_phi[New][i][j][k]);
      }

      }
    }
  }
  
  if(maxE > 1.0e12) exit(0);
  // output "maxE" //
  switch(flag){
    case 1:
    std::cout << "max  E_r[" << I << "][" << J << "][" << K <<"] = " << maxE << std::endl;
    break;

    case 2:
    std::cout << "max  E_th[" << I << "][" << J << "][" << K << "] = " << maxE << std::endl;
    break;

    case 3:
    std::cout << "max  E_phi[" << I << "][" << J << "][" << K << "] = " << maxE << std::endl;
    break;

  }

}
