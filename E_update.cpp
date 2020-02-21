#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void E_update(double**** E_r, double**** E_theta, double**** E_phi,
              double**** D_r, double**** D_theta, double**** D_phi,
              int NEW, int OLD, double*** sigma_car, double*** sigma_car_r,
              double ****Cmat_r, double ****Fmat_r, double ****Cmat_th, double ****Fmat_th,
              double ****Cmat_phi, double ****Fmat_phi)
{
  double theta(0.0), phi(0.0);
  double maxE(0.0);
  double interpol_Er(0.0), interpol_Eth(0.0), interpol_Eph(0.0);
  double interpol_nDr(0.0), interpol_nDth(0.0), interpol_nDph(0.0),
  interpol_oDr(0.0), interpol_oDth(0.0), interpol_oDph(0.0);
  int flag(0);
  int iono_flag(0);
  int I, J, K;

  for(int i = 0; i < Nr - ion_L; i++){
    for(int j = 1; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
        E_r[NEW][i][j][k] = E_r[OLD][i][j][k] + (D_r[NEW][i][j][k] - D_r[OLD][i][j][k])/EPS0;
        
        if(maxE < std::abs(E_r[NEW][i][j][k])){
          flag = 1;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_r[NEW][i][j][k]);
      }

      }
    }
  }

   iono_flag = 1;
  for(int i = Nr - ion_L; i < Nr; i++){
    for(int j = 1; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
          phi = k*delta_phi;
          
          interpol_Eth = (E_theta[OLD][i][j][k] + E_theta[OLD][i][j - 1][k]
           + E_theta[OLD][i + 1][j][k] + E_theta[OLD][i + 1][j - 1][k])/4.0;
          interpol_nDth = (D_theta[NEW][i][j][k] + D_theta[NEW][i][j - 1][k]
           + D_theta[NEW][i + 1][j][k] + D_theta[NEW][i + 1][j - 1][k])/4.0;
          interpol_oDth = (D_theta[OLD][i][j][k] + D_theta[OLD][i][j - 1][k]
           + D_theta[OLD][i + 1][j][k] + D_theta[OLD][i + 1][j - 1][k])/4.0;

          interpol_Eph = (E_phi[OLD][i][j][k] + E_phi[OLD][i][j][k - 1]
           + E_phi[OLD][i + 1][j][k] + E_phi[OLD][i + 1][j][k - 1])/4.0;
          interpol_nDph = (D_phi[NEW][i][j][k] + D_phi[NEW][i][j][k - 1]
           + D_phi[NEW][i + 1][j][k] + D_phi[NEW][i + 1][j][k - 1])/4.0;
          interpol_oDph = (D_phi[OLD][i][j][k] + D_phi[OLD][i][j][k - 1]
           + D_phi[OLD][i + 1][j][k] + D_phi[OLD][i + 1][j][k - 1])/4.0;
           
         E_r[NEW][i][j][k] = E_update_iono(sigma_car_r[i - Nr + ion_L], E_r[OLD][i][j][k], interpol_Eth, interpol_Eph,
         D_r[NEW][i][j][k], interpol_nDth, interpol_nDph, D_r[OLD][i][j][k], interpol_oDth, interpol_oDph,
         iono_flag, Cmat_r[i - Nr + ion_L][j][k], Fmat_r[i - Nr + ion_L][j][k]);

          if(maxE < std::abs(E_r[NEW][i][j][k])){
          flag = 1;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_r[NEW][i][j][k]);
      }

      }
    }
  }

  for(int i = 1; i < Nr - ion_L; i++){
    for(int j = 0; j < Ntheta; j++){
      for(int k = 1; k < Nphi; k++){
        E_theta[NEW][i][j][k] = E_theta[OLD][i][j][k] + (D_theta[NEW][i][j][k] - D_theta[OLD][i][j][k])/EPS0;

         if(maxE < std::abs(E_theta[NEW][i][j][k])){
          flag = 2;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_theta[NEW][i][j][k]);
      }

      }
    }
  }

   iono_flag = 2;
  for(int i = Nr - ion_L; i < Nr; i++){
    for(int j = 0; j < Ntheta; j++){
      theta = th(j + 0.5);
      for(int k = 1; k < Nphi; k++){
          phi = k*delta_phi;

          interpol_Er = (E_r[OLD][i][j][k] + E_r[OLD][i - 1][j][k]
           + E_r[OLD][i][j + 1][k] + E_r[OLD][i - 1][j + 1][k])/4.0;
          interpol_nDr = (D_r[NEW][i][j][k] + D_r[NEW][i - 1][j][k]
           + D_r[NEW][i][j + 1][k] + D_r[NEW][i - 1][j + 1][k])/4.0;
          interpol_oDr = (D_r[OLD][i][j][k] + D_r[OLD][i - 1][j][k]
           + D_r[OLD][i][j + 1][k] + D_r[OLD][i - 1][j + 1][k])/4.0;

          interpol_Eph = (E_phi[OLD][i][j][k] + E_phi[OLD][i][j][k - 1]
           + E_phi[OLD][i][j + 1][k] + E_phi[OLD][i][j + 1][k - 1])/4.0;
          interpol_nDph = (D_phi[NEW][i][j][k] + D_phi[NEW][i][j][k - 1]
           + D_phi[NEW][i][j + 1][k] + D_phi[NEW][i][j + 1][k - 1])/4.0;
          interpol_oDph = (D_phi[OLD][i][j][k] + D_phi[OLD][i][j][k - 1]
           + D_phi[OLD][i][j + 1][k] + D_phi[OLD][i][j + 1][k - 1])/4.0;

          E_theta[NEW][i][j][k] = E_update_iono(sigma_car[i - Nr + ion_L], interpol_Er, E_theta[OLD][i][j][k], interpol_Eph,
          interpol_nDr, D_theta[NEW][i][j][k], interpol_nDph, interpol_oDr, D_theta[OLD][i][j][k], interpol_oDph,
          iono_flag, Cmat_th[i - Nr + ion_L][j][k], Fmat_th[i - Nr + ion_L][j][k]);

          if(maxE < std::abs(E_theta[NEW][i][j][k])){
          flag = 2;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_theta[NEW][i][j][k]);
      }

      }
    } 
  }

  for(int i = 1; i < Nr - ion_L; i++){
    for(int j = 1; j < Ntheta; j++){
      for(int k = 0; k < Nphi; k++){
        E_phi[NEW][i][j][k] = E_phi[OLD][i][j][k] + (D_phi[NEW][i][j][k] - D_phi[OLD][i][j][k])/EPS0;

        if(maxE < std::abs(E_phi[NEW][i][j][k])){
          flag = 3;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_phi[NEW][i][j][k]);
      }

      }
    }
  }

  iono_flag = 3;
  for(int i = Nr - ion_L; i < Nr; i++){
    for(int j = 1; j < Ntheta; j++){
      theta = th(j);
      for(int k = 0; k < Nphi; k++){
          phi = (k + 0.5)*delta_phi;

          interpol_Er = (E_r[OLD][i][j][k] + E_r[OLD][i - 1][j][k]
           + E_r[OLD][i][j][k + 1] + E_r[OLD][i - 1][j][k + 1])/4.0;
          interpol_nDr = (D_r[NEW][i][j][k] + D_r[NEW][i - 1][j][k]
           + D_r[NEW][i][j][k + 1] + D_r[NEW][i - 1][j][k + 1])/4.0;
          interpol_oDr = (D_r[OLD][i][j][k] + D_r[OLD][i - 1][j][k]
           + D_r[OLD][i][j][k + 1] + D_r[OLD][i - 1][j][k + 1])/4.0;

          interpol_Eth = (E_theta[OLD][i][j][k] + E_theta[OLD][i][j - 1][k]
           + E_theta[OLD][i][j][k + 1] + E_theta[OLD][i][j - 1][k + 1])/4.0;
          interpol_nDth = (D_theta[NEW][i][j][k] + D_theta[NEW][i][j - 1][k]
           + D_theta[NEW][i][j][k + 1] + D_theta[NEW][i][j - 1][k + 1])/4.0;
          interpol_oDth = (D_theta[OLD][i][j][k] + D_theta[OLD][i][j - 1][k]
           + D_theta[OLD][i][j][k + 1] + D_theta[OLD][i][j - 1][k + 1])/4.0;

          E_phi[NEW][i][j][k] = E_update_iono(sigma_car[i - Nr + ion_L], interpol_Er, interpol_Eth, E_phi[OLD][i][j][k],
          interpol_nDr, interpol_nDth, D_phi[NEW][i][j][k], interpol_oDr, interpol_oDth, D_phi[OLD][i][j][k],
          iono_flag, Cmat_phi[i - Nr + ion_L][j][k], Fmat_phi[i - Nr + ion_L][j][k]);

          if(maxE < std::abs(E_phi[NEW][i][j][k])){
          flag = 3;
          I = i;
          J = j;
          K = k;
          maxE = std::abs(E_phi[NEW][i][j][k]);
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
