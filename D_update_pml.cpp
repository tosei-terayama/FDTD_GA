#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"
#include "pml.h"

void D_update_pml(double*** newD_r, double*** newD_theta, double*** newD_phi,
                  double*** H_r, double*** H_theta, double*** H_phi, 
                  double**** Dr_theta1, double**** Dr_theta2, double**** Dr_phi, 
                  double**** Dtheta_phi, double**** Dtheta_r, 
                  double**** Dphi_r, double**** Dphi_theta,
                  double* sigma_theta, double* sigma_phi,
                  pml* idx_Dr, pml* idx_Dth, pml* idx_Dphi)
{
  double ri_1, ri_2, ri_3;
  double thetaj_1, thetaj_2;
  
  //Update Dr using Dr_theta1, Dr_theta2, Dr_phi//
  for(int area = 0; area < 4; area++){
    for(int i = 0; i < Nr; i++){
      ri_2 = dist(i + 0.5);
      for(int j = idx_Dr[area].j1(); j <= idx_Dr[area].j2(); j++){
        int j_area = j - idx_Dr[area].j1();
        thetaj_1 = th(j);
        for(int k = idx_Dr[area].k1(); k <= idx_Dr[area].k2(); k++){
          int k_area = k - idx_Dr[area].k1();
          Dr_theta1[area][i][j_area][k_area] = C_1(sigma_theta[j])*Dr_theta1[area][i][j_area][k_area]
            + C_2(ri_2, sigma_theta[j])*(H_phi[i][j][k] - H_phi[i][j-1][k]);

          Dr_theta2[area][i][j_area][k_area] = Dr_theta2[area][i][j_area][k_area]
            + C_3(ri_2, thetaj_1)*(H_phi[i][j][k] + H_phi[i][j-1][k]);

          Dr_phi[area][i][j_area][k_area] = C_1(sigma_phi[k])*Dr_phi[area][i][j_area][k_area]
            - C_4(ri_2, thetaj_1, sigma_phi[k])*(H_theta[i][j][k] - H_theta[i][j][k-1]);

          newD_r[i][j][k] = Dr_theta1[area][i][j_area][k_area]
            + Dr_theta2[area][i][j_area][k_area] + Dr_phi[area][i][j_area][k_area];
        }
      }
    }
  }
  
  //Update Dtheta using Dtheta_phi, Dtheta_r//
  for(int area = 0; area < 4; area++){
    for(int i = 1; i < Nr; i++){
      ri_1 = dist(i - 0.5);
      ri_2 = dist(i);
      ri_3 = dist(i + 0.5);
      for(int j = idx_Dth[area].j1(); j <= idx_Dth[area].j2(); j++){
        int j_area = j - idx_Dth[area].j1();
        thetaj_2 = th(j + 0.5);
        for(int k = idx_Dth[area].k1(); k <= idx_Dth[area].k2(); k++){
          int k_area = k - idx_Dth[area].k1();
          Dtheta_phi[area][i][j_area][k_area] = C_1(sigma_phi[k])*Dtheta_phi[area][i][j_area][k_area]
          + C_4(ri_2, thetaj_2, sigma_phi[k])*(H_r[i][j][k] - H_r[i][j][k-1]);

          Dtheta_r[area][i][j_area][k_area] = Dtheta_r[area][i][j_area][k_area]
          - C_5(ri_2)*(ri_3*H_phi[i][j][k] - ri_1*H_phi[i-1][j][k]);

          newD_theta[i][j][k] = Dtheta_phi[area][i][j_area][k_area] + Dtheta_r[area][i][j_area][k_area];
        }
      }
    }
  }
  
  //Update Dphi using Dphi_r, Dphi_theta//
  for(int area = 0; area < 4; area++){
    for(int i = 1; i < Nr; i++){
      ri_1 = dist(i - 0.5);
      ri_2 = dist(i);
      ri_3 = dist(i + 0.5);
      for(int j = idx_Dphi[area].j1(); j <= idx_Dphi[area].j2(); j++){
        int j_area = j - idx_Dphi[area].j1();
        for(int k = idx_Dphi[area].k1(); k <= idx_Dphi[area].k2(); k++){
          int k_area = k - idx_Dphi[area].k1();
          Dphi_r[area][i][j_area][k_area] = Dphi_r[area][i][j_area][k_area]
            + C_5(ri_2)*(ri_3*H_theta[i][j][k] - ri_1*H_theta[i-1][j][k]);
          
          Dphi_theta[area][i][j_area][k_area] = C_1(sigma_theta[j])*Dphi_theta[area][i][j_area][k_area]
            - C_6(ri_2, sigma_theta[j])*(H_r[i][j][k] - H_r[i][j-1][k]);
          
          newD_phi[i][j][k] = Dphi_r[area][i][j_area][k_area] + Dphi_theta[area][i][j_area][k_area];
        }
      }
    }
  }

}