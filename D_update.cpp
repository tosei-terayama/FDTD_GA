#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void D_update(double**** D_r, double**** D_theta, double**** D_phi, 
              double*** H_r, double*** H_theta, double*** H_phi, 
              int NEW, int OLD)
{
  double val_1, val_2;
  double ri_1, ri_2, ri_3;
  double sin_th1, sin_th2, sin_th3;
  
  //D update (outside PML)//
  for(int i = 0; i < Nr; i++){
    ri_2 = dist(i + 0.5);
    for(int j = L + 1; j < Ntheta - L; j++){
      sin_th1 = std::sin(th(j));
      sin_th2 = std::sin(th(j + 0.5));
      sin_th3 = std::sin(th(j - 0.5));
      for(int k = L + 1; k < Nphi - L; k++){
        val_1 = Dt/ri_2/sin_th1/delta_theta;
        val_2 = Dt/ri_2/sin_th1/delta_phi;

        D_r[NEW][i][j][k] = D_r[OLD][i][j][k] + val_1*(sin_th2*H_phi[i][j][k] - sin_th3*H_phi[i][j-1][k])
          - val_2*(H_theta[i][j][k] - H_theta[i][j][k-1]);
      }
    }
  }
  
  for(int i = 1; i < Nr - 1; i++){
    ri_1 = dist(i);
    ri_2 = dist(i + 0.5);
    ri_3 = dist(i - 0.5);
    for(int j = L; j < Ntheta - L; j++){
      sin_th2 = std::sin(th(j + 0.5));
      for(int k = L + 1; k < Nphi - L; k++){
        val_1 = Dt/ri_1/sin_th2/delta_phi;
        val_2 = Dt/ri_1/delta_r;

        D_theta[NEW][i][j][k] = D_theta[OLD][i][j][k] + val_1*(H_r[i][j][k] - H_r[i][j][k-1])
          - val_2*(ri_2*H_phi[i][j][k] - ri_3*H_phi[i-1][j][k]);
      }
    }
  }
  
  for(int i = 1; i < Nr - 1; i++){
    ri_1 = dist(i);
    ri_2 = dist(i + 0.5);
    ri_3 = dist(i - 0.5);
    for(int j = L + 1; j < Ntheta - L; j++){
      for(int k = L; k < Nphi - L; k++){
        val_1 = Dt/ri_1/delta_r;
        val_2 = Dt/ri_1/delta_theta;

        D_phi[NEW][i][j][k] = D_phi[OLD][i][j][k] + val_1*(ri_2*H_theta[i][j][k] - ri_3*H_theta[i-1][j][k])
          - val_2*(H_r[i][j][k] - H_r[i][j-1][k]); 
      }
    }
  }
  
}
