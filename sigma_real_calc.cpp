#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

#include "nrlmsise-00.h"
#include "fdtd3d.h"

void sig_real_calc(double *Nh, double *ny, double *Nh_h, double *ny_h, 
                    double ***sigma_re, double ***sigma_re_r)
{

  double height(0.0);

  std::complex <double> ***sigma_com
    = memory_allocate3cd(ion_L, 3, 3, std::complex <double> (0.0, 0.0));
  
  std::complex <double> ***sigma_com_r
    = memory_allocate3cd(ion_L, 3, 3, std::complex <double> (0.0, 0.0));

  std::complex <double> re(1.0, 0.0);
  std::complex <double> zj(0.0, 1.0);

  std::complex <double> omega_cyc(0.0, 0.0);

  std::complex <double> omega_plasma(0.0, 0.0);
  std::complex <double> omega_dash(0.0, 0.0);
  std::complex <double> omega_diff(0.0, 0.0);
  std::complex <double> coeff(0.0, 0.0);

  std::complex <double> omega_plasma_r(0.0, 0.0);
  std::complex <double> omega_dash_r(0.0, 0.0);
  std::complex <double> omega_diff_r(0.0, 0.0);
  std::complex <double> coeff_r(0.0, 0.0);
  
  std::ofstream ofs_1("Nh.dat");
  std::ofstream ofs_2("ny.dat");

  //calclation inverse matrix in LHS yields//
  for(int i = 0; i < ion_L; i++){

    height = Nr - ion_L + i;

    //cyclotron frequency//
    omega_cyc = re*E_Q*B_abs/E_M;

    ofs_1 << Nh[i] << " " << height << std::endl;
    ofs_2 << ny[i] << " " << height << std::endl;

    omega_dash = re*omega - zj*ny[i];
    omega_plasma = re*std::sqrt(Nh[i]*E_Q*E_Q/E_M/EPS0);
    omega_diff = (omega_cyc*omega_cyc) - (omega_dash*omega_dash);

    coeff = zj*EPS0*omega_plasma*omega_plasma;

    sigma_com[i][0][0] = omega_dash/omega_diff;
    sigma_com[i][0][1] = zj*omega_cyc/omega_diff;
    sigma_com[i][1][0] = -zj*omega_cyc/omega_diff;
    sigma_com[i][1][1] = omega_dash/omega_diff;
    sigma_com[i][2][2] = -1.0/omega_dash;

    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
        sigma_com[i][j][k] = coeff*sigma_com[i][j][k];
      }
    }

    omega_dash_r = re*omega - zj*ny_h[i];
    omega_plasma_r = re*std::sqrt(Nh_h[i]*E_Q*E_Q/E_M/EPS0);
    omega_diff_r = (omega_cyc*omega_cyc) - (omega_dash_r*omega_dash_r);

    coeff_r = zj*EPS0*omega_plasma_r*omega_plasma_r;

    sigma_com_r[i][0][0] = omega_dash_r/omega_diff_r;
    sigma_com_r[i][0][1] = zj*omega_cyc/omega_diff_r;
    sigma_com_r[i][1][0] = -zj*omega_cyc/omega_diff_r;
    sigma_com_r[i][1][1] = omega_dash_r/omega_diff_r;
    sigma_com_r[i][2][2] = -1.0/omega_dash_r;

    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
        sigma_com_r[i][j][k] = coeff_r*sigma_com_r[i][j][k];
      }
    }

  }

  //get real part//
  for(int i = 0; i < ion_L; i++){
    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
        sigma_re[i][j][k] = sigma_com[i][j][k].real();
        sigma_re_r[i][j][k] = sigma_com_r[i][j][k].real();
      }
    }
  }

  ofs_1.close();
  ofs_2.close();

}