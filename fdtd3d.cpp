#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <chrono>
#include <Eigen/Core>

#include "fdtd3d.h"
//#include <mpi.h>

//The number of R, Theta, Phi element//
const int Nr{100};
const int Ntheta{100};
const int Nphi{600};

constexpr double R_r{100.0e3};

//Minute R, Theta, Phi, Time//
const double delta_r = R_r/double(Nr);
const double delta_theta = 1.0e3/double(R0);
const double delta_phi = 1.0e3/double(R0);
const double Dt = 0.99/C0/std::sqrt(1.0/delta_r/delta_r
 + 1.0/R0/R0/delta_theta/delta_theta
 + 1.0/R0/R0/std::sin(THETA0)/std::sin(THETA0)/delta_phi/delta_phi);
const double inv_Dt = 1.0/Dt;
const double sigma_t = 7.0*Dt;
const double t0 = 6.0*sigma_t;

//center_point//
const int i_0 = Nr/2;
const int j_0 = Ntheta/2;
const int k_0 = Nphi/2;

//Source_point//
const int i_s{1};
const int j_s{50};
const int k_s{50};

//Receive_point//
const int i_r{1};
const int j_r{50};
const int k_r{Nphi - 50};

//PML info//
const int L{10};
const double M{3.5};
const double R{1.0e-6};

const double sigma_th_max = -(M + 1.0)*C0*std::log(R)/2.0/double(L)/delta_theta/R0;
const double sigma_phi_max = -(M + 1.0)*C0*std::log(R)/2.0/double(L)/delta_phi/R0;

//Ionosphere info//
constexpr double Alt_lower_ionosphere{60.0e3};
const int ion_L = int((R_r - Alt_lower_ionosphere)/delta_r);
const double freq{22.2e3};
const double omega = 2.0*M_PI*freq;

//Geomagnetic info//
const double B_abs{4.6468e-5};
const double Dec{-7.0*M_PI/180.0};
const double Inc{49.0*M_PI/180.0};
const double Azim{61.0*M_PI/180.0};

int main(int argc, char** argv)
{
  int flag(0);
  int time_step = 1200;
  double t;
  double J;
  double time_1, time_2, total_time;
  int NEW;
  int OLD;
  std::complex <double> zj(0.0, 1.0);

  double*** Hr, ***Htheta, ***Hphi;
  Hr = memory_allocate3d(Nr + 1, Ntheta, Nphi, 0.0);
  Htheta = memory_allocate3d(Nr, Ntheta + 1, Nphi, 0.0);
  Hphi = memory_allocate3d(Nr, Ntheta, Nphi + 1, 0.0);
  
  double**** Er, **** Etheta, **** Ephi;
  Er = memory_allocate4d(2, Nr, Ntheta + 1, Nphi + 1, 0.0);
  Etheta = memory_allocate4d(2, Nr + 1, Ntheta, Nphi + 1, 0.0);
  Ephi = memory_allocate4d(2, Nr + 1, Ntheta + 1, Nphi, 0.0);

  double**** Dr, **** Dtheta, **** Dphi;
  Dr = memory_allocate4d(2, Nr, Ntheta + 1, Nphi + 1, 0.0);
  Dtheta = memory_allocate4d(2, Nr + 1, Ntheta, Nphi + 1, 0.0);
  Dphi = memory_allocate4d(2, Nr + 1, Ntheta + 1, Nphi, 0.0);

  double**** Dr_theta1, **** Dr_theta2, **** Dr_phi;
  double**** Dtheta_phi, **** Dtheta_r;
  double**** Dphi_r, **** Dphi_theta;

  Dr_theta1 = new double***[4];
  Dr_theta2 = new double***[4];
  Dr_phi = new double***[4];
  Dtheta_phi = new double***[4];
  Dtheta_r = new double***[4];
  Dphi_r = new double***[4];
  Dphi_theta = new double***[4];
  
  double**** Hr_theta1, ****Hr_theta2, ****Hr_phi;
  double**** Htheta_phi, ****Htheta_r;
  double**** Hphi_r, ****Hphi_theta;

  Hr_theta1 = new double***[4];
  Hr_theta2 = new double***[4];
  Hr_phi = new double***[4];
  Htheta_phi = new double***[4];
  Htheta_r = new double***[4];
  Hphi_r = new double***[4];
  Hphi_theta = new double***[4];

  //PML region (Theta direction)//
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
  
  double *sigma_theta, *sigma_phi, *sigma_theta_h, *sigma_phi_h;
  sigma_theta = new double[Ntheta + 1];
  sigma_phi = new double[Nphi + 1];
  sigma_theta_h = new double[Ntheta + 1];
  sigma_phi_h = new double[Nphi + 1];

  sigma_calc(sigma_theta, sigma_phi, sigma_theta_h, sigma_phi_h);

  //Geomagnetic field//
  double *geo_B = new double[3];
  double *sph_B = new double[3];
  double B_th(0.0), B_phi(0.0);

  geo_mag(geo_B, sph_B);

  B_th = std::acos(-sph_B[1]/B_abs);
  B_phi = std::atan2(sph_B[2], sph_B[0]);

  std::cout << "B_theta = " << B_th << "\tB_phi = " << B_phi << std::endl;

  //Ne, nyu//
  double *Nh = new double[ion_L];
  double *ny = new double[ion_L];

  double *Nh_h = new double[ion_L];
  double *ny_h = new double[ion_L];

  double *Re = new double[ion_L];
  double *Re_h = new double[ion_L];

  //Rotation matrix//
  double **R2_1 = memory_allocate2d(3, 3, 0.0);
  double **invR1_2 = memory_allocate2d(3, 3, 0.0);

  //Corrdinate_transform_matrix//
  double ****Cmat_r = memory_allocate4d(ion_L, Ntheta, Nphi, 9, 0.0);
  double ****Fmat_r = memory_allocate4d(ion_L, Ntheta, Nphi, 9, 0.0);
  double ****Cmat_th = memory_allocate4d(ion_L, Ntheta, Nphi, 9, 0.0);
  double ****Fmat_th = memory_allocate4d(ion_L, Ntheta, Nphi, 9, 0.0);
  double ****Cmat_phi = memory_allocate4d(ion_L, Ntheta, Nphi, 9, 0.0);
  double ****Fmat_phi = memory_allocate4d(ion_L, Ntheta, Nphi, 9, 0.0);

  //sigma_matrix(complex_number, real_part)//
  double ***sigma_real
    = memory_allocate3d(ion_L, 3, 3, 0.0);

  double ***sigma_real_r
    = memory_allocate3d(ion_L, 3, 3, 0.0);

  double ***sigma_cartesian
    = memory_allocate3d(ion_L, 3, 3, 0.0);

  double ***sigma_cartesian_r
    = memory_allocate3d(ion_L, 3, 3, 0.0);

  double b_th(0.0), b_phi(0.0);

  Ne_allocate(Nh, Nh_h, Re, Re_h);

  ny_allocate(ny, ny_h, Re, Re_h);

  make_rot_mat(R2_1, invR1_2, B_th, B_phi);
  
  sig_real_calc(Nh, ny, Nh_h, ny_h, sigma_real, sigma_real_r);

  sig_car_calc(sigma_cartesian, sigma_real, R2_1, invR1_2);
  sig_car_calc(sigma_cartesian_r, sigma_real_r, R2_1, invR1_2);

  coordinate_trans(Cmat_r, Fmat_r, Cmat_th, Fmat_th, Cmat_phi, Fmat_phi,
                  sigma_cartesian, sigma_cartesian_r);

  //calculate surface impedance//
  std::complex <double> Z(0.0, 0.0);
  double Z_real, Z_imag;

  Z = surface_impe(zj);

  //get realpart imaginaly part//
  Z_real = Z.real();
  Z_imag = Z.imag()/omega;

  std::ofstream ofs_1;
  ofs_1.open("./dat_file/E0.dat");

  std::ofstream ofs_2;
  ofs_2.open("./dat_file/receive.dat");

  std::ofstream ofs_3;
  ofs_3.open("./dat_file/serve.dat");

  std::ofstream ofs_4;
  ofs_4.open("./dat_file/fourie.dat");

  std::ofstream ofs_5;
  ofs_5.open("./dat_file/gain_graph.dat");

  std::ofstream ofs_j;
  ofs_j.open("./dat_file/J_value.dat");

  std::complex <double> *E_famp;
  E_famp = new std::complex<double> [Nphi + 1];
  for(int i = 0; i <= Nphi; i++) E_famp[i] = (0.0, 0.0);

  t = Dt*0;

  for(int k = L; k < Nphi - L; k++){
    double Phi = R0*k*delta_phi/1000.0;
    for(int i = 0; i < Nr; i++){
      double R = i*delta_r/1000.0;
      ofs_1 << Phi << " " << R << " " << Etheta[NEW][i][j_s][k] << std::endl;
    }
    ofs_1 << std::endl;
  }

  ofs_1.close();

  //fourie//
  for(int k = 0; k < Nphi; k++){
    E_famp[k] += Etheta[NEW][1][Ntheta/2][k]*std::exp(-zj*omega*0.0)*Dt;
  }

  ofs_2 << 0 << " " << Etheta[NEW][i_r][j_r][k_r] << std::endl;
  ofs_3 << 0 << " " << Etheta[NEW][i_s][j_s][k_s] << std::endl;

  std::cout << "R : " << dist(Nr) << " θ : " << R0*delta_theta*Ntheta << " φ : " << R0*delta_phi*Nphi << std::endl;
  std::cout << "time_step : " << time_step << " Dt : " << Dt << std::endl;
  std::cout << "_______________________________________" << std::endl;
  
  ////////計測開始////////
  std::chrono::system_clock::time_point start
    = std::chrono::system_clock::now();

  //FDTD_update//
  for(int n = 1; n < time_step + 1; n++){
    
    NEW = n%2;
    OLD = (n+1)%2;
    
    //t = (double(n) - 0.5)*Dt;
    t = n*Dt;
    
    //Forced current//
    J = -((t - t0)/sigma_t/sigma_t/delta_r/(dist(i_s + 0.5)*delta_theta)/(dist(i_s + 0.5)*delta_phi))
      *std::exp(-(t - t0)*(t - t0)/2.0/sigma_t/sigma_t);
    //if(t < t0) J = std::exp(-(t - t0)*(t - t0)/2.0/sigma_t/sigma_t)*std::sin(2.0*M_PI*freq*t);
    //else J = std::sin(2.0*M_PI*freq*t);
    std::cout << " J = " << J << std::endl;

    ofs_j << t << " " << J << std::endl;

    Etheta[OLD][i_s][j_s][k_s] = Etheta[OLD][i_s][j_s][k_s] + J;
    
    // std::cout << "Etheta[" << i_s << "][" << j_s << "][" << k_s << "] = " << Etheta[NEW][i_s][j_s][k_s] << std::endl;
    //std::cout << "Etheta[" << Nr - ion_L << "][" << j_0 << "][" << k_0 << "] = " << Etheta[NEW][Nr - ion_L][j_0][k_0] << std::endl;
    std::cout << "E[1][" << Ntheta/2 << "][" << 30 << "] = " << Etheta[NEW][1][Ntheta/2][30] << std::endl;

    /////   D, E update   /////
    //outside PML//
    D_update(Dr, Dtheta, Dphi, Hr, Htheta, Hphi, NEW, OLD);
    
    //inside PML//
    D_update_pml(Dr[NEW], Dtheta[NEW], Dphi[NEW], Hr, Htheta, Hphi, 
    Dr_theta1, Dr_theta2, Dr_phi, Dtheta_phi, Dtheta_r, Dphi_r, Dphi_theta, 
    sigma_theta, sigma_phi);
    
    //update E using D//
    E_update( Er, Etheta, Ephi, Dr, Dtheta, Dphi, NEW, OLD,
             sigma_cartesian, sigma_cartesian_r, Cmat_r, Fmat_r,
             Cmat_th, Fmat_th, Cmat_phi, Fmat_phi);

    /////   H update   /////
    //outside PML//
    H_update(Er[NEW], Etheta[NEW], Ephi[NEW], Hr, Htheta, Hphi);
    
    //inside PML//
    H_update_pml(Er[NEW], Etheta[NEW], Ephi[NEW], Hr, Htheta, Hphi, 
    Hr_theta1, Hr_theta2, Hr_phi, Htheta_phi, Htheta_r, Hphi_r, Hphi_theta, 
    sigma_theta_h, sigma_phi_h);

    //surface Ground//
    surface_H_update(Er[NEW][0], Etheta[NEW][1], Ephi[NEW][1], Htheta[0], Hphi[0],
                    Z_real, Z_imag);
    
    std::string fn = "./dat_file/E" + std::to_string(n) + ".dat";
    ofs_1.open(fn);
    std::ofstream ofs_1(fn.c_str());

    for(int k = L; k < Nphi - L; k++){
      double Phi = R0*k*delta_phi/1000.0;
      for(int i = 0; i < Nr; i++){
        double R = i*delta_r/1000.0;
        ofs_1 << Phi << " " << R << " " << Etheta[NEW][i][j_s][k] << std::endl;
      }
      ofs_1 << std::endl;
    }
    
    ofs_1.close();

    ofs_2 << t << " " << Etheta[NEW][i_r][j_r][k_r] << std::endl;
    ofs_3 << t << " " << Etheta[NEW][i_s][j_s][k_s] << std::endl;

    for(int k = 0; k < Nphi; k++){
      E_famp[k] += Etheta[NEW][1][Ntheta/2][k]*std::exp(-zj*omega*t)*Dt;
    }
    
    std::cout << n << " / " << time_step << std::endl << std::endl;
    
  }
  
  std::chrono::system_clock::time_point end
    = std::chrono::system_clock::now();
  ///////計測終了///////

  //MPI::Finalize();
  
  total_time = std::chrono::duration_cast <std::chrono::milliseconds>
    (end - start).count();
  
  std::cout << "elapsed_time = " << total_time*1.0e-3 << " [sec]"<< std::endl;

  for(int k = k_s; k < k_r; k++){
    ofs_4 << R0*k*delta_phi/1000.0 << " " << std::abs(E_famp[k]) << std::endl;
  }

  for(int k = k_s; k < k_r + 1; k++){
    ofs_5 << k << " " << 20.0*std::log10(std::abs(E_famp[k]/E_famp[k_s])) << std::endl;
  }

  ofs_1.close();
  ofs_2.close();
  ofs_3.close();
  ofs_4.close();
  ofs_5.close();
  ofs_j.close();

  delete [] Er;
  delete [] Etheta;
  delete [] Ephi;
  delete [] Dr;
  delete [] Dtheta;
  delete [] Dphi;
  delete [] Dr_theta1;
  delete [] Dr_theta2;
  delete [] Dr_phi;
  delete [] Dtheta_r;
  delete [] Dtheta_phi;
  delete [] Dphi_r;
  delete [] Dphi_theta;
  delete [] Hr;
  delete [] Htheta;
  delete [] Hphi;
  delete [] Hr_theta1;
  delete [] Hr_theta2;
  delete [] Htheta_r;
  delete [] Htheta_phi;
  delete [] Hphi_r;
  delete [] Hphi_theta;
  
  return 0;
  
}



