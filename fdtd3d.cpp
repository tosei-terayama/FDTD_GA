#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <chrono>
#include <Eigen/Core>
#include "fdtd3d.h"
#include "pml.h"
#include "geocoordinate.h"
#include "perturbation.h"
#include "date.h"
//#include <mpi.h>

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

int main(void)
{
  int time_step = 1700;
  double t;
  double J;
  double total_time;
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

  PML_field_initialize(
    Dr_theta1, Dr_theta2, Dr_phi,
    Dtheta_phi, Dtheta_r,
    Dphi_r, Dphi_theta,
    Hr_theta1, Hr_theta2, Hr_phi,
    Htheta_phi, Htheta_r,
    Hphi_r, Hphi_theta
  );

  pml* idx_Dr = new pml[4];
  pml* idx_Dth = new pml[4];
  pml* idx_Dphi = new pml[4];
  pml* idx_Hr = new pml[4];
  pml* idx_Hth = new pml[4];
  pml* idx_Hphi = new pml[4];

  PML_idx_initialize(
    idx_Dr, idx_Dth, idx_Dphi,
    idx_Hr, idx_Hth, idx_Hphi
  );
  
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

  // Geo class //
  geocoordinate lla_info;
  lla_info.set_point(32.0, 135.0, 60.0);

  // Date class (UT)//
  date ymd;
  ymd.set_ymd(2016, 3 ,1);
  ymd.set_h(9.0);

  //Ne, nyu//
  double *Nh = new double[ion_L+1];
  double *ny = new double[ion_L+1];
  double *Re = new double[ion_L+1];

  //iri_profile(ymd, lla_info, Nh, Re);
  Ne_allocate(Nh, Re);
  ny_allocate(ymd, lla_info, ny, Re);

  double *****Cmat = memory_allocate5d(ion_L+1, Ntheta + 1, Nphi + 1, 3, 3, 0.0);
  double *****Fmat = memory_allocate5d(ion_L+1, Ntheta + 1, Nphi + 1, 3, 3, 0.0);
  
  double*** noise_Nh = memory_allocate3d(ion_L + 1, Ntheta + 1, Nphi + 1, 0.0);
  
  perturbation P_info;

  // Set Perturbation Information //
  /*P_info.set_center(74, Ntheta/2, Nphi/2);
  P_info.set_alpha(10.0);
  P_info.set_sigma(2.0e3, 50.0e3);*/

  P_info.set_center(77, 12, 529);
  P_info.set_alpha(5.0);
  P_info.set_sigma(1142.86, 75.0e3);

  set_perturbation(P_info, noise_Nh, Nh);
  set_matrix(zj, Cmat, Fmat, noise_Nh, ny);

  //calculate surface impedance//
  std::complex <double> Z(0.0, 0.0);
  double Z_real, Z_imag;

  Z = surface_impe(zj);

  //get realpart imaginaly part//
  Z_real = Z.real();
  Z_imag = Z.imag()/omega;

  std::ofstream ofs_1;
  ofs_1.open("./dat_file/E0.dat");
  std::ofstream ofs_receive;
  ofs_receive.open("./dat_file/receive.dat");
  std::ofstream ofs_serve;
  ofs_serve.open("./dat_file/serve.dat");
  std::ofstream ofs_j;
  ofs_j.open("./dat_file/J_value.dat");
  std::ofstream ofs_Nphi;
  ofs_Nphi.open("./dat_file/obs_Nphi.dat");
  std::ofstream ofs_NphidB;
  ofs_NphidB.open("./dat_file/target.dat");
  std::ofstream ofs_servedNphi;
  ofs_servedNphi.open("./dat_file/obs_3dNphi.dat");
  std::ofstream ofs_servedNphidB;
  ofs_servedNphidB.open("./dat_file/obs_3dNphidB.dat");

  // output analyze model //
  output_model();
  output_profile(P_info, Nh, noise_Nh);

  t = Dt*0.0;

  for(int k = 0; k < Nphi; k++){
    double Phi = R0*ph(k)/1000.0;
    for(int i = 0; i < Nr; i++){
      double R = i*delta_r/1000.0;
      ofs_1 << Phi << " " << R << " " << Etheta[0][i][j_s][k] << std::endl;
    }
    ofs_1 << std::endl;
  }

  ofs_1.close();
  ofs_receive << 0 << " " << Etheta[0][i_r][j_r][k_r] << std::endl;
  ofs_serve << 0 << " " << Etheta[0][i_s][j_s][k_s] << std::endl;

  std::cout << "R : " << dist(Nr) << " θ : " << R0*delta_theta*Ntheta << " φ : " << R0*ph(Nphi) << std::endl;
  std::cout << "time_step : " << time_step << " Dt : " << Dt << std::endl << std::endl;
  std::cout << "Perturbation r0 : " << P_info.r0() << " th0 : " << P_info.th0() << " phi0 : " << P_info.phi0() << std::endl;
  std::cout << "enhance : " << P_info.alpha() << std::endl;
  std::cout << "_______________________________________" << std::endl;

  ////主経路電波強度観測/////
  int Num_obs = (Nphi - 2*L) - k_s;
  double *Magnitude = new double[Num_obs + 1];

  //fourie//
  std::complex <double>* E_famp = new std::complex <double> [Num_obs + 1];
  std::complex <double>** E_famp3d = memory_allocate2cd(Ntheta + 1, Nphi + 1, std::complex <double> (0.0, 0.0));

  geocoordinate *obs_p = new geocoordinate[Num_obs + 1];
  geocoordinate **obs_p3d = new geocoordinate*[Ntheta + 1];
  for(int j = 0; j <= Ntheta; j++){
    obs_p3d[j] = new geocoordinate[Num_obs + 1];
  }

  obs_ini(obs_p, obs_p3d, Num_obs);  // Initialize observation point //

  for(int k = 0; k < Num_obs; k++){
    E_famp[k] += Er[0][obs_p[k].i()][obs_p[k].j()][obs_p[k].k()]*std::exp(-zj*omega*t)*Dt;
  }

  for(int j = L; j <= Ntheta - L; j++){
    for(int k = 0; k < Num_obs; k++){
      E_famp3d[j][k] += Er[0][obs_p3d[j][k].i()][obs_p3d[j][k].j()][obs_p3d[j][k].k()]*std::exp(-zj*omega*t)*Dt;
    }
  }
  
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
      *std::exp(-std::pow(t - t0, 2.0)/2.0/std::pow(sigma_t, 2.0));

    std::cout << " J = " << J << std::endl;
    
    ofs_j << t << " " << J << std::endl;

    Etheta[OLD][i_s][j_s][k_s] = Etheta[OLD][i_s][j_s][k_s] + J;
    
    std::cout << "E[1][50][20] = " << Er[NEW][1][50][20] << std::endl;

    /////   D, E update   /////
    //outside PML//
    D_update(
      Dr, Dtheta, Dphi, Hr, Htheta, Hphi, NEW, OLD);
    
    //inside PML//
    D_update_pml(
      Dr[NEW], Dtheta[NEW], Dphi[NEW], Hr, Htheta, Hphi, 
      Dr_theta1, Dr_theta2, Dr_phi, Dtheta_phi, Dtheta_r, Dphi_r, Dphi_theta, 
      sigma_theta, sigma_phi, idx_Dr, idx_Dth, idx_Dphi);
    
    //update E using D//
    E_update( 
      Er, Etheta, Ephi, Dr, Dtheta, Dphi, NEW, OLD,
      Cmat, Fmat);

    /////   H update   /////
    //outside PML//
    H_update(
      Er[NEW], Etheta[NEW], Ephi[NEW], Hr, Htheta, Hphi);
    
    //surface Ground//
    surface_H_update(
      Er[NEW][0], Etheta[NEW][1], Ephi[NEW][1], 
      Htheta[0], Hphi[0], Z_real, Z_imag);
    
    //inside PML//
    H_update_pml(
      Er[NEW], Etheta[NEW], Ephi[NEW], Hr, Htheta, Hphi, 
      Hr_theta1, Hr_theta2, Hr_phi, Htheta_phi, Htheta_r, Hphi_r, Hphi_theta, 
      sigma_theta_h, sigma_phi_h, idx_Hr, idx_Hth, idx_Hphi);

    std::string fn = "./dat_file/E" + std::to_string(n) + ".dat";
    ofs_1.open(fn);
    std::ofstream ofs_1(fn.c_str());

    for(int k = 0; k < Nphi; k++){
      double Phi = R0*ph(k)/1000.0;
      for(int i = 0; i < Nr; i++){
        double R = i*delta_r/1000.0;
        ofs_1 << Phi << " " << R << " " << Etheta[NEW][i][j_s][k] << std::endl;
      }
      ofs_1 << std::endl;
    }
    
    ofs_1.close();

    ofs_receive << t << " " << Etheta[NEW][i_r][j_r][k_r] << std::endl;
    ofs_serve << t << " " << Etheta[NEW][i_s][j_s][k_s] << std::endl;

    for(int k = 0; k < Num_obs; k++){
      E_famp[k] += Er[NEW][obs_p[k].i()][obs_p[k].j()][obs_p[k].k()]*std::exp(-zj*omega*t)*Dt;
    }

    for(int j = L; j <= Ntheta - L; j++){
      for(int k = 0; k < Num_obs; k++){
        E_famp3d[j][k] += Er[NEW][obs_p3d[j][k].i()][obs_p3d[j][k].j()][obs_p3d[j][k].k()]*std::exp(-zj*omega*t)*Dt;
      }
    }
    
    std::cout << n << " / " << time_step << std::endl << std::endl;
    
  }
  
  std::chrono::system_clock::time_point end
    = std::chrono::system_clock::now();
  ///////計測終了///////
  
  total_time = std::chrono::duration_cast <std::chrono::milliseconds>
    (end - start).count();
  
  std::cout << "elapsed_time = " << total_time*1.0e-3 << " [sec]"<< std::endl;

  for(int k = 0; k < Num_obs; k++){
    Magnitude[k] = 20.0*std::log10(std::abs(E_famp[k]/E_famp[0]));
    ofs_Nphi << k << " " << std::log10(std::abs(E_famp[k])) << std::endl;
    ofs_NphidB << k << " " << Magnitude[k] << std::endl;
  }

  for(int k = 0; k < Num_obs; k++){
    for(int j = L; j <= Ntheta - L; j++){
      ofs_servedNphi << k << " " << j - L << " " << std::abs(E_famp3d[j][k]) << std::endl;
      ofs_servedNphidB << k << " " << j - L << " " << 20.0*std::log10(std::abs(E_famp3d[j][k]/E_famp3d[j_s][0])) << std::endl;
    }
    ofs_servedNphi << std::endl;
    ofs_servedNphidB << std::endl;
  }

  ofs_1.close();
  ofs_receive.close();
  ofs_serve.close();
  ofs_j.close();
  ofs_Nphi.close();
  ofs_NphidB.close();
  ofs_servedNphi.close();
  ofs_servedNphidB.close();

  delete_5d(Cmat, ion_L + 1, Ntheta + 1, Nphi + 1, 3);
  delete_5d(Fmat, ion_L + 1, Ntheta + 1, Nphi + 1, 3);
  delete_4d(Er, 2, Nr, Ntheta + 1);
  delete_4d(Etheta, 2, Nr + 1, Ntheta);
  delete_4d(Ephi, 2, Nr + 1, Ntheta + 1);
  delete_4d(Dr, 2, Nr, Ntheta + 1);
  delete_4d(Dtheta, 2, Nr + 1, Ntheta);
  delete_4d(Dphi, 2, Nr + 1, Ntheta + 1);
  delete_3d(Hr, Nr + 1, Ntheta);
  delete_3d(Htheta, Nr, Ntheta + 1);
  delete_3d(Hphi, Nr, Ntheta);
  delete_3d(noise_Nh, ion_L + 1, Ntheta + 1);

  delete_PML(Dr_theta1, Dr_theta2, Dr_phi,
            Dtheta_phi, Dtheta_r,
            Dphi_r, Dphi_theta,
            Hr_theta1, Hr_theta2, Hr_phi,
            Htheta_phi, Htheta_r,
            Hphi_r, Hphi_theta);

  delete[] idx_Dr;
  delete[] idx_Dth;
  delete[] idx_Dphi;
  delete[] idx_Hr;
  delete[] idx_Hth;
  delete[] idx_Hphi;
  delete[] sigma_theta;
  delete[] sigma_phi;
  delete[] sigma_theta_h;
  delete[] sigma_phi_h;
  delete[] geo_B;
  delete[] sph_B;
  delete[] Nh;
  delete[] ny;
  delete[] Re;

  delete [] E_famp;
  delete [] Magnitude;
  delete [] E_famp;
  delete [] E_famp3d;
  
  return 0;
  
}
