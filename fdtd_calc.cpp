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
#include "GA_agent.h"
//#include <mpi.h>

const double R_r{ 100.0e3 };

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

void fdtd_calc(perturbation P_info, date ymd, geocoordinate lla_info,
            int Num_obs, geocoordinate* obs_p, double* Magnitude, int myrank)
{
  for(int i = 0; i < Num_obs; i++) Magnitude[i] = 0.0;

  int time_step = 2000;
  double t;
  double J;
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
  //double B_th(0.0), B_phi(0.0);

  geo_mag(geo_B, sph_B);

  //B_th = std::acos(-sph_B[1]/B_abs);
  //B_phi = std::atan2(sph_B[2], sph_B[0]);

  //Ne, nyu//
  double *Nh = new double[ion_L+1];
  double*** noise_Nh = memory_allocate3d(ion_L+1, Ntheta+1, Nphi+1, 0.0);
  double *ny = new double[ion_L+1];
  double *Re = new double[ion_L+1];

  //iri_profile(ymd, lla_info, Nh, Re);
  Ne_allocate(Nh, ny);
  ny_allocate(ymd, lla_info, ny, Re);

  double *****Cmat = memory_allocate5d(ion_L+1, Ntheta, Nphi, 3, 3, 0.0);
  double *****Fmat = memory_allocate5d(ion_L+1, Ntheta, Nphi, 3, 3, 0.0);

  set_perturbation(P_info, noise_Nh, Nh);

  set_matrix(zj, Cmat, Fmat, noise_Nh, ny);

  //calculate surface impedance//
  std::complex <double> Z(0.0, 0.0);
  double Z_real, Z_imag;

  Z = surface_impe(zj);

  //get realpart imaginaly part//
  Z_real = Z.real();
  Z_imag = Z.imag()/omega;
  
  t = Dt*0.0;

  //fourie//
  std::complex <double>* E_famp = new std::complex <double> [Num_obs + 1];
  
  for(int k = 0; k < Num_obs; k++){
    E_famp[k] += Er[0][obs_p[k].i()][obs_p[k].j()][obs_p[k].k()]*std::exp(-zj*omega*t)*Dt;
  }

  int n = 1;
  
  //FDTD_update//
  for(n = 1; n < time_step + 1; n++){

    NEW = n%2;
    OLD = (n + 1)%2;
    
    if(myrank == 1 && n % 100 == 0){
      std::cout << n << " / " << time_step << "  Er(50, 50, 500)" << Er[NEW][50][50][500] << std::endl;
    }
    
    //t = (double(n) - 0.5)*Dt;
    t = n*Dt;
    
    //Forced current//
    J = -((t - t0)/sigma_t/sigma_t/delta_r/(dist(i_s + 0.5)*delta_theta)/(dist(i_s + 0.5)*delta_phi))
      *std::exp(-std::pow(t - t0, 2.0)/2.0/std::pow(sigma_t, 2.0));

    Etheta[OLD][i_s][j_s][k_s] = Etheta[OLD][i_s][j_s][k_s] + J;

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

    for(int k = 0; k < Num_obs; k++){
      E_famp[k] += Er[NEW][obs_p[k].i()][obs_p[k].j()][obs_p[k].k()]*std::exp(-zj*omega*t)*Dt;
    }
    
  }
  
  for(int k = 0; k < Num_obs; k++){
    Magnitude[k] = 20.0*std::log10(std::abs(E_famp[k]/E_famp[0]));
  }

  delete_5d(Cmat, ion_L + 1, Ntheta, Nphi, 3);
  delete_5d(Fmat, ion_L + 1, Ntheta, Nphi, 3);
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

}
