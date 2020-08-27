#ifndef FDTD_H_
#define FDTD_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include "pml.h"
#include "geocoordinate.h"
#include "perturbation.h"
#include "date.h"

//Physical quantity//
#define C0 (3.0e8)
#define MU0 (4.0*M_PI*1.0e-7)
#define EPS0 (1.0/MU0/C0/C0)
#define EPSR (10.0)
#define Z0 (std::sqrt(MU0/EPS0))
#define R0 (6370e3)
#define THETA0 (M_PI*0.5 - std::atan(50e3/R0))
#define E_Q (1.6e-19)
#define E_M (9.11e-31)
#define SIGMA_PEC (1.0e7)
#define SIGMA_SEA (1.0)
#define SIGMA_WET_GROUND (1.0e-2)
#define SIGMA_DRY_GROUND (1.0e-3)
#define SIGMA_VERY_DRY_GROUND (1.0e-4)
#define SIGMA_FRESH_WATER_ICE (1.0e-5)

//Minute R, Theta, Phi, Time//
extern const double delta_r;
extern const double delta_theta;
extern const double delta_phi;
extern const double Dt;
extern const double inv_Dt;
extern const double sigma_t;
extern const double t0;

//The number of R, Theta, Phi element//
const int Nr{ 100 };
const int Ntheta{ 100 };
const int Nphi{ 1000 };

extern const double R_r;

// Center Point //
const int i_0 = Nr/2;
const int j_0 = Ntheta/2;
const int k_0 = Nphi/2;

// Source Point //
const int i_s{ 1 };
const int j_s{ 50 };
const int k_s{ 100 };       // phi = 100 //

// Receive Point //
const int i_r{ 1 };
const int j_r{ 50 };
const int k_r{ Nphi - 50 };       // phi = 950 //

//PML information//
extern const int L;
extern const double M;
extern const double R;
extern const double sigma_th_max;
extern const double sigma_phi_max;

//Ionosphere info//
extern const double Alt_lower_ionosphere;
extern const int ion_L;
extern const double freq;
extern const double omega;

//Geomagnetic info//
extern const double B_abs;
extern const double Dec;
extern const double Inc;
extern const double Azim;

//function in main//
double** memory_allocate2d(int, int, double);
double*** memory_allocate3d(int, int, int, double);
double**** memory_allocate4d(int, int, int, int, double);
double***** memory_allocate5d(int, int, int, int, int, double);
std::complex <double>*** memory_allocate3cd(int, int, int, std::complex<double>);
std::complex <double>** memory_allocate2cd(int, int, std::complex<double>);

void delete_2d(double** , int);
void delete_3d(double***, int, int);
void delete_4d(double****, int, int, int);
void delete_5d(double*****, int, int, int, int);
void delete_2cd(std::complex <double>**, int);
void delete_3cd(std::complex <double>***, int, int);
void delete_PML(double**** Dr_theta1, double**** Dr_theta2, double**** Dr_phi,
                double**** Dtheta_phi, double**** Dtheta_r, double**** Dphi_r, double**** Dphi_theta,
                double**** Hr_theta1, double**** Hr_theta2, double**** Hr_phi,
                double**** Htheta_phi, double**** Htheta_r, double**** Hphi_r, double**** Hphi_theta);

void output_profile(
    perturbation perturbation_information, double* Nh, double*** Nh_with_noise);
void output_model(void);

void fdtd_calc(
    perturbation* Perturbation_info, date ymd, geocoordinate lla_info,
    int Num_obs, geocoordinate* obs_p, double* Magnitude);

void sigma_calc(
    double* sigma_theta, double* sigma_phi, 
    double* sigma_theta_h, double* sigma_phi_h);

void D_update(
    double**** Dr, double**** Dtheta, double**** Dphi,
    double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_n_minus_halfDt,
    int num_of_new, int num_of_old);

void D_update_pml(
    double*** Dr_at_nDt, double*** Dtheta_at_nDt, double*** Dphi_at_nDt,
	double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_n_minus_halfDt,
	double**** Dr_theta1_n_minus_oneDt, double**** Dr_theta2_n_minus_halfDt, double**** Dr_phi_n_minus_oneDt,
	double**** Dtheta_phi_at_n_minus_oneDt, double**** Dtheta_r_at_n_minus_oneDt,
	double**** Dphi_r_at_n_minus_oneDt, double**** Dphi_theta_n_minus_oneDt,
    double* sigma_theta, double* sigma_phi, pml* index_of_Dr, pml* index_of_Dth, pml* index_of_Dphi);

void E_update(
    double**** Er, double**** Etheta, double**** Ephi, double**** Dr, double**** Dtheta, double**** Dphi,
    int num_of_new, int num_of_old, double***** Cmatrix, double***** Fmatrix);

double E_update_iono(
    double** Sigma_cartesian, double Value_of_Er, double Value_of_Etheta, double Value_of_Ephi,
    double Value_of_NewDr, double Value_of_NewDtheta, double Value_of_NewDphi,
    double Value_of_OldDr, double Value_of_OldDtheta, double Value_of_OldDphi,
    int Flag, double* Cmatrix, double* Fmatrix);

void H_update(
    double*** Er_at_nDt, double*** Etheta_at_nDt, double*** Ephi_at_nDt,
    double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_at_n_minus_halfDt);

void H_update_pml(
    double*** Er_at_nDt, double*** Etheta_at_nDt, double*** Ephi_at_n_Dt,
	double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_at_n_minus_halfDt,
	double**** Hr_theta1_at_n_minus_halfDt, double**** Hr_theta2_at_n_minus_halfDt, double**** Hr_phi_at_n_minus_halfDt,
	double**** Htheta_phi_at_n_minus_halfDt, double**** Htheta_r_at_n_minus_halfDt,
	double**** Hphi_r_at_n_minus_halfDt, double**** Hphi_theta_at_n_minus_halfDt,
	double* Sigma_theta_half, double* Sigma_phi_half, pml* index_of_Hr, pml* index_of_Hth, pml* index_of_Hphi);

void Ne_allocate(double *Electron_Density, double* Electron_Temperature);

void ny_allocate(date, geocoordinate, double* Colision_frequency, double* Electron_Temprature);

void geo_mag(double* Geomagnetic_in_geographic, double* Geomagnetic_in_spherical);

std::complex <double> surface_impe(std::complex <double> comp);

void surface_H_update(double **newEr, double **newEth, double **newEph, 
                      double **Hth_r0, double **Hphi_r0, double impedance_real, double impedance_image);

void PML_field_initialize(double**** Dr_theta1, double**** Dr_theta2, double**** Dr_phi,
                          double**** Dtheta_phi, double**** Dtheta_r, double**** Dphi_r, double**** Dphi_theta,
                          double**** Hr_theta1, double**** Hr_theta2, double**** Hr_phi,
                          double**** Htheta_phi, double**** Htheta_r, double**** Hphi_r, double**** Hphi_theta);

void PML_idx_initialize(pml*, pml*, pml*, pml*, pml*, pml*);

void set_matrix(
    std::complex <double> zj, double***** C_matrix, double***** F_matrix,
    double*** Iono_density, double* Collision_frequency);

void set_perturbation(perturbation P_infomation, double*** Disturbance_ion, double* Nh);

void obs_ini(geocoordinate*, geocoordinate**, int Num_of_observation_point);

void iri_profile(date, geocoordinate, double* Electron_Density, double* Electron_Tempreture);

//inline function//
inline double dist(double i){return R0 + i*delta_r;};
inline double th(double j){return THETA0 + j*delta_theta;};
inline double ph(double k){return k*delta_phi;};

inline double C_1(double sig){return ((inv_Dt - sig/2.0)/(inv_Dt + sig/2.0));};
inline double C_2(double r, double sig){return 1.0/r/delta_theta/(inv_Dt + sig/2.0);};
inline double C_3(double r, double theta){return Dt*std::cos(theta)/std::sin(theta)/2.0/r;};
inline double C_4(double r, double theta, double sig){return 1.0/r/std::sin(theta)/delta_phi/(inv_Dt + sig/2.0);};
inline double C_5(double r){return Dt/r/delta_r;};
inline double C_6(double r, double sig){return 1.0/(inv_Dt + sig/2.0)/r/delta_theta;};

#endif  // FDTD_H_ //