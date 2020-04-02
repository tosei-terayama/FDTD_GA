#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include "pml.h"

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


//The number of R, Theta, Phi element//
extern const int Nr;
extern const int Ntheta;
extern const int Nphi;

//Minute R, Theta, Phi, Time//
extern const double delta_r;
extern const double delta_theta;
extern const double delta_phi;
extern const double Dt;
extern const double inv_Dt;

//PML information//
extern const int L;
extern const double M;
extern const double R;
extern const double sigma_th_max;
extern const double sigma_phi_max;

//Ionosphere info//
extern const int ion_L;
extern const double freq;
extern const double omega;
extern const double B_abs;
extern const double Dec;
extern const double Inc;
extern const double Azim;

//function in main//
double** memory_allocate2d(int, int, double);
double*** memory_allocate3d(int, int, int, double);
double**** memory_allocate4d(int, int, int, int, double);
std::complex <double>*** memory_allocate3cd(int, int, int, std::complex<double>);

void sigma_calc(double* sigma_theta, double* sigma_phi, 
                double* sigma_theta_h, double* sigma_phi_h);

void D_update(double**** Dr, double**** Dtheta, double**** Dphi,
              double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_n_minus_halfDt,
              int num_of_new, int num_of_old);

void D_update_pml(double*** Dr_at_nDt, double*** Dtheta_at_nDt, double*** Dphi_at_nDt,
		              double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_n_minus_halfDt,
		              double**** Dr_theta1_n_minus_oneDt, double**** Dr_theta2_n_minus_halfDt, double**** Dr_phi_n_minus_oneDt,
		              double**** Dtheta_phi_at_n_minus_oneDt, double**** Dtheta_r_at_n_minus_oneDt,
		              double**** Dphi_r_at_n_minus_oneDt, double**** Dphi_theta_n_minus_oneDt,
                  double* sigma_theta, double* sigma_phi, pml* index_of_Dr, pml* index_of_Dth, pml* index_of_Dphi);

void E_update(double**** Er, double**** Etheta, double**** Ephi, double**** Dr, double**** Dtheta, double**** Dphi,
              int num_of_new, int num_of_old, double*** Sigma_cartesian, double*** Sigma_cartesian_r,
              double**** Cmatrix_in_r_comp, double**** Fmatrix_in_r_comp, double**** Cmatrix_in_th_comp, double**** Fmatrix_in_th_comp,
              double**** Cmatrix_in_phi_comp, double**** Fmatrix_in_phi_comp);

double E_update_iono(double** Sigma_cartesian, double Value_of_Er, double Value_of_Etheta, double Value_of_Ephi,
                     double Value_of_NewDr, double Value_of_NewDtheta, double Value_of_NewDphi,
                     double Value_of_OldDr, double Value_of_OldDtheta, double Value_of_OldDphi,
                     int Flag, double* Cmatrix, double* Fmatrix);

void H_update(double*** Er_at_nDt, double*** Etheta_at_nDt, double*** Ephi_at_nDt,
              double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_at_n_minus_halfDt);

void H_update_pml(double*** Er_at_nDt, double*** Etheta_at_nDt, double*** Ephi_at_n_Dt,
		              double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_at_n_minus_halfDt,
		              double**** Hr_theta1_at_n_minus_halfDt, double**** Hr_theta2_at_n_minus_halfDt, double**** Hr_phi_at_n_minus_halfDt,
		              double**** Htheta_phi_at_n_minus_halfDt, double**** Htheta_r_at_n_minus_halfDt,
		              double**** Hphi_r_at_n_minus_halfDt, double**** Hphi_theta_at_n_minus_halfDt,
		              double* Sigma_theta_half, double* Sigma_phi_half, pml* index_of_Hr, pml* index_of_Hth, pml* index_of_Hphi);

void Ne_allocate(double* Electron_density, double* Electron_Density_half, 
                  double* Electron_Temperature, double* Electron_Temperature_half);

/*void Ne_generate(double* Electron_density, double* Electron_density_half,
                  double* Electron_Temperature, double* Electron_Temperature_half);*/

void ny_allocate(double* Colision_frequency, double* Colision_frequency_half,
                  double* Electron_Temerature, double* Electron_Temperature_half);

void sig_real_calc(double* Electron_density, double* Colision_Frequency, 
                    double* Electron_density_half, double* Colision_Frequency_half,
                    double*** Sigma_realpart, double*** Sigma_realpart_r);

void make_rot_mat(double** MultipledR2_R1, double** InverseR1_2, double B_theta, double B_phi);
 
void sig_car_calc(double*** Sigma_cartesian, double*** Sigma_real, double** R1, double** R2);

void coordinate_trans(double**** Cmatrix_in_rcomp, double**** Fmatrix_in_r_comp, double**** Cmatrix_in_thcomp, double**** Fmatrix_in_thcomp,
                      double**** Cmatrix_in_phicomp, double**** Fmatrix_in_phicomp, double*** Sigma_cartesian, double*** Sigma_cartesian_r_comp);

void geo_mag(double* Geomagnetic_in_geographic, double* Geomagnetic_in_spherical);

std::complex <double> surface_impe(std::complex <double> comp);

void surface_H_update(double **newEr, double **newEth, double **newEph, 
                      double **Hth_r0, double **Hphi_r0, double impedance_real, double impedance_image);

void PML_field_initialize(double**** Dr_theta1, double**** Dr_theta2, double**** Dr_phi,
                          double**** Dtheta_phi, double**** Dtheta_r, double**** Dphi_r, double**** Dphi_theta,
                          double**** Hr_theta1, double**** Hr_theta2, double**** Hr_phi,
                          double**** Htheta_phi, double**** Htheta_r, double**** Hphi_r, double**** Hphi_theta);

void PML_idx_initialize(pml*, pml*, pml*, pml*, pml*, pml*);


//inline function//
inline double dist(double i){return R0 + i*delta_r;};
inline double th(double j){return THETA0 + j*delta_theta;};

inline double C_1(double sig){return ((inv_Dt - sig/2.0)/(inv_Dt + sig/2.0));};
inline double C_2(double r, double sig){return 1.0/r/delta_theta/(inv_Dt + sig/2.0);};
inline double C_3(double r, double theta){return Dt*std::cos(theta)/std::sin(theta)/2.0/r;};
inline double C_4(double r, double theta, double sig){return 1.0/r/std::sin(theta)/delta_phi/(inv_Dt + sig/2.0);};
inline double C_5(double r){return Dt/r/delta_r;};
inline double C_6(double r, double sig){return 1.0/(inv_Dt + sig/2.0)/r/delta_theta;};
