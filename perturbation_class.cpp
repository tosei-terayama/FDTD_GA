#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "perturbation.h"

void perturbation::set_geo(double lati, double longi, double alti){
    Lati = lati;
    Longi = longi;
    Alti = alti;
}

void perturbation::set_center(int r, int th, int phi){
    P_r0 = r;
    P_th0 = th;
    P_phi0 = phi;
}

void perturbation::set_alpha(double alpha){
    N0 = alpha;
}

void perturbation::set_sigma(double sig_r, double sig_h){
    Sigma_r = sig_r;
    Sigma_h = sig_h;
}