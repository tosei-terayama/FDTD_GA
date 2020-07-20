#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "perturbation.h"

void perturbation::set_geo(double lati, double longi, double alti){
    Lati = lati;
    Longi = longi;
    Alti = alti;
}

void perturbation::set_center(double r, double th, double phi){
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

void perturbation::set_range(int range_r, int range_th, int range_phi){
    Range_r = range_r;
    Range_th = range_th;
    Range_phi = range_phi;
}