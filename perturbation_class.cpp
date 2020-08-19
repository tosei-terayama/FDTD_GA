#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "perturbation.h"

void perturbation::set_alpha(double alpha){
    N0 = alpha;
}

void perturbation::set_alt(double r){
    P_r0 = r;
}

void perturbation::set_th(double th){
    P_th0 = th;
}

void perturbation::set_phi(double phi){
    P_phi0 = phi;
}

void perturbation::set_sigr(double sig_r){
    Sigma_r = sig_r;
}

void perturbation::set_sigh(double sig_h){
    Sigma_h = sig_h;
}

void perturbation::set_geo(double lati, double longi, double alti){
    Lati = lati;
    Longi = longi;
    Alti = alti;
}

void perturbation::set_center(int r, int th, int phi){
    set_alt(r);
    set_th(th);
    set_phi(phi);
}

void perturbation::set_sigma(double sig_r, double sig_h){
    set_sigr(sig_r);
    set_sigh(sig_h);
}