#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

#include "fdtd3d.h"

void geo_mag(double *geo_B, double *sph_B)
{
    //Geomagnetic_field_in_geographic//
    geo_B[0] = -B_abs*std::sin(Inc);
    geo_B[1] = -B_abs*std::cos(Inc)*std::cos(Dec);
    geo_B[2] = B_abs*std::cos(Inc)*std::sin(Dec);

    //Geomagnetic_field_in_spherical_coordinate//
    sph_B[0] = geo_B[0];
    sph_B[1] = B_abs*(-std::cos(Inc)*std::cos(Dec)*std::cos(Azim)
     + std::cos(Inc)*std::sin(Dec)*std::sin(Azim));
    sph_B[2] = B_abs*(std::cos(Inc)*std::cos(Dec)*std::sin(Azim)
     + std::cos(Inc)*std::sin(Dec)*std::cos(Azim));
     
}