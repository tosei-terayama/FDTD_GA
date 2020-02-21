#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include <Eigen/Core>

#include "fdtd3d.h"

std::complex <double> surface_impe(std::complex <double> zj)
{
  double conduct = SIGMA_WET_GROUND;

    std::complex <double> z = Z0*std::pow(EPSR - zj*conduct/EPS0/omega, -0.5);

    std::cout << "Conductivity : " << conduct;
    
    return z;
}
