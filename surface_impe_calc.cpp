#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>

#include "fdtd3d.h"

std::complex <double> surface_impe(std::complex <double> zj)
{
  double conduct = SIGMA_SEA;

    std::complex <double> z = Z0*std::pow(EPSR - zj*conduct/EPS0/omega, -0.5);

    std::cout << "Conductivity : " << conduct << std::endl;
    
    return z;
}
