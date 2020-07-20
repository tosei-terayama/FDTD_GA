#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include "fdtd3d.h"

std::complex <double> surface_impe(std::complex <double> zj)
{
  double conduct = SIGMA_VERY_DRY_GROUND;
  //double conduct = SIGMA_DRY_GROUND;

  std::complex <double> z = Z0/std::sqrt(EPSR - (zj*conduct/EPS0/omega));
  //std::complex <double> z = Z0*std::pow(EPSR - (zj*conduct/(EPS0*omega)), -0.5);

  std::cout << z << std::endl;

  std::cout << "Conductivity : " << conduct << std::endl;
    
  return z;
}
