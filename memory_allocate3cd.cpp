#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>

std::complex<double>*** memory_allocate3cd(int m, int n, int o, std::complex<double> ini)
{

  std::complex<double>*** array;
  array = new std::complex<double> **[m];
  for(int i = 0; i < m; i++){
    array[i] = new std::complex<double> *[n];
    for(int j = 0; j < n; j++){
      array[i][j] = new std::complex<double> [o];
      for(int k = 0; k < o; k++){
        array[i][j][k] = ini;
      }
    }
  }
  return array;  
}
