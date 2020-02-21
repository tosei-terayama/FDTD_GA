#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

double*** memory_allocate3d(int m, int n, int o, double ini)
{
  double*** array;
  array = new double**[m];
  for(int i = 0; i < m; i++){
    array[i] = new double*[n];
    for(int j = 0; j < n; j++){
      array[i][j] = new double[o];
      for(int k = 0; k < o; k++){
        array[i][j][k] = ini;
      }
    }
  }
  return array;
}



