#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include "fdtd3d.h"

std::complex<double>** memory_allocate2cd(int m, int n, std::complex <double> ini)
{
    std::complex<double>** array;
    array = new std::complex<double> *[m];

    for(int i = 0; i < m; i++){
        array[i] = new std::complex <double>[n];
        for(int j = 0; j < n; j++){
            array[i][j] = ini;
        }
    }
    return array;
}