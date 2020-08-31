#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include "fdtd3d.h"

void delete_2cd(std::complex <double>** array, int m)
{
    for(int i = 0; i < m; i++){
        delete[] array[i];
    }

    delete[] array;
    
}