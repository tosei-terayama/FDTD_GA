#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

double***** memory_allocate5d(int n, int o, int p, int q, int r, double ini)
{
    double *****array;
    array = new double****[n];

    for(int i = 0; i < n; i++){
        array[i] = new double***[o];
        for(int j = 0; j < o; j++){
            array[i][j] = new double**[p];
            for(int k = 0; k < p; k++){
                array[i][j][k] = new double*[q];
                for(int l = 0; l < q; l++){
                    array[i][j][k][l] = new double[r];
                    for(int m = 0; m < r; m++){
                        array[i][j][k][l][m] = ini;
                    }
                }
            }
        }
    }

    return array;
}