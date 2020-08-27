#include <iostream>
#include "fdtd3d.h"

void delete_4d(double**** array, int m, int n, int o)
{
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < o; k++){
                delete[] array[i][j][k];
            }
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    
    delete[] array;

}