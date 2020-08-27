#include <iostream>
#include "fdtd3d.h"

void delete_3d(double*** array, int m, int n)
{
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            delete[] array[i][j];
        }
        delete[] array[i];
    }

    delete array;
}