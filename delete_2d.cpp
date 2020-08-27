#include <iostream>
#include "fdtd3d.h"

void delete_2d(double** array, int m)
{
    for(int i = 0; i < m; i++){
        delete[] array[i];
    }

    delete[] array;
    
}