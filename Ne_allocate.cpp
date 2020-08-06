#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "fdtd3d.h"

void Ne_allocate(double *Nh, double *Re)
{
    double trash[100];

    double height(0.0);
    int bufsize = 256;
    char buf[bufsize];
    int s = 25;

    std::ifstream ifs("./IRI-2016_for_cpp/fort.7");

    if(!ifs){
        std::cout << "error" << std::endl;
        exit(0);
    }

    // ??? //
    for(int i = 0; i < s; i++){
        ifs.getline(buf, bufsize);
    }

    for(int i = 0; i <= ion_L; i++){
        ifs >> height >> Nh[i];
        if(Nh[i] == -1) Nh[i] = 0.0;
        else Nh[i] = Nh[i]*1.0e6;

        for(int m = 0; m < 3; m++) ifs >> trash[m];
        ifs >> Re[i];
        Re[i] = Re[i]/300.0;

        for(int m = 0; m < 24; m++) ifs >> trash[m];

    }

    ifs.close();
    
}