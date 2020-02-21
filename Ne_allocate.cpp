#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "fdtd3d.h"

void Ne_allocate(double *Nh, double *Nh_h, double *Re, double *Re_h)
{
    double trash[13];
    double Te(0.0), Te_h(0.0);

    double height(0.0);
    int bufsize = 256;
    char buf[bufsize];  //空読み用//
    int s = 25;    //読み込み初めの行//
    int e = s + ion_L;    //読み終わりの行//
    
    std::ifstream ifs("./IRI-2016_for_cpp/fort.7");

    if(!ifs){
        std::cout << "error" << std::endl;
        exit(0);
    }

    //空読み//
    for(int i = 0; i < s; i++){
        ifs.getline(buf, bufsize);
    }

    for(int i = 0; i < ion_L; i++){

        ifs >> height >> Nh[i];

        Nh[i] = Nh[i]*1.0e6;

        for(int m = 0; m < 3; m++) ifs >> trash[m];

        ifs >> Te;

        Re[i] = Te/300.0;

        for(int m = 0; m < 9; m++) ifs >> trash[m];

        if(Nh[i] < 0.0) Nh[i] = 0.0;

        ifs >> height >> Nh_h[i];

        Nh_h[i] = Nh_h[i]*1.e6;

        for(int m = 0; m < 3; m++) ifs >> trash[m];

        ifs >> Te_h;
        
        Re_h[i] = Te_h/300.0;

        for(int m = 0; m < 9; m++) ifs >> trash[m];

        if(Nh_h[i] < 0.0) Nh_h[i] = 0.0;

    }

    ifs.close();
    
}