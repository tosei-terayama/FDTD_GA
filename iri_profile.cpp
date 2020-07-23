#include <iostream>
#include <fstream>
#include "./IRI-2016_for_cpp/iri.h"
#include "fdtd3d.h"

void iri_profile(date date, geocoordinate lla, double* Nh, double* Re){
    
    float *outf = new float [20*1000];
    float xlat = (float)lla.lati();
    float xlon = (float)lla.longi();
    int iy = date.iy();
    int imd = date.imd();
    float hour = (float)date.ih();

    float alt_beg = (float)lla.alt();
    float alt_end = (float)alt_beg + ion_L;
    float alt_step = 1.0f;

    iri_(&xlat, &xlon, &iy, &imd, &hour, &alt_beg, &alt_end, &alt_step, outf);

    for(int i = 0; i <= ion_L; i++){
        if(outf[i] == -1){
            Nh[i] = 0.0;
        }
        else Nh[i] = (double)outf[i]*1.0e6;

        if(outf[i] == -1){
            Re[i] = 0.0;
        }
        else Re[i] = (double)outf[3000 + i]/300.0;
    }
    
}