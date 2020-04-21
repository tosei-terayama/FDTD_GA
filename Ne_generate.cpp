#include <iostream>
#include <fstream>

#include "./IRI-2016_for_cpp/iri.h"
#include "fdtd3d.h"

void Ne_generate(double* Nh, double* Nh_h, double* Re, double* Re_h)
{
    float *outf = new float [20*1000];
    float xlat = 35.0f, xlon = 135.0f;  /* latitude & Longitude */
    int iy = 2016, imd = 229;  /* Year -> 2016  Month/Date -> 2/29 in Grenidge */
    float hour = 21.5f; /* Universal time ( LT 06:30 )*/

    /* Altitude : 60 ~ 100 [km]  every 0.5 km */
    float alt_beg = 60.0f, alt_end = 100.0f, alt_step = 0.5;

    iri_(&xlat, &xlon, &iy, &imd, &hour, &alt_beg, &alt_end, &alt_step, outf);

    for(int i = 0; i < ion_L; i += 2){
        Nh[i] = outf[i]*1.0e6;
        Nh_h[i] = outf[i + 1]*1.0e6;

        Re[i] = outf[3000 + i]/300.0;
        Re_h[i] = outf[3000 + i + 1]/300.0;
    }
}