#include <iostream>
#include <fstream>

#include "fdtd3d.h"
#include "iri.h"

int Ne_iri(void){
    float *outf = new float[20*1000];
    float xlat = 35.0f, xlon = 135.0f;
    int iy = 2016, imd = 301;
    float hour = 3.0f;

    float alt_beg = 60.0f, alt_end = 100.0f, alt_step = 1.0f;

    iri_(&xlat, &xlon, &iy, &imd, &hour, &alt_beg, &alt_end, &alt_step, outf);
    

}