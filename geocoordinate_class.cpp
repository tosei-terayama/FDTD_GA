#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"
#include "geocoordinate.h"

//GEOCOORDINATE class//
void geocoordinate::geo_ijk(double lati, double longi, double alti){
    double unit_lati = (2.0 * M_PI * R0 / 1000.0) * (lati - 135.0)/360.0;    /* 緯度1°の距離 (えびの起点) */

    //I = ;
    //J = ;
    Geo_K = alti/1.0e3;
}

void geocoordinate::set_obs(int i, int j, int k){
    Obs_I = i;
    Obs_J = j;
    Obs_K = k;
}