#include "fdtd3d.h"
#include "geocoordinate.h"

void obs_ini(geocoordinate* obs2d, geocoordinate** obs3d, int num_obs){

    for(int k = 0; k <= num_obs; k++){
        obs2d[k].set_obs(0, 50, k + k_s);
    }

    for(int j = 0; j <= Ntheta; j++){
        for(int k = 0; k <= num_obs; k++){
            obs3d[j][k].set_obs(0, j, k + k_s);
        }
    }

}