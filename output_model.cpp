#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include "fdtd3d.h"

void output_model(void)
{
    std::ofstream ofs_model;
    ofs_model.open("./profile/model.dat");

    double radius[2];
    double zx_angle[2];
    double xy_angle[2];

    radius[0] = dist(0.0);
    radius[1] = dist(Nr);
    zx_angle[0] = th(0.0);
    zx_angle[1] = th(Ntheta);
    xy_angle[0] = ph(0.0);
    xy_angle[1] = ph(Nphi);

    for(int s = 0; s < 4; s++){

         for(int k = 0; k <= Nphi; k++){
            ofs_model << radius[s/2]*std::sin(zx_angle[s%2])*std::cos(ph(k)) << " " << radius[s/2]*std::sin(zx_angle[s%2])*std::sin(ph(k))
                        << " " << radius[s/2]*std::cos(zx_angle[s%2]) << std::endl;
        }
        ofs_model << std::endl;

        for(int j = 0; j <= Ntheta; j++){
            ofs_model << radius[s%2]*std::sin(th(j))*std::cos(xy_angle[s/2]) << " " << radius[s%2]*std::sin(th(j))*std::sin(xy_angle[s/2])
                        << " " << radius[s%2]*std::cos(th(j)) << std::endl;
        }
        ofs_model << std::endl;

        for(int i = 0; i <= Nr; i++){
            ofs_model << dist(i)*std::sin(zx_angle[s%2])*std::cos(xy_angle[s/2]) << " " << dist(i)*std::sin(zx_angle[s%2])*std::sin(xy_angle[s/2]) 
                     << " " << dist(i)*std::cos(zx_angle[s%2]) << std::endl;
        }
        ofs_model << std::endl;

    }

    ofs_model.close();

}