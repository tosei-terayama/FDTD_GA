#include <iostream>
#include <fstream>

#include "perturbation.h"
#include "fdtd3d.h"

void output_profile(perturbation P_info, double* Nh, double*** noise_Nh){
    
    std::ofstream ofs_Nh;
    ofs_Nh.open("./profile/Nh.dat");
    std::ofstream ofs_noiseNh;
    ofs_noiseNh.open("./profile/Nh_noise.dat");

    std::ofstream ofs_noiseNh_r_th;
    ofs_noiseNh_r_th.open("./profile/Nh_r_th.dat");
    std::ofstream ofs_noiseNh_th_phi;
    ofs_noiseNh_th_phi.open("./profile/Nh_th_phi.dat");
    std::ofstream ofs_noiseNh_phi_r;
    ofs_noiseNh_phi_r.open("./profile/Nh_phi_r.dat");

    int p_r0 = (int)P_info.r0();
    int p_th0 = (int)P_info.th0();
    int p_phi0 = (int)P_info.phi0();
    int lower_r = (int)Alt_lower_ionosphere/1.0e3;

    // original density profile //
    for(int i = 0; i < ion_L; i++){
        ofs_Nh << Nh[i] << " " << i + lower_r << std::endl;
    }

    for(int i = 0; i < ion_L; i++){
        for(int k = 0; k < Nphi; k++){
            ofs_noiseNh_phi_r << k << " " << i + lower_r << " " << noise_Nh[i][p_th0][k] << std::endl;
        }
        ofs_noiseNh_phi_r << std::endl;
    }

    for(int i = 0; i < ion_L; i++){
        for(int j = 0; j < Ntheta; j++){
            ofs_noiseNh_r_th << j << " " << i + lower_r << " " << noise_Nh[i][j][p_phi0] << std::endl;
        }
        ofs_noiseNh_r_th << std::endl;
    }

    for(int j = 0; j < Ntheta; j++){
        for(int k = 0; k < Nphi; k++){
            ofs_noiseNh_th_phi << k << " " << j << " " << noise_Nh[p_r0 - lower_r][j][k] << std::endl;
        }
        ofs_noiseNh_th_phi << std::endl;
    }

    for(int i = 0; i < ion_L; i++){
        ofs_noiseNh << noise_Nh[i][p_th0][p_phi0] << " " << i + lower_r << std::endl;
    }
    
    ofs_Nh.close();
    ofs_noiseNh_phi_r.close();
    ofs_noiseNh_r_th.close();
    ofs_noiseNh_th_phi.close();
    ofs_noiseNh.close();
}