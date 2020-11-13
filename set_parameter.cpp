#define _USE_MATH_DEFINES
#include <cmath>
#include "GA_agent.h"
#include "perturbation.h"

void set_parameter(perturbation* P_info, bool* chrom){

    bool bit1[Nbit_enhance];
    bool bit2[Nbit_alt];
    bool bit3[Nbit_th];
    bool bit4[Nbit_phi];
    bool bit5[Nbit_sigr];
    bool bit6[Nbit_sigh];

    int count { 0 };

    for(int i = 0; i < Num_Individual; i++){

        for(int j = 0; j < Nbit_enhance; j++){
            bit1[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_alpha( b2d(bit1, Nbit_enhance, param1_min, param1_step) );
        count += Nbit_enhance;

        for(int j = 0; j < Nbit_alt; j++){
            bit2[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_alt( b2d(bit2, Nbit_alt, param2_min, param2_step) );
        count += Nbit_alt;

        for(int j = 0; j < Nbit_th; j++){
            bit3[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_th( b2d(bit3, Nbit_th, param3_min, param3_step) );
        count += Nbit_th;

        for(int j = 0; j < Nbit_phi; j++){
            bit4[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_phi( b2d(bit4, Nbit_phi, param4_min, param4_step) );
        count += Nbit_phi;

        for(int j = 0; j < Nbit_sigr; j++){
            bit5[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_sigr( b2d(bit5, Nbit_sigr, param5_min, param5_step) );
        count += Nbit_sigr;

        for(int j = 0; j < Nbit_sigh; j++){
            bit6[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_sigh( b2d(bit6, Nbit_sigh, param6_min, param6_step) );
        count = 0;

    }

    /*bool bit1[Nbit_enhance];
    bool bit2[Nbit_alt];
    bool bit3[Nbit_th];
    bool bit4[Nbit_phi];
    bool bit5[Nbit_sigr];
    bool bit6[Nbit_sigh];*/

    /*for(int i = 0; i < Num_Individual; i++){
        
        for(int j = 0; j < Nbit_enhance; j++){
            bit1[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_alpha(b2d(bit1, Nbit_enhance, param1_min, param1_step));
        count += Nbit_enhance;

        for(int j = 0; j < Nbit_alt; j++){
            bit2[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_alt(b2d(bit2, Nbit_alt, param2_min, param2_step));
        count += Nbit_alt;

        for(int j = 0; j < Nbit_th; j++){
            bit3[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_th(b2d(bit3, Nbit_th, param3_min, param3_step));
        count += Nbit_th;

        for(int j = 0; j < Nbit_phi; j++){
            bit4[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_phi(b2d(bit4, Nbit_phi, param4_min, param4_step));
        count += Nbit_phi;

        for(int j = 0; j < Nbit_sigr; j++){
            bit5[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_sigr(b2d(bit5, Nbit_sigr, param5_min, param5_step));
        count += Nbit_sigr;

        for(int j = 0; j < Nbit_sigh; j++){
            bit6[j] = chrom[i*Nbit_total + j + count];
        }
        P_info[i].set_sigh(b2d(bit6, Nbit_sigh, param6_min, param6_step));
        count = 0;
    }*/

}
