#ifndef GA_AGENT_H_
#define GA_AGENT_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

double fdtd3d(void);

void sort_Individual(int Num_Individual, double* score, bool* chromosome);

double* make_roulette(int Num_Individual, double* score);

void Cross_over(int Head_Index, int* Index_of_Individual, 
                bool* Child_Chromosome, bool* Parent_Chromosome);

void Mutation(int Num_of_Elete, int Num_of_Individual,
                double Incidence_of_Mutation, bool* Child_Chromosome);

int b2i(bool* Binary_Array, int Length_of_Array);

double b2d(bool* Binary_Array, int Length_of_Array,
            double Minimum_of_Parameter, double Parameter_Step);

constexpr int Num_patameter { 7 };

/* パラメタのビット数 (location, density, disturbance) */
//適当 後でちゃんと決める//
constexpr int Nbit_altitude { 20 };
constexpr int Nbit_latitude { Ntheta };
constexpr int Nbit_longitude { Nphi };
constexpr int Nbit_sig_r { 2 };
constexpr int Nbit_sig_th { 2 };
constexpr int Nbit_sig_ph { 2 };
constexpr int Nbit_disturbance { 2 };
constexpr int Nbit_total {
    Nbit_altitude + Nbit_latitude + Nbit_longitude +
    Nbit_sig_r + Nbit_sig_th + Nbit_sig_ph + 
    Nbit_disturbance
};

/* 適正なのかは分からず */
constexpr double param1_min { Alt_lower_ionosphere };
constexpr double param1_max { (double)Nphi };

constexpr double param2_min { 0.0 };
constexpr double param2_max { (double)Ntheta };

constexpr double param3_min { 0.0 };
constexpr double param3_max { (double)Nphi };

constexpr double param4_min { 1.0e8 };
constexpr double param4_max { 1.0e11 };

constexpr double param5_min { 1.0e8 };
constexpr double param5_max { 1.0e11 };

constexpr double param6_min { 1.0e8 };
constexpr double param6_max { 1.0e11 };

constexpr double param7_min { 1.0e8 };
constexpr double param7_max { 1.011 };

class GA_agent{
public:
    bool chrom[Nbit_total];
    double score;
    void calc_fdtd_score(void);
};

#endif