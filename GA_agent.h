#ifndef GA_AGENT_H_
#define GA_AGENT_H_

double fdtd3d(void);

void sort_Individual(int Num_Individual, double* score, bool* chromosome);
double* make_roulette(int Num_Individual, double* score);
void Cross_over(int Head_Index, int* Index_of_Individual, 
                bool* Child_Chromosome, bool* Parent_Chromosome)

constexpr int Num_patameter { 7 };

/* パラメタのビット数 (location, density, disturbance) */
//適当 後でちゃんと決める//
constexpr int Nbit_altitude { 20 };
constexpr int Nbit_latitude { Ntheta };
constexpr int Nbit_longitude { Nphi };
constexpr int Nbit_sigx { 2 };
constexpr int Nbit_sigy { 2 };
constexpr int Nbit_sigz { 2 };
constexpr int Nbit_disturbance { 2 };
constexpr int Nbit_total {
    Nbit_altitude + Nbit_latitude + Nbit_longitude +
    Nbit_sigx + Nbit_sigy + Nbit_sigz +
    Nbit_disturbance
};

/* 適当　後で決める*/
/*constexpr double param1_min { 0.0 };
constexpr double param1_max { 1.0 };

constexpr double param2_min { 0.0 };
constexpr double param2_max { 1.0 };

constexpr double param3_min { 0.0 };
constexpr double param3_max { 1.0 };

constexpr double param4_min { 0.0 };
constexpr double param4_max { 1.0 };

constexpr double param5_min { 0.0 };
constexpr double param5_max { 1.0 };

constexpr double param6_min { 0.0 };
constexpr double param6_max { 1.0 };

constexpr double param7_min { 0.0 };
constexpr double param7_max { 1.0 };*/

class GA_agent{
public:
    bool chrom[Nbit_total];
    double score;
    void calc_fdtd_score(void);
};

#endif