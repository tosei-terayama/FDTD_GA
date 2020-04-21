#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include "GA_agent.h"

extern const double rnd_max;

void sort_Individual(int Num, double* score, bool* chromosome){
    double score_tmp;
    bool chrom_tmp;
    for(int i = 0; i < Num - 1; i++){
        for(int j = i + 1; j < Num; j++){
            if(score[i] < score[j]){
                /* swap score */
                score_tmp = score[i];
                score[i] = score[j];
                score[j] = score_tmp;

                /* swap chromosome */
                for(int k = 0; k < Nbit_total; k++){
                    chrom_tmp = chromosome[i * Nbit_total + k];
                    chromosome[i * Nbit_total + k] = chromosome[j * Nbit_total + k];
                    chromosome[j * Nbit_total + k] = chrom_tmp;
            }
        }    
    }
}

double* make_roulette(int Num, double* score){
    double sum(0.0);
    double* roulette = new double[Num];

    for(int i = 0; i < Num; i++){
        sum += score[i];
    }

    for(int i = 0; i < Num; i++){
        roulette[i] = score[i]/sum;
    }

    return roulette;
}

void Cross_over(int Head_idx, int *Ind_idx, bool* Parent_chrom, bool* Child_chrom){
    std::random_device seed;
    std::mt19937 engine( seed() );

    for(int i = 0; i < Nbit_total; i++){
        if( engine()/rnd_max > 0.5){
            /*
                cross over
            */
        }
        else {
            /*
                copy
            */
        }

    }
}