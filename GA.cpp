#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include "GA_agent.h"

const double rnd_max{ std::pow(2, 32) };

void sort_Individual(int Num, double* score, bool* chromosome){
    double score_tmp;
    bool chrom_tmp;

    for(int i = 0; i < Num - 1; i++){
        for(int j = i + 1; j < Num; j++){
            if(score[i] < score[j]){
                // swap score //
                score_tmp = score[i];
                score[i] = score[j];
                score[j] = score_tmp;

                // swap chromosome //
                for(int k = 0; k < Nbit_total; k++){
                    chrom_tmp = chromosome[i * Nbit_total + k];
                    chromosome[i * Nbit_total + k] = chromosome[j * Nbit_total + k];
                    chromosome[j * Nbit_total + k] = chrom_tmp;
                }
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
    
    std::random_device seed_gen;
    std::mt19937 engine( seed_gen() );

    for(int j = 0; j < Nbit_total; j++){
        if( engine()/rnd_max > 0.5){
            /* Cross over */
            Child_chrom[Head_idx * Nbit_total + j]
                = Parent_chrom[Ind_idx[1] * Nbit_total + j];
            Child_chrom[(Head_idx+1) * Nbit_total + j]
                = Parent_chrom[Ind_idx[0] * Nbit_total + j];
        }

        else {
            Child_chrom[Head_idx * Nbit_total + j]
             = Parent_chrom[Ind_idx[0] * Nbit_total + j];
            Child_chrom[(Head_idx+1) * Nbit_total + j]
             = Parent_chrom[Ind_idx[1] * Nbit_total + j];
        }
    }
    
}

void Mutation(int Num_Elete, int Num_Ind, double Mutation_rate, bool* Child_chrom){
    
    std::random_device seed;
    std::mt19937 engine( seed() );

    for(int i = Num_Elete; i < Num_Ind; i++){
        for(int j = 0; j < Nbit_total; j++){
            /* Bit flip */
            if( engine()/rnd_max < Mutation_rate){
                Child_chrom[i * Nbit_total + j]
                = !Child_chrom[i * Nbit_total + j];
            }
        }
    }

}