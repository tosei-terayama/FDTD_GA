#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include <mpi.h>

#include "GA_agent.h"

constexpr int Num_Individual { 56 };  // Number of individuals
constexpr int Num_Generation { 30 };  // Number of generations to repeat
constexpr int Num_Elete { 2 };  //  Number of elete
constexpr double rnd_max { std::pow(2, 32) };  //   Max of mersenne twister (32 bit)
constexpr double Mutation { 0.03 };  // Mutation rate

int main(int argc, char** argv){

    /* Target score */
    /*
     FDTDでスコア生成 
    */
    
    MPI::Init(argc, argv);
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();
    const int assigned_num = Num_Individual / size; // Assignement to processor

    std::random_device seed;
    std::mt19937 engine( seed() );    //mersenne twister engine

    GA_agent Individual[2][Num_Individual];

    /* chromosome & score array for sorting */
    bool chromosome[2][Num_Individual * Nbit_total];
    double* score = new double[Num_Individual];

    /* Initialize chromosomes */
    for(int i = 0; i < Num_Individual; i++){
        for(int j = 0; j < Nbit_total; j++){
            if( engine()/rnd_max < 0.5 ) {
                Individual[0][i].chrom[j] = true;
            }
            else {
                Individual[0][i].chrom[j] = false;
            }

            chromosome[0][i*Nbit_total + j] = Individual[0][i].chrom[j];
        }
    }

    /* Assignment individual to each processor */
    int* start_idx = new int[size];
    int* end_idx = new int[size];
    for(int myrank = 0; myrank < size; myrank++){
        start_idx[myrank] = myrank * assigned_num;
        end_idx[myrank] = (myrank + 1) * assigned_num;
    }

    /* GA programming */
    for(int gen = 0; gen < Num_Generation; gen++){
        const int PARENT { gen % 2 };
        const int CHILD { (gen + 1) % 2 };

        /* Calculate FDTD & Score (PE n)*/
        for(int i = start_idx[rank]; i < end_idx[rank]; i++){
            /*
                FDTD simulate
            */
           score[i] = Individual[PARENT][i].score;
        }

        /* Merging scores */
        if( rank != 0){
            for(int myrank = 1; myrank < size; myrank++){
                MPI::COMM_WORLD.Send(score + start_idx[myrank], assigned_num,
                                    MPI::DOUBLE, 0, 0);
                MPI::COMM_WORLD.Recv(score + start_idx[myrank], assigned_num,
                                    MPI::DOUBLE, myrank, 0);
            }
        }

        /* Genetic Algorithm */
        if(rank == 0){
            /* Sort */
            sort_Individual(Num_Individual, score, chromosome[PARENT]);

            /* Elete strategy */
            for(int i = 0; i < Num_Elete; i++){
                for(int j = 0; j < Nbit_total; j++){
                    chromosome[CHILD][i*Nbit_total + j]
                        = chromosome[PARENT][i*Nbit_total + j];
                }
            }

            /* Make roulette */
            double* Roulette = make_roulette(Num_Individual, score);

            /* Selection & Cross over ( 28×2 - 2 ) */
            for(int i = Num_Elete; i < Num_Individual; i += 2){
                int ind_idx[2];

                /* Roulette selection */
                for(int j = 0; j < 2; j++){
                    int k = 0;
                    double sum = Roulette[k];
                    double rnd_num = engine() / rnd_max;
                    while( sum < rnd_num){
                        k++;
                        sum += Roulette[k]; 
                    }
                    ind_idx[j] = k;
                }

                /* Uniform cross over */
                Cross_over(i, ind_idx, chromosome[PARENT], chromosome[CHILD]);
            }
        }

    }

    MPI::Finalize();

    return 0;
}