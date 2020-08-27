#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <mpi.h>

#include "GA_agent.h"
#include "fdtd3d.h"

/*Target Information
/////////////////////////////////////////////
perturbation target_Pinfo;
target_Pinfo.set_alpha(10.0);
target_Pinfo.set_center(74, Ntheta/2, Nphi/2);
target_Pinfo.set_sigma(2.0e3, 60.0e3);
/////////////////////////////////////////////
*/

constexpr int Num_Individual { 6 };  // Number of individuals
constexpr int Num_Generation { 30 };  // Number of generations to repeat
constexpr int Num_Elete { 2 };  //  Number of elete
constexpr double rnd_max { std::pow(2, 32) };  //   Max of mersenne twister (32 bit)
constexpr double Mutation_rate { 0.03 };  // Mutation incidence

int main(int argc, char** argv){

    //double total_time;
    
    std::ofstream ofs;
    std::ofstream ofs_score;
    std::ofstream ofs_score2;
    std::ofstream ofs_score3;

    ofs.open("./result/magnitude.dat");
    ofs_score.open("./result/score.dat");
    ofs_score2.open("./result/score2.dat");
    ofs_score3.open("./result/score3.dat");
    
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
    // ここでP_infoの初期値(ランダム)を設定する？ //
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

    // boolean -> parameter //
    perturbation P_info[Num_Individual];

    set_parameter(P_info, chromosome[0]);

    /* Initialize parameter ( date/perturbation/geocoordinate ) */
    date ymd;
    ymd.set_ymd(2016, 3, 1);
    ymd.set_h(9.0);

    geocoordinate lla_info;
    lla_info.set_point(32.0, 135.0, 60.0);

    // Observation points on propagation path //
    int Num_obs = (Nphi - 2*L) - k_s;
    geocoordinate* obs_p = new geocoordinate[Num_obs + 1];

    // Magnitude //
    double **Magnitude = memory_allocate2d(Num_Individual, Num_obs + 1, 0.0); 
    double *Target_Magnitude = new double[Num_obs + 1];

    std::ifstream ifs;
    ifs.open("./target.dat");

    for(int i = 0; i < Num_obs; i++){
        ifs >> Target_Magnitude[i];
    }
    ifs.close();

    int child{ 0 };
    double judge{1.0e-2};
    bool flag = false;

    /*if(rank == 0) {
        std::chrono::system_clock::time_point start
        = std::chrono::system_clock::now();
        }*/

    /* GA programming(本体) */
    for(int gen = 0; gen < Num_Generation; gen++){

        if(rank == 0){
            std::cout << gen << " Generation " << std::endl;
        }

        const int PARENT { gen % 2 };
        const int CHILD { (gen + 1) % 2 };
        child = CHILD;

        /* Calculate FDTD & Score (PE n) */
        // problem section //
        for(int i = start_idx[rank]; i < end_idx[rank]; i++){
                fdtd_calc(P_info[i], ymd, lla_info, Num_obs, obs_p, Magnitude[i], rank);
                std::cout << "Mag(0) : " << Magnitude[0] << " Mag(750) : " << Magnitude[750] << std::endl; 
                Individual[PARENT][i].score = calc_score(Magnitude[i], Target_Magnitude, Num_obs);
                score[i] = Individual[PARENT][i].score;
                std::cout << "Individual.score : " << i << " " << Individual[PARENT][i].score <<
                 "  score : " << score[i] << std::endl;
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

        if(rank == 0){
            for(int i = 0; i < Num_Individual; i++){
                if(judge > score[i]){
                    std::cout << "best score : " << i << " " << score[i] << std::endl;
                    flag = true;
                }
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

            /* 3% Mutation */
            Mutation(Num_Elete, Num_Individual, Mutation_rate, chromosome[CHILD]);

            ofs_score << gen << " " << score[0] << std::endl;
            ofs_score2 << gen << " " << score[1] << std::endl;
            ofs_score3 << gen << " " << score[2] << std::endl;
        }

        /* Assign Chromosome to each agent */
        if( rank == 0 ){
            for(int i = 1; i < size; i++){
                MPI::COMM_WORLD.Send(chromosome[CHILD], 
                    Num_Individual*Nbit_total, MPI::BOOL, i, 1);
            }
        }   else{
            MPI::COMM_WORLD.Recv(chromosome[CHILD],
                Num_Individual*Nbit_total, MPI::BOOL, 0, 1);
        }

        for(int i = start_idx[rank]; i < end_idx[rank]; i++){
            for(int j = 0; j < Nbit_total; j++){
                Individual[CHILD][i].chrom[j]
                    = chromosome[CHILD][i*Nbit_total + j];
            }
        }
        
        if( rank == 0){
            set_parameter(P_info, chromosome[CHILD]);
        }

        if(flag == true) break;

    }

    if(rank == 0){
        double best_score = Individual[child][0].score;
        int best_ind = 0;

        for(int i = 0; i < Num_Individual; i++){
            if(best_score < Individual[child][i].score){
                best_ind = i;
            }
        }

        for(int i = 0; i < Num_obs + 1; i++){
            ofs << i << " " << Magnitude[best_ind][i] << std::endl;
        }
    }

    /*if(rank == 0){
        std::chrono::system_clock::time_point end
            = std::chrono::system_clock::now();
    }*/

    MPI::Finalize();

    /*total_time = std::chrono::duration_cast <std::chrono::milliseconds>
        (end - start).count();
    
    std::cout << "elapsed time : " << total_time*1.0e-3 << " [sec]" << std::endl;*/
    ofs.close();

    return 0;
}