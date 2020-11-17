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
Target_param.set_alpha( 10.0 );
Target_param.set_center(74, 25, Nphi/2);
Target_param.set_sigma(2.0e3, 30.0e3);
/////////////////////////////////////////////
*/

constexpr int Num_Individual { 24 };  // Number of individuals
constexpr int Num_Generation { 2 };  // Number of generations to repeat
constexpr int Num_Elete { 2 };  //  Number of elete
constexpr double rnd_max { std::pow(2, 32) };  //   Max of mersenne twister (32 bit)
constexpr double Mutation_rate { 0.03 };  // Mutation incidence

int main(int argc, char** argv){

    //double total_time;

    MPI::Init(argc, argv);
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();
    const int assigned_num = Num_Individual / size; // Assignement to processor
    int name_length = 256;
    char* name = new char[name_length];
    MPI::Get_processor_name(name, name_length);

      std::ofstream ofs;
      std::ofstream ofs_gen;
      std::ofstream ofs_param;
      std::ofstream ofs_score0;
      std::ofstream ofs_score1;
      std::ofstream ofs_score2;
      std::ofstream ofs_score3;
      std::ofstream ofs_score4;
      std::ofstream ofs_score5;
      std::ofstream ofs_ave;
      std::ofstream ofs_test;

    if( rank == 0 ){
      ofs.open("./result/magnitude.dat");
      ofs_param.open("./result/param.dat");
      ofs_score0.open("./result/score0.dat");
      ofs_score1.open("./result/score1.dat");
      ofs_score2.open("./result/score2.dat");
      ofs_score3.open("./result/score3.dat");
      ofs_score4.open("./result/score4.dat");
      ofs_score5.open("./result/score5.dat");
      ofs_ave.open("./result/score_average.dat");
      ofs_test.open("./result/test.dat");
    }
    
    std::random_device seed;
    std::mt19937 engine( seed() );    //mersenne twister engine

    GA_agent Individual[2][Num_Individual];

    /* chromosome & score array for sorting */
    bool chromosome[2][Num_Individual * Nbit_total];
    double* score = new double[Num_Individual];

    /* Initialize chromosomes */
    for( int i = 0; i < Num_Individual; i++ ){
        for( int j = 0; j < Nbit_total; j++ ){
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

    if(rank == 0){
        for(int i = 0; i < Num_Individual; i++){
            std::cout << i << " :   ";
            for(int j = 0; j < Nbit_total; j++){
                std::cout << chromosome[0][i * Nbit_total + j] << " ";
            }
            std::cout << std::endl;
        }
        set_parameter(P_info, chromosome[0]);

        for(int i = 0; i < Num_Individual; i++){
            std::cout << i << " " << P_info[i].alpha() << " " << P_info[i].r0() << " " << P_info[i].th0() << " " << P_info[i].phi0()
                    << " " << P_info[i].sig_r() << " " << P_info[i].sig_h() << std::endl;
        }
    }

    /* Initialize parameter ( date/perturbation/geocoordinate ) */
    date ymd;
    ymd.set_ymd(2016, 3, 1);
    ymd.set_h(9.0);

    geocoordinate lla_info;
    lla_info.set_point(32.0, 135.0, (Alt_lower_ionosphere/1.0e3) );

    // Observation points on propagation path //
    int Num_obs = (Nphi - 2*L) - k_s + 1;
    geocoordinate* obs_p = new geocoordinate[Num_obs];
    for( int k = 0; k < Num_obs; k++ ){
        obs_p[k].set_obs(0, 50, k + k_s);
    }

    // Magnitude //
    double *Target_Magnitude = new double[Num_obs];
    double **Magnitude = memory_allocate2d(Num_Individual, Num_obs, 0.0);

    std::ifstream ifs;
    ifs.open("./target.dat");

    int buff;

    for( int i = 0; i < Num_obs; i++ ){
        ifs >> buff;
        ifs >> Target_Magnitude[i];
    }
    ifs.close();

    if( rank == 0 ) {
        ofs_param << " # Ind   alpha  r  the  phi  sigma_r   sigma_h   score #" << std::endl;
    }

    double judge{ 1.0e3 };
    bool flag = false;

     std::cout << name << "  rank : " << rank << std::endl;
     ofs_test << name << " " << rank << std::endl;

     ofs_test.close();

    /* GA programming(本体) */
    for(int gen = 0; gen <= Num_Generation; gen++){

        if(rank == 0){
            std::cout << gen << " Generation " << std::endl;
            ofs_param << " #   " << gen << " generation   #" << std::endl;
        }

        const int PARENT { gen % 2 };
        const int CHILD { (gen + 1) % 2 };

        /* Calculate FDTD & Score (PE n) */
        // problem section //
        for(int i = start_idx[rank]; i < end_idx[rank]; i++){

          fdtd_calc(P_info[i], ymd, lla_info, Num_obs, obs_p, Magnitude[i]);
          score[i] = calc_score(Magnitude[i], Target_Magnitude, Num_obs, i);
          Individual[PARENT][i].score = score[i];

        }

        std::cout << rank << " is completed. " << std::endl;

        /* Merging scores */
        if( rank != 0){
            MPI::COMM_WORLD.Send(score + start_idx[rank], assigned_num,
                                MPI::DOUBLE, 0, 0);
            }
        else{  /* rank0 : 計算結果の受信 */
            for(int i = 1; i < size; i++){
                MPI::COMM_WORLD.Recv(score + start_idx[i], assigned_num, 
                                    MPI::DOUBLE, i, 0);
            }
        }

        if( rank == 0 ){
            for( int i = 0; i < Num_Individual; i++ ) std::cout << i << "  " << score[i] << " " << Individual[PARENT][i].score << std::endl;
        }

        /* Sync All Process */
        //MPI::COMM_WORLD.Barrier();

        if(rank == 0){
            std::string fn = "./result/gen" + std::to_string(gen) + ".dat";
            ofs_gen.open(fn);
            std::ofstream ofs_gen(fn.c_str());

            ofs_gen << gen << " generation. " << std::endl;
            ofs_gen.close();

            for( int i = 0; i < Num_Individual; i++ ){
                ofs_param << i << "  " << P_info[i].alpha() << "   " << P_info[i].r0() << "   " << P_info[i].th0() << "   "
                            << P_info[i].phi0() << "   " << P_info[i].sig_r() << "   " << P_info[i].sig_h() << std::endl;
            }
            ofs_param << std::endl;

            double score_ave{ 0.0 };
            for(int i = 0; i < Num_Individual; i++) score_ave += score[i];
            score_ave = score_ave/Num_Individual;
            ofs_ave << gen << " " << score_ave << std::endl; 
        }

        //MPI::COMM_WORLD.Barrier();

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
                    double rnd_num = engine() / rnd_max;
                    while( Roulette[k] < rnd_num){
                        k++;
                    }
                    ind_idx[j] = k;
                }

                /* Uniform cross over */
                Cross_over(i, ind_idx, chromosome[PARENT], chromosome[CHILD]);
            }

            /* 3% Mutation */
            Mutation(Num_Elete, Num_Individual, Mutation_rate, chromosome[CHILD]);

            for(int i = 0; i < Num_Individual; i++){
                std::cout << i << " score :  " << score[i] << std::endl;
            }

            if( rank == 0 && gen == Num_Generation ){
                for( int i = 0; i < Num_obs; i++ ){
                    ofs << i << " " << Magnitude[0][i] << std::endl;
                }
            }

            ofs_score0 << gen << " " << score[0] << std::endl;
            ofs_score1 << gen << " " << score[1] << std::endl;
            ofs_score2 << gen << " " << score[2] << std::endl;
            ofs_score3 << gen << " " << score[3] << std::endl;
            ofs_score4 << gen << " " << score[4] << std::endl;
            ofs_score5 << gen << " " << score[5] << std::endl;
            std::cout << "GA complete" << std::endl;
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

        /* Sync All Process */
        //MPI_Barrier(MPI_COMM_WORLD);

        for(int i = start_idx[rank]; i < end_idx[rank]; i++){
            for(int j = 0; j < Nbit_total; j++){
                Individual[CHILD][i].chrom[j]
                    = chromosome[CHILD][i*Nbit_total + j];
            }
        }

        if(rank == 0 && score[0] >= judge) flag = true;

        if(flag == true) break;

        set_parameter(P_info, chromosome[CHILD]);

    }

    MPI::Finalize();

    /*set_parameter(P_info, chromosome[child]);

    for(int i = 0; i < Num_Individual; i++) {
        std::cout << "///////////////////////////////////////////////" << std::endl; 
        std::cout << i << "alpha : " << P_info[i].alpha() << " r0 : " << P_info[i].r0()
                << " th0 : " << P_info[i].th0() << " phi0 : " << P_info[i].phi0() << std::endl
                << " sig_r : " << P_info[i].sig_r() << " sig_h : " << P_info[i].sig_h() << std::endl;
        std::cout << "///////////////////////////////////////////////" << std::endl;
        }*/

    ofs.close();
    ofs_param.close();
    ofs_score0.close();
    ofs_score1.close();
    ofs_score2.close();
    ofs_score3.close();
    ofs_score4.close();
    ofs_score5.close();

    return 0;
}





