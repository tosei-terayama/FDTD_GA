#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>

const int Num_Individual{ 3 };
const int N_num1{ 2 };
const int N_num2{ 3 };
const int N_num3{ 4 };
const int N_num4{ 5 };
const int Nbit_total{ N_num1 + N_num2 + N_num3 + N_num4 };

const double rnd_max { std::pow(2, 32) };

void make_random(bool* A, int size);

int main(void)
{
    /*std::random_device seed;
    std::mt19937 engine( seed() );  //mersenne twister//

    bool Nbit[ Num_Individual*Nbit_total ];
    bool Nbit1[N_num1];
    bool Nbit2[N_num2];
    bool Nbit3[N_num3];
    bool Nbit4[N_num4];

    for(int i = 0; i < Num_Individual; i++){
        for(int j = 0; j < Nbit_total; j++){

            if( engine()/rnd_max < 0.5){
                Nbit[i*Num_Individual + j] = true;
            }
            else Nbit[i*Num_Individual + j] = false;
        }
    }

    int count{0};

    for(int i = 0; i < Num_Individual; i++){

        std::cout << std::endl << i << "\tNbit1: " << std::endl;
        std::cout << "____________________" << std::endl;
        for(int j = 0; j < N_num1; j++){
            Nbit1[j] = Nbit[i*Nbit_total + j];
            std::cout << j << " " << Nbit1[j] << std::endl;
            //std::cout << i*Nbit_total + j << std::endl;
        }
        count += N_num1;

        std::cout << std::endl  << "Nbit2: " << std::endl;

        for(int j = 0; j < N_num2; j++){
            Nbit2[j] = Nbit[i*Nbit_total + j + count];
            std::cout << j << " " << Nbit2[j] << std::endl;
            //std::cout << i*Nbit_total + j + count << std::endl;
        }
        count += N_num2;

        std::cout << std::endl  << "Nbit3: " << std::endl;
        for(int j = 0; j < N_num3; j++){
            Nbit3[j] = Nbit[i*Nbit_total + j + count];
            std::cout << j << " " << Nbit3[j] << std::endl;
            //std::cout << i*Nbit_total + j + count << std::endl;
        }
        count += N_num3;

        std::cout << std::endl  << "Nbit4: " << std::endl;
        for(int j = 0; j < N_num4; j++){
            Nbit4[j] = Nbit[i*Nbit_total + j + count];
            std::cout << j << " " << Nbit4[j] << std::endl;
            //std::cout << i*Nbit_total + j + count << std::endl;
        }
        std::cout << "____________________" << std::endl << std::endl;
        count = 0;
    }

    for(int i = 0; i < Num_Individual * Nbit_total; i++){
        std::cout << i << " " << Nbit[i] << std::endl;
    }*/

    int size = 12;
    bool* senpai = new bool[size];

    make_random(senpai, size);

    return 0;
}