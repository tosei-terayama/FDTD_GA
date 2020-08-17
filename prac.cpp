#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>

int main(void)
{
    constexpr int Nbit_enhance { 4 };
    constexpr int Nbit_alt{ 3 };
    constexpr int Nbit_th{ 7 };
    constexpr int Nbit_phi{ 10 };
    constexpr int Nbit_sigr{ 3 };
    constexpr int Nbit_sigth{ 6 };
    constexpr int Nbit_sigphi{ 6 };
    constexpr int Nbit_total
    { Nbit_enhance + Nbit_alt + Nbit_th
    + Nbit_th + Nbit_sigr + Nbit_sigth + Nbit_sigphi };
    constexpr double rnd_max { std::pow(2, 32) };
    int Num_Individual { 20 };

    int Individual[2][Num_Individual];
    int chromosome[2][Num_Individual*Nbit_total];

    std::random_device seed;
    std::mt19937 engine( seed() );

    for(int i = 0; i < Num_Individual; i++){
        for(int j = 0; j < Nbit_total; j++){
            if( engine()/rnd_max < 0.5){
                Individual
            }
        }
    }


    return 0;
}