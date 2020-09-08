#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>

const double rnd_max{ std::pow(2, 32) };

void make_random(bool* A, int size){
    std::random_device seed_gen;
    std::mt19937 engine( seed_gen() );

    for(int i = 0; i < size; i++){
        if( engine()/rnd_max > 0.5 ){
            std::cout << "1 " << std::endl;
        }
        else std::cout << "2 " << std::endl;
    }

}