#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>

int main(void){
    std::ifstream ifs("./alt_study/amb_alt80dB.dat");
    double max(0.0);
    double min(0.0);
    int trash;
    int max_i;
    int min_i;
    double normal(0.0);
    double joran(0.0);
    double value;

    ifs >> trash >> normal >> joran;
    ifs >> trash >> normal >> joran;
    
    max = std::abs(joran/normal);
    min = max;
    max_i = 0;
    min_i = 0;
    
    for(int i = 1; i < 1000; i++){
        ifs >> trash >> normal >> joran;
        //std::cout << trash << " " << normal << " " << joran << std::endl;
        value = std::abs(joran/normal);
        if(max <= value){
            max_i = i;
            max = value;
        }
        if(min >= value){
            min_i = i;
            min = value;
        }
    }

    std::cout << "max_i = " << max_i << " max_value = " << max << std::endl;
    std::cout << "min_i = " << min_i << " min_value = " << min << std::endl;

    return 0; 
}