#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "GA_agent.h"

double calc_score(double *Magnitude, double* Target_Magnitude, int Num)
{
    double sum{ 0.0 };
    double score{ 0.0 };
    double diff{ 0.0 };

    for(int i = 0; i < Num; i++){
        diff = Magnitude[i] - Target_magnitude[i];
        sum += std::pow( diff , 2.0);
    }

    score = 1.0 / sum;
    
    return score;

}
