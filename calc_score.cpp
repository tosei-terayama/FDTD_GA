#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "GA_agent.h"

double calc_score(double *Magnitude, double* Target_Magnitude, int Num)
{
    double score{ 0.0 };

    for(int i = 0; i < Num; i++){
        score += std::abs(Magnitude[i] - Target_Magnitude[i]);
    }

    score = std::pow(score, 2.0);
    return score;

}