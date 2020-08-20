#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "GA_agent.h"

double calc_score(double *Magnitude, double* Target_Magnitude, int size)
{
    double score{ 0.0 };

    for(int i = 0; i < size; i++){
        score += std::pow(Magnitude[i] - Target_Magnitude[i], 2.0);
    }

    return score;

}