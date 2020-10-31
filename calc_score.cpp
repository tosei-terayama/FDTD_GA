#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "GA_agent.h"

double calc_score(double *Magnitude, double* Target_Magnitude, int Num, int rank)
{
    std::ofstream ofs_check;
    ofs_check.open("./result/check.dat");
    double sum{ 0.0 };
    double score{ 0.0 };
    double diff{ 0.0 };

    for(int i = 0; i < Num; i++){
      diff = Magnitude[i] - Target_Magnitude[i];

      if( std::isnan( diff ) ){
        ofs_check << "rank : " << rank << " obs : " << i << " Mag " << Magnitude[i]
                  << " diff : " << diff << std::endl;
      }
      sum += std::pow( diff , 2.0);
    }

    score = 1.0 / sum;
    
    return score;

}
