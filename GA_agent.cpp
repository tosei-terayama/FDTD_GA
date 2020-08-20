#define _USE_MATH_DEFINES
#include <cmath>
#include "GA_agent.h"
#include "perturbation.h"

/* Binaly num -> Decimal num */
int b2i(bool* b, int length){
    int d{ 0 };
    for(int i = 0; i < length; i++){
        d += b[i] * std::pow(2.0, i);
    }
    return d;
}

/* Decimel num -> Parameter */
double b2d(bool* p, int length, double min, double step){
    double param{0.0};
    param = min + b2i(p, length) * step;
    return param;
}
