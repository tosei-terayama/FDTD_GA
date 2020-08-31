#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>

int main(void)
{
    std::ifstream ifs1;
    std::ifstream ifs2;
    ifs1.open("./target.dat");
    ifs2.open("./magnitude.dat");
    double mag_1[880];
    double mag_2[880];
    double score{0.0};

    int buff;

    for(int i = 0; i < 880;i++){
        ifs1 >> buff;
        ifs1 >> mag_1[i];
        ifs2 >> buff;
        ifs2 >> mag_2[i];
    }
    ifs1.close();
    ifs2.close();

    for(int i = 0; i < 880; i++){
        score += std::pow(mag_1[i] - mag_2[i], 2.0);
    }

    score = 1.0/score;

    std::cout << "score : " << score << std::endl;

    return 0;
}