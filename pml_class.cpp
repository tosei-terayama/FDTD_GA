#include <iostream>
#include <cmath>
#include "pml.h"

//PML class//
void pml::set_point_1(int py, int pz){
  p_j1 = py;
  p_k1 = pz;
}

void pml::set_point_2(int py, int pz){
  p_j2 = py;
  p_k2 = pz;
}

void pml::set_point(int vy1, int vy2, int vz1, int vz2){
  set_point_1(vy1, vz1);
  set_point_2(vy2, vz2);
}