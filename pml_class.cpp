#define _USE_MATH_DEFINES
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

void pml::set_point(int v_y1, int v_y2, int v_z1, int v_z2){
  set_point_1(v_y1, v_z1);
  set_point_2(v_y2, v_z2);
}