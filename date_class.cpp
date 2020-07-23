#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "date.h"
#include "geocoordinate.h"

void date::set_y(int year){
    Year = year;
}

void date::set_m(int month){
    Month = month;
}

void date::set_d(int day){
    Day = day;
}

void date::set_h(double hour){
    Hour = hour;
}

void date::set_ymd(int year, int month, int day){
    set_y(year);
    set_m(month);
    set_d(day);

    MD = month*100 + day;

}