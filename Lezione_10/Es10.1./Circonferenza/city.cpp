#include <iostream>
#include <fstream>
#include <cmath>

#include "city.h"

using namespace std;

void city :: Set(double t, double r){
 theta=t;
 x=r*cos(t);
 y=r*sin(t);
}

double city :: GetX(){
 return x;
}

double city :: GetY(){
 return y;
}

double city :: GetAngle(){
 return theta;
}



