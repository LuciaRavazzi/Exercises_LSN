#include <iostream>
#include <fstream>
#include <cmath>

#include "city.h"

using namespace std;

void city :: Set(Random* rnd){
 x=rnd->Rannyu(-1.0,1.0);
 y=rnd->Rannyu(-1.0,1.0);
}

double city :: GetX(){
 return x;
}

double city :: GetY(){
 return y;
}

double city :: SetX(double x_p){
 x=x_p;
}

double city :: SetY(double y_p){
 y=y_p;
}




