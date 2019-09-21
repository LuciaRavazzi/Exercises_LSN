#include <iostream>
#include <fstream>
#include <cmath>

#include "costFunction.h"

using namespace std;

double costFunction :: cost(percorso per,city* c, int N_city){
 double sum=0.;
 int a,b;
 a=0; b=0;
 
 //cout << "Nella cost function " << endl;
 //for(int k=0; k<N_city; k++) cout << per.Get(k) << " ";
 //cout << endl;
//La funzione costo scelta è la somma dei moduli delle differenze tra i vettori posizione.
 for(int i=0; i<N_city-1; i++){ 
  a=per.Get(i+1);
  b=per.Get(i);
  //sum += (pow(c[a].GetX()-c[b].GetX(),2)+pow(c[a].GetY()-c[b].GetY(),2));
  sum += pow(pow(c[a].GetX()-c[b].GetX(),2)+pow(c[a].GetY()-c[b].GetY(),2),0.5);
 }

 a=per.Get(N_city-1);
 b=per.Get(0);
 //sum += (pow(c[a].GetX()-c[b].GetX(),2)+pow(c[a].GetY()-c[b].GetY(),2)); //Voglio tornare a casa.
 sum += pow(pow(c[a].GetX()-c[b].GetX(),2)+pow(c[a].GetY()-c[b].GetY(),2),0.5);

return sum;
}
