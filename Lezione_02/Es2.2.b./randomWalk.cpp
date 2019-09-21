#include <iostream>
#include <cmath>
#include "randomWalk.h"
#include "random.h"
using namespace std;

randomWalk :: randomWalk(double _passo){
 sumx=0.;
 sumy=0.;
 sumz=0.;
 passo=_passo;
}

void randomWalk :: azzero(){
 sumx=0.;
 sumy=0.;
 sumz=0.;
}

void randomWalk :: setPasso(double a){
 passo=a;
}

double randomWalk :: move(int k){
 azzero();           
 double theta, phi;	

 for(int i=0; i<k; i++){ //faccio k passi.
   //estraggo un theta e phi casuale.
   theta=rnd.Rannyu(0.,M_PI);
   phi=rnd.Rannyu(0,2*M_PI);
   sumx+=passo*sin(theta)*cos(phi);
   sumy+=passo*sin(theta)*sin(phi);
   sumz+=passo*cos(theta);
 }
	
return (pow(sumx,2)+pow(sumy,2)+pow(sumz,2));	
}

