#include <iostream>
#include <cmath>
#include "randomWalk.h"
#include <fstream>

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
 double c;	
 for(int i=0; i<k; i++){ //faccio k passi.
   c=rnd.Rannyu(0,3);
   if(c>=0 && c<1){
     if(rnd.Rannyu() >= 0.5){
       sumx+=passo;
     } else {
	sumx=sumx-passo;
     }
   }

   if(c>=1 && c<2){
     if(rnd.Rannyu() >= 0.5){
       sumy+=passo;
     } else {
	sumy=sumy-passo;
     }
   }

   if(c>=2 && c<=3){
     if(rnd.Rannyu() >= 0.5){
       sumz+=passo;
     } else {
       sumz=sumz-passo;
     }
    }
  }

	
return (pow(sumx,2)+pow(sumy,2)+pow(sumz,2));	
}

