#include <iostream>
#include <cmath>
#include <fstream>
#include "ago.h"
#include "random.h"

using namespace std;

void Ago :: setParametri(int a, int b, int p){ 
 L=a;
 D=b;
 lanci=p;
}

Ago :: Ago(){
 rnd=new Random();
}

Ago :: ~Ago(){}


double Ago :: pigreco(){
 double x,theta;
 double counter=0;
		
 for(int i=0; i<lanci; i++){
  //Lancio l'ago.
  x=rnd->Rannyu(0,D);
  theta=rnd->Uniforme();
  //Ha intersecato una riga del pavimento?
  if( x > D/2.){
   if( x + (L/2.)*abs(cos(theta)) >= D) counter++;
   } else {
     if( x - (L/2.)*abs(cos(theta)) <= 0) counter++;
   }	
  }			


return (double)(2.*L*lanci)/(double)((double)counter*D);
}

