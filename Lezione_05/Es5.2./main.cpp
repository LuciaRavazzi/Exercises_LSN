#include <iostream>
#include "random.h"
#include <cmath>
#include "main.h"

using namespace std;

int main(){
  Remove();
  Input(); //Condizioni iniziali.

  for(int j=0; j<N; j++){	
    sum=0;
    for(int i=0; i<L; i++) Move();
    Accumulate(j);
  }	

  Statistica();
  cout << "Accettanza " << (double)counter/M << endl;

return 0;
}

double prob2p(double x, double y, double z){
	double r=sqrt(x*x+y*y+z*z);
	return (pow(rb,-5)/(M_PI*32.0))*pow(z,2)*exp(-r/rb);
}

double error(double* sum_prog,double* su2_prog,int i){
  if(i==0){
    return 0;
  } else {
    return sqrt((su2_prog[i]-pow(sum_prog[i],2))/i);
  }
}

void Input(){
 x=1.;
 y=1.;
 z=1.;

 delta=3.0;
}

void Move(){
 fileout.open("orbitale2p.dat", ios::app); 
 xnew=x+delta*(2*rnd.Rannyu()-1);
 ynew=y+delta*(2*rnd.Rannyu()-1);
 znew=z+delta*(2*rnd.Rannyu()-1);

 w=(double)prob2p(xnew,ynew,znew)/prob2p(x,y,z);
 if(w >= 1)
 {
   fileout << xnew << " " << ynew << " " << znew << endl;
   x=xnew;
   y=ynew;
   z=znew;
   sum+=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
   counter++;
 } else {
   if(rnd.Rannyu() <= w)
   { 
     fileout << xnew << " " << ynew << " " << znew << endl;		
     x=xnew;
     y=ynew;
     z=znew;
     sum+=sqrt(pow(x,2)+pow(y,2)+pow(z,2));	
     counter++;
    } else {//non mi muovo.
     sum+=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    }
 }
 fileout.close();
} 

void Remove(){
 remove("orbitale2p.dat");
}

void Accumulate(int j){
 media[j]=sum/L;	
 media2[j]=pow(media[j],2);
}

void Statistica(){
 fileout2.open("Statistica");
 for(int i=0; i<N; i++){
   for(int j=0; j<i+1; j++){
     media_prog[i]+=media[j]; //sommo le medie
     media2_prog[i]+=media2[j]; //sommo le medie al quadrato
    }

  media_prog[i]=media_prog[i]/(i+1); //medie progressive	
  media2_prog[i]=media2_prog[i]/(i+1); //media progressive sui quadrati	
  err_prog[i]=error(media_prog,media2_prog,i); //incertezza statistica
  fileout2 << i << " " << media_prog[i] << " " << err_prog[i] << endl;
  cout << i << " " << media_prog[i] << " " << err_prog[i] << endl;
  }
 
 fileout2.close();
}




