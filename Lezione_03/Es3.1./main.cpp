#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include "main.h"

using namespace std;

int main (int argc, char *argv[]){

  Random rnd;

  for(int j=0; j<N; j++){ 
    sum=0;
    for(int i=0; i<L; i++) sum+=exp(-r*T)*fmax(-(Scont(T,rnd.Gauss(0,1))-K));
      call[j]=sum/L;
  }

  Statistica();
		
return 0;
}


double error(double* sum_prog,double* su2_prog,int i){
 if(i==0){
   return 0;
 } else {
   return sqrt((su2_prog[i]-pow(sum_prog[i],2))/i);
 }
}

double Scont(double t, double c){//t è il tempo finale, c è il numero generato secondo una distribuzione normale.
 return 100.*exp(t*(r-((vol*vol)/2.))+vol*sqrt(t)*c);
}

double fmax(double a){
 if(a>0){
   return a; //Se ho profitto, ritorna il profitto.
 } else {
   return 0; //Non ho profitto.
 }	
}	


void Statistica(){
  ofstream fileout;
  fileout.open("put1.dat");	
  for(int i=0; i<N; i++){
    for(int j=0; j<i+1; j++){
      sum_prog[i]+=call[j]; //sommo le medie
      su2_prog[i]+=call[j]*call[j]; //sommo le medie al quadrato
    }

    sum_prog[i]=sum_prog[i]/(i+1); //medie progressive	
    su2_prog[i]=su2_prog[i]/(i+1); //media progressive sui quadrati	
    err_prog[i]=error(sum_prog,su2_prog,i); //incertezza statistica
    fileout << i << "   " << sum_prog[i] << "   " << err_prog[i] << endl;
    cout << i << "   " << sum_prog[i] << "   " << err_prog[i] << endl;
  }

  fileout.close();
}
