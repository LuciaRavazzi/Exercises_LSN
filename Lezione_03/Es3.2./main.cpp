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
    for(int k=0; k<L; k++){
      t=0;
      p=siniz;
      for(int i=0; i<appo; i++){//trovo il prezzo al tempo T discretizzando il tempo.
        s=Sdisc(p,t,rnd.Gauss(0,1),passo);
	t=t+passo;
	p=s;//memorizzo il valore precedente in p.
      }
      sum+=exp(-r*T)*fmax(-s+K);//call	
    }
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

double Sdisc(double sprima,double t, double c,double passo){
  return sprima*exp((r-((vol*vol)/2))*passo+vol*c*sqrt(passo));
}

double fmax(double a){
  if(a>0){
    return a; //se ho profitto, ritorna il profitto.
  } else {
    return 0; //se non ho profitto, ritorna zero.
  }	
}	


void Statistica(){
  ofstream fileout;
  fileout.open("put2.dat");
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

