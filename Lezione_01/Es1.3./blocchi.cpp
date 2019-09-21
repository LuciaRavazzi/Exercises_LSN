#include <iostream>
#include <cmath>
#include <fstream>
#include "blocchi.h"
#include "random.h"

using namespace std;

Blocchi :: Blocchi(){} //costruttore.

Blocchi :: ~Blocchi(){//distruttore.
	delete [] media;
	delete [] media2;
	delete [] sum_prog;
	delete [] su2_prog;
	delete [] err_prog;
}

void Blocchi::setParametri(int a, int b){ //step uno: inserire i parametri.
	M=a;
	N=b;
	L=(int) M/N;
}

void Blocchi::mediaBlocco(double sum, int i){ 
	media[i]=(double)sum/L;
	media2[i]=pow(media[i],2);
}


void Blocchi::erroreBlocco(){
	for(int i=0; i<N; i++){
	for(int j=0; j<i+1; j++){
		sum_prog[i]=sum_prog[i]+media[j]; //sommo le medie
		su2_prog[i]=su2_prog[i]+media2[j]; //sommo le medie al quadrato
	}

  sum_prog[i]=sum_prog[i]/(i+1); //medie progressive	
  su2_prog[i]=su2_prog[i]/(i+1); //media progressive sui quadrati	
  err_prog[i]=error(sum_prog,su2_prog,i); //incertezza statistica
 	}
}


double Blocchi::error(double* sum_prog,double* su2_prog,int i){
	if(i==0){
		return 0;
	} else {
		return sqrt((su2_prog[i]-pow(sum_prog[i],2))/i);
	}
}


void Blocchi::fileFinale(){
	ofstream fileout;
	fileout.open("risultati.dat");
	
	for(int i=0; i<N; i++){
		fileout << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
		cout << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
	}
	fileout.close();
}






