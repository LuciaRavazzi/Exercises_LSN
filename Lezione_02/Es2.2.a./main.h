#include <fstream>
#include <iostream>

using namespace std;

 int M=10000;	
 int N=100;
 int L=M/N;
 int passi=100;
 double sum=0;

 double* media=new double[N];
 double* sum_prog=new double[N];
 double* su2_prog=new double[N];
 double* err_prog=new double[N];


 double error(double*, double*,int);
 void ResetVettori();
 void Statistica(int);


void ResetVettori(){ 
  for(int l=0; l<N; l++){
    media[l]=0;
    sum_prog[l]=0;
    su2_prog[l]=0;
    err_prog[l]=0;
  } 
}

double error(double* sum_prog,double* su2_prog,int i){
 if(i==0){
  return 0;
 } else {
  return sqrt((su2_prog[i]-pow(sum_prog[i],2))/i);
 }
}

void Statistica(int k){	
 for(int i=0; i<N; i++){
    for(int j=0; j<i+1; j++){
      sum_prog[i]+=media[j]; //sommo le medie
      su2_prog[i]+=media[j]*media[j]; //sommo le medie al quadrato
    }	

    sum_prog[i]=sum_prog[i]/(i+1); //medie progressive	
    su2_prog[i]=su2_prog[i]/(i+1); //media progressive sui quadrati
    if(k == 0){	
       err_prog[i]=0.;
    } else {
      err_prog[i]=error(sum_prog,su2_prog,i); //incertezza statistica
      if(i >= 1) err_prog[i]=(0.5*err_prog[i]*sqrt(sum_prog[i]))/sum_prog[i];
    }
  // cout << i << " " << sqrt(sum_prog[i]) << " " << err_prog[i] << endl;
 }

 
}


