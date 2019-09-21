#include <iostream>
#include <fstream>
#include "random.h"

using namespace std;
//Variabili.
 Random rnd;
 double a;
 int lanci;
//Funzioni.
 int* Reset(int,int*);
 double lorentziana(double ,double, double);
 double esponenziale(double ,double);
 double distribuzione(int, Random*,int);
 void Istogramma(double,int,int*);


double lorentziana(double T,double x, double m){ //varibile distribuita secondo una distribuzione lorentziana
 return m+T*(tan(M_PI*(x-0.5)));
}
	
double esponenziale(double lambda, double x){ //variabile distribuita secondo una distribuzione esponenziale
 return -(1./lambda)*log(1.-x);
}

int* Reset(int nBins,int* conta){
 for(int k=0; k<nBins; k++) conta[k]=0;
 return conta;
}

void Istogramma(double b,int nBins,int* conta){
 for(int k=0; k<nBins; k++){
   if(b>((double)k*a) && b<=((double)(k+1)*a)) conta[k]+=1;
 }
}

double distribuzione(int lan, Random* rnd,int quale_distribuzione){

 if(quale_distribuzione == 0){
  a=0.01;
  lanci=100000;
  int const nBins=100;
  int* conta=new int[nBins];
  double sum=0.;
  ofstream fileout;
  conta=Reset(nBins,conta);
  fileout.open("uniforme"+to_string(lan)+".dat");
  for(int i=0; i<lanci; i++)
  {
   sum=0;
   for(int k=0; k<lan; k++) sum+=rnd->Rannyu();
   Istogramma(sum/(double)lan,nBins,conta);
  }  
  for(int k=0; k<nBins; k++) fileout << k*a << " " << conta[k] << endl;
  fileout.close();
 }

 if(quale_distribuzione == 1){
  a=0.05;
  lanci=100000;
  int const nBins=100;
  int* conta=new int[nBins];
  double sum=0.;
  ofstream fileout;
  conta=Reset(nBins,conta);
  fileout.open("esponenziale"+to_string(lan)+".dat");
  for(int i=0; i<lanci; i++)
  {
   sum=0;
   for(int k=0; k<lan; k++) sum+=esponenziale(1,rnd->Rannyu());
   Istogramma(sum/(double)lan,nBins,conta);
  }  
  for(int k=0; k<nBins; k++) fileout << k*a << " " << conta[k] << endl;
  fileout.close();
 }

 if(quale_distribuzione == 2){
  a=0.1;
  double p=0;
  lanci=100000;
  int const nBins=50;
  int* conta=new int[2*nBins];
  double sum=0.;
  ofstream fileout;
  for(int i=0; i< 2*nBins; i++) conta[i]=0;
  fileout.open("lorentziana"+to_string(lan)+".dat");
  for(int i=0; i<lanci; i++)
  {
   sum=0;
   for(int k=0; k<lan; k++) sum+=lorentziana(1,rnd->Rannyu(),0);
   sum=sum/(double)lan;

    for(int k=0; k<2*nBins; k++) if(sum>((double)(k-nBins)*a) && sum<=((double)(k-nBins+1)*a)) conta[k]+=1;
   }
    for(int k=0; k<2*nBins; k++){ 
    fileout << (k-nBins)*a << " " << conta[k] << endl;
    //cout << (k-nBins)*a << " " << conta[k] << endl;
  }
   fileout.close();
  }
 
}

