#include <iostream>

using namespace std;

//Variabili.
 Random rnd;

 int M=100000; //lanci totali
 int N=100; //numero blocchi
 int L=M/N; //quanti numeri casuali in un blocco
 double sum=0;

 double* inte=new double[N];
 double* inte_prog=new double[N]; //integrale in ogni blocco
 double* inte2_prog=new double[N];
 double* err_prog=new double[N];


//Dichiarazione funzioni.
 double error(double* ,double* ,int );
 double function(double);
 void PrintStatistica(int);
 void ResetVettori();


//Definizione funzioni.
double error(double* sum_prog,double* su2_prog,int i){
  if(i==0){
    return 0;
  } else {
    return sqrt((su2_prog[i]-pow(sum_prog[i],2))/i);
  }
}

double function(double x){
 double k=M_PI/2.;
 return (k)*cos(k*x);
}

double functionI(double x){
  double k=M_PI/2.;
  return (k*cos(k*x))/(-2.*x+2.);
}

void PrintStatistica(int j){
 ofstream fileout;

 for(int i=0; i<N; i++){
  for(int j=0; j<i+1; j++){
   inte_prog[i]+=inte[j]; //sommo le medie
   inte2_prog[i]+=inte[j]*inte[j]; //sommo le medie al quadrato
  }
  inte_prog[i]=inte_prog[i]/(i+1); //medie progressive	
  inte2_prog[i]=inte2_prog[i]/(i+1); //media progressive sui quadrati	
  err_prog[i]=error(inte_prog,inte2_prog,i); //incertezza statistica
 }	

 fileout.open("Statistica" + to_string(j) + ".dat");
 for(int i=0; i<N; i++) fileout << i << " " << inte_prog[i] << " " << err_prog[i] << endl;   
 fileout.close();
}

void ResetVettori(){
  for(int i=0; i<N; i++){
    inte[i]=0;
    inte_prog[i]=0;
    inte2_prog[i]=0;
    err_prog[i]=0;
  }
}

