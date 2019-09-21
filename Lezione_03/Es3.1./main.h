#include <iostream>
#include <fstream>

using namespace std;

 double sum=0;
 double c=0;	

 const double siniz=100.; //Prezzo iniziale.
 const double T=1.;       //Tempo dopo la quale scade l'opzione.
 const double K=100.;     //Prezzo concordato.
 const double r=0.1;      //Tasso d'interesse.
 const double vol=0.25;   //Volatilità.

 int M=100000;//numero lanci totali
 int N=100;//numero blocchi

 double L=(double) M/N;
 double* call=new double[N];
 double* sum_prog=new double[N];
 double* su2_prog=new double[N];
 double* err_prog=new double[N];

//Dichiarazione funzioni.
double error(double* ,double* ,int );
double Scont(double , double );
void Statistica();
double fmax(double);
double Scont(double, double);



