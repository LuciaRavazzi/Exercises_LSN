#include <iostream>
#include <fstream>
#include "random.h"
#include <cmath>

using namespace std;

//Dati del problema.
 double sigma,mu; //Varianza e media dell'autofunzione di prova.
 double delta; //Passo.

//Variabili.
 Random rnd;
 double x;
 double xnew;
 double w;
 double sum=0.;
 int counter=0;
 double accettanza=0.;
 double pezzo1, pezzo2, pezzo3, pezzo4;

 int M=10000;
 int const N=100;
 int L=(int)M/N;
 double media[N];
 double media_prog[N];
 double media2_prog[N];
 double err_prog[N];

//Parametri termodinamici.
 double etot=0.;

//Funzioni.
 void Input(void);
 void Equilibrio(void);
 void Reset(void);
 void Move(void);
 double Density(double);
 double Wave_function(double);
 double Ekin(double);
 double Pot(double);
 void Measure(void);
 void Accumulate(void);
 void Averages(int);
 void Statistica(void);
 double error(int);
 void Print(void);


