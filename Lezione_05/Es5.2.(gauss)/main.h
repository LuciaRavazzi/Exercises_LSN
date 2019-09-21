#include <iostream>
#include <fstream>
#include "random.h"
#include <cmath>
#include <ostream>
#include <iomanip>

using namespace std;

//Variabili.
 double const rb=1;
 double x,y,z;
 double xnew, ynew, znew;
 double w;
 Random rnd;
 double delta;

 double sum=0;
 int counter=0;
 ofstream fileout,fileout2;

//Variabili per il metodo a blocchi.
 int M=1000000;
 int const N=100;	
 int L=(int)M/N;
 double media[N];
 double media2[N];
 double media_prog[N];
 double media2_prog[N];
 double err_prog[N]; 	

//Dichiarazione funzioni.
double prob2p(double, double);
double error(double* ,double* ,int );
void Input(void);
void Move(void);
void Remove(void);
void Accumulate(int);
void Statistica(void);



