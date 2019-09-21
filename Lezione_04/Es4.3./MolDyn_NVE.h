/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "random.h"
#include <iostream>
using namespace std;

Random rnd;

//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp,stima_press;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut,boolean, cicli;
	

// simulation
int nstep, iprint, seed;
double delta;

//Variabili per la statistica.
 int M=10000;
 int const N=100;
 int L=(int)M/N;

 double* ekin=new double[N];
 double* epot=new double [N];
 double* etot=new double[N];
 double* T=new double[N];
 double* press=new double[N];
 double walker[5];

 double* sum_prog=new double[N];
 double* sum2_prog=new double[N];
 double* err_prog=new double[N];
//Funzioni per la simulazione.
void ScalaVel(void);
void Input(void);
void Move(void);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Clean(void);
void Equilibrazione();
double Distanza(int,int);
//Funzioni per la statistica.
void Reset(void);
void Accumulate(void);
void Averages(int);
double Error(double*,double*,int);
void Statistica(double*,string);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
