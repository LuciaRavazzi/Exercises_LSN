/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow

using namespace std;


//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut,boolean, cicli;
int i=0; //Indica quanti cicli sono stati fatti.	

// simulation
int nstep, iprint, seed;
double delta;

//functions
void ScalaVel(void);
void Input(void);
void Input2(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
double Force(int, int);
double Pbc(double);
void Clean(void);
void Equilibrazione(void);
void Reset(int);
void Measure(void); 
void Accumulate(void);
double Error(double sum, double sum2, int iblk);
void Averages(int);

//parameters, observables
double vtail,ptail,bin_size;
int nblk, L;
int iv, iw,n_props,igofr;
const int m_props=1000;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,err_pot,err_press;
const int nbins=100;
double err_gdir[nbins];
double stima_column[nbins];

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
