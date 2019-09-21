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
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;


int main(){ 

  Clean();
  Input(); //Inizialization.
  Equilibrazione();
	
  return 0;
}


void Clean(){
  remove("output_epotVR.dat");
  remove("output_ekinVR.dat");
  remove("output_tempVR.dat");
  remove("output_etotVR.dat");
  remove("output_epot.dat");
  remove("output_ekin.dat");
  remove("output_temp.dat");
  remove("output_etot.dat");
}

void Input(){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;
  
  cout << "	Classic Lennard-Jones fluid        " << endl;
  cout << "	Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "	Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "	The program uses Lennard-Jones units " << endl;
  
  ReadInput.open("input.dat"); //Read input
  ReadInput >> temp;
  ReadInput >> npart;
  cout << "	Number of particles = " << npart << endl;
  ReadInput >> rho;
  cout << "	Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "	Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "	Edge of the simulation box = " << box << endl;
  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> boolean;
  ReadInput >> cicli;

  cout << "	The program integrates Newton equations with the Verlet method " << endl;
  cout << "	Time step = " << delta << endl;
  cout << "	Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "	Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "	Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

void Equilibrazione(){
 int i=0; //Indica quanti cicli sono stati fatti.
  do{
    cout << "***** 	QUESTO È IL CICLO NUMERO " << i+1 << "  *****" << endl;
    
    for(int istep=1; istep <= nstep; ++istep){ 
      Move(); //Move particles with Verlet algorithm	
      if(istep % 1000 == 0) cout << "          Step " << istep << endl;
    }
  
/*Punto 1.
Dopo l'ultimo Move() del precedente ciclo, posso identificare r(t-dt) con xold mentre r(t) con x. Attraverso Verlet mi muovo fino a r(t+dt).
Dopo il prossimo Move(), r(t+dt) diventerà x mentre r(t) sarà xold; perderò informazione su r(t-dt).
LO SCOPO DI QUESTO ESERCIZIO È DI TROVARE DELLE COORDIANTE r(t-dt) CHE POSSANO RAPPRESENTARE MEGLIO LA CONDIZIONE INIZIALE. */
  Move();
	
//Con r⃗(t+dt) e r⃗ (t) devo calcolare v⃗ (t+dt/2). 
  for(int k=0; k<npart; k++){
    vx[k] = Pbc(x[k] - xold[k])/(delta);
    vy[k] = Pbc(y[k] - yold[k])/(delta);
    vz[k] = Pbc(z[k] - zold[k])/(delta);
  }	

//Con le nuove velocità calcolo la temperatura T(t+dt/2). 
  double EKin=0.;
  for (int k=0; k<npart; k++)  EKin += 0.5 * (vx[k]*vx[k] + vy[k]*vy[k] + vz[k]*vz[k]);

  cout << "Configurazione al passo r(t+dt/2)" << endl;
  cout << "Energia cinetica: " << EKin << endl;
  cout << "Energia cinetica per part: " << EKin/npart << endl;
  cout << "Temperatura: " << ((2./3.)*EKin)/(double)npart << endl;

//Devo scalare le velocità.
  ScalaVel();

//Adesso le coordinate Old definite in ScalaVel() saranno proprio le nuove r(t) poichè le x[i] sono le r(t+dt).
//Posso finalmente calcolarmi r(t+2*dt) con la funzione Move e corrisponderà alla configurazione iniziale della mia nuova simulazione.
    if (i == cicli) break;
	i++; 
	
 } while(boolean == 1);
}

void ScalaVel(){
  double sumv[3] = {0.0, 0.0, 0.0};

  for (int i=0; i<npart; ++i){
    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }

  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;

  for (int i=0; i<npart; ++i){
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }

  sumv2 /= (double)npart;

  fs = sqrt(3 * temp / sumv2);   //Sto praticamente immergendo il sistema nel bagno termico: riscalo le velocità rispetto alla temperatura target.
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = Pbc(x[i] - vx[i] * delta); //Sono le coordinate al tempo t.
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    //cout << "*** V: " << vx[i] << " " << vy[i] << " " << vz[i] << endl;

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(int loop){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
  ofstream EpotVR, EkinVR, EtotVR, TempVR;

  Epot.open("output_epot.dat",ios::app); //ios::app serve per inserire un nuovo dato in un file che è già esistente.
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i){	
	 t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	}
   
   stima_pot = v/(double)npart; //Potential energy per particle
   stima_kin = t/(double)npart; //Kinetic energy per particle
   stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
   stima_etot = (t+v)/(double)npart; //Total energy per particle
    
   Epot << stima_pot  << endl;
   Ekin << stima_kin  << endl;
   Temp << stima_temp << endl;
   Etot << stima_etot << endl;

     EpotVR.close();
     EkinVR.close();
     TempVR.close();
     EtotVR.close();

    return;
}



void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl; 
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
       WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
   //    WriteXYZ << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
