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

  Input(); //Inizialization of Molecular dynamics.
  Equilibrazione(); //Porto il sistema all'equilibrio immergendolo in un bagno alla stessa T fissata.
	
  Clean(); //I file nei quali devo scrivere le mie misure saranno puliti.
  Input2(); //Inizialization of the second part of simulation.
 
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= L; ++istep)
    {
      Move(); //with Verlet alghoritm
      Measure(); //with Monte Carlo code, not Molecular dynamics code.
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
 

	
  return 0;
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



void Clean(){
    remove("output.epot.dat");
    remove("output.pres.dat");
    remove("output.gofr.dat");
    remove("output.gave.dat");
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

	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	for (int i=0; i<npart; ++i){
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;

		xold[i] = Pbc(x[i] - vx[i] * delta);
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
	}
}




void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
 
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> boolean;
  ReadInput >> cicli;
  ReadInput >> nblk;
  ReadInput >> L;


  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();


//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

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


void Input2(){
  //Tail corrections for potential energy and pressure
  vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));

//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
 
  n_props = 2; //Number of observables

//measurement of g(r)
  igofr = 2;
  //nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;
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


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props+nbins; ++i)
       {
           glob_av[i] = 0; //Somme progressive.
           glob_av2[i] = 0;
           if(i < nbins){ stima_column[i]=0.; err_gdir[i]=0; }
       }  
   }

   for(int i=0; i<n_props+nbins; ++i)
   {
     blk_av[i] = 0; //Media in ogni blocco.
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Measure()
{
  int bin;
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
     if(dr < box/2.){
      for(int j=0; j<nbins; j++){
        if(dr > j*bin_size && dr <= (j+1)*bin_size)  walker[j+igofr] += 2; //Ogni interazione ad una certa distanza conta due volte.
      }	   				  
     }
   
     if(dr < rcut)
     {
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

// contribution to energy and virial
       v += vij;
       w += wij;
     }
    }          
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;

//Normalizzo l'istogramma.
 double deltaV=0.;
 for (int j=0; j<nbins; ++j){
   deltaV=((4.*M_PI)/(3.))*(pow((j+1)*bin_size,3)-pow(j*bin_size,3));	
   walker[j+igofr]=walker[j+igofr]/(deltaV*rho*npart);
 }
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props+nbins; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void Averages(int iblk) //Print results for current block
{
   
   double r, gdir;
   ofstream Gofr, Gave, Epot, Pres, GaveFinale, GofrFinale;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("output.epot.dat",ios::app);
    Pres.open("output.pres.dat",ios::app);
    Gofr.open("output.gofr.dat",ios::app);
    Gave.open("output.gave.dat",ios::app);
    GofrFinale.open("output.gofrFinale.dat");
    GaveFinale.open("output.gaveFinale.dat");
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

    for(int i=0; i<nbins; i++){
	stima_column[i]=blk_av[i+igofr]/(double)blk_norm; //Average of g(r).
	cout << iblk << " stima_column[i] " << stima_column[i] << endl;
        cout << iblk << " " << stima_column[i] << endl;
	Gave << iblk <<  " " << i*bin_size << " " << stima_column[i] << endl; //Salvo le medie.
    }

    for(int i=0; i<nbins; i++){
	 glob_av[i+igofr] += stima_column[i];
	 glob_av2[i+igofr] += stima_column[i]*stima_column[i];
         err_gdir[i]=Error(glob_av[i+igofr],glob_av2[i+igofr],iblk);
         if(iblk == nblk)  GaveFinale << i <<  " " << stima_column[i] << endl;
    }

//Potential energy per particle
    Epot << iblk  << " " << glob_av[iv]/(double)iblk << " " << err_pot << endl;
//Pressure
    Pres << iblk << " " << glob_av[iw]/(double)iblk << " " << err_press << endl;
//g(r)
    for(int i=0; i<nbins; i++){
     Gofr  <<  iblk <<  " " << i*bin_size << " " << glob_av[i+igofr]/(double)iblk << " " << err_gdir[i] << endl; //Salvo somme progressive e errori.
     if(iblk == nblk) GofrFinale  <<  i*bin_size <<  " " << glob_av[i+igofr]/(double)iblk << " " << err_gdir[i] << endl;
    }

    
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
    GofrFinale.close();
    GaveFinale.close();
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
