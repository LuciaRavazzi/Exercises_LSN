/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_NVT.h"

using namespace std;

int main()
{ 

//Trovo il passo per il Metropolis tale da avere accettanza circa del 50%.
 double incr=0.005;
 delta=0.005;
 double a=0.;
 do{
   attempted=0;
   accepted=0;
   Input(); 
   for(int i=0; i<100; i++) Move();
   a=(double)accepted/attempted;
   delta+=incr;
 }while(a > 0.5);

  delta=delta-2*incr; //Il ciclo si accorge di aver trovato il valore giusto dopo aver incrementato il passo due volte.

//Inizio la simulazione.
  attempted=0;
  accepted=0;
  Input(); //Inizialization
  Equilibration();
  Clean();

  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }

  ConfFinal(); //Write final configuration


  return 0;
}

void Input(void)
{
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl; 
  double appo;
  ReadInput >> appo;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
 
  n_props = 2; //Number of observables

//measurement of g(r)
  igofr = 2;
  //nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();
  
//Evaluate potential energy and virial of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;

}


void Move(void)
{
  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;


  for(int i=0; i<npart; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(rnd.Rannyu()*npart);

  //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];

    energy_old = Boltzmann(xold,yold,zold,o);

  //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    energy_new = Boltzmann(xnew,ynew,znew,o);

  //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu())  
    {
    //Update
       x[o] = xnew;
       y[o] = ynew;
       z[o] = znew;
    
       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
  }
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
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
// cout << "Istogramma" << endl;
 double deltaV=0.;
 for (int j=0; j<nbins; ++j){
   deltaV=((4.*M_PI)/(3.))*(pow((j+1)*bin_size,3)-pow(j*bin_size,3));	
   walker[j+igofr]=walker[j+igofr]/(deltaV*rho*npart);
  // cout << walker[j+igofr] << endl;
 }
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


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props+nbins; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
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
    glob_av2[iw] +=stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

    for(int i=0; i<nbins; i++){
	stima_column[i]=blk_av[i+igofr]*pow((double)blk_norm,-1.); //Media delle componenti di g_i(r).
	cout << iblk << " stima_column[i] " << stima_column[i] << endl;
	Gave << iblk <<  " " << stima_column[i] << endl; 
    }

    for(int i=0; i<nbins; i++){
	 glob_av[i+igofr] += stima_column[i]; //Somme progressive.
	 //cout << "glob_av[i+igofr]/(double)iblk " <<  glob_av[i+igofr]/(double)iblk << " glob_av[i+igofr] " << glob_av[i+igofr] << endl;
         glob_av2[i+igofr] += stima_column[i]*stima_column[i];
         //cout << iblk << " " << bin_size*i << " glob_av[i+igofr] " << glob_av[i+igofr] << " glob_av2[i+igofr] " << glob_av2[i+igofr] << endl;
         err_gdir[i]=Error(glob_av[i+igofr],glob_av2[i+igofr],iblk);
	 //cout << bin_size*i << " err_gdir[i] " << err_gdir[i] << endl;
         if(iblk == nblk)  GaveFinale << i <<  " " << stima_column[i] << endl;
    }

//Potential energy per particle
    Epot << iblk << " "  << glob_av[iv]/(double)iblk << " " << err_pot << endl;
//Pressure
    Pres << iblk <<  " " <<  glob_av[iw]/(double)iblk << " " << err_press << endl;
//g(r)
    for(int i=0; i<nbins; i++){
     cout  <<  iblk << " " << bin_size*i <<  " " << glob_av[i+igofr]/(double)iblk << " " << err_gdir[i] << endl; //Salvo somme progressive e errori.
     Gofr  <<  iblk <<  " " << glob_av[i+igofr]/(double)iblk << " " << err_gdir[i] << endl; //Salvo somme progressive e errori.
     if(iblk == nblk) GofrFinale  <<  bin_size*i <<  " " << glob_av[i+igofr]/(double)iblk << " " << err_gdir[i] << endl;
    }

    
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
    GofrFinale.close();
    GaveFinale.close();
}



void ConfFinal(void)
{
  ofstream WriteConf, WriteSeed;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    //cout << "sum " << sum << endl;
    //cout << "sum2 " << sum2 << endl;
    //cout << "sum2/(double)iblk " << sum2/(double)iblk << endl;
    //cout << "pow(sum/(double)iblk,2) " << pow(sum/(double)iblk,2) << endl;
    cout << "Err " << pow(((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk),0.5) << endl;
    return pow(((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk),0.5);
}

void Equilibration(){
 ifstream filein;
 filein.open("TempiDiCorrelazione");

 if (filein.is_open()){
   filein >> tcPot; //Tempo di correlazione del potenziale.
   filein >> tcPress; //Tempo di correlazione della pressione.
   filein.close();
   int tc;

   cout << "Equilibro il sistema" << endl;

   if(tcPot > tcPress){
    tc=tcPot;
   } else {
    tc=tcPress;
   }

  for(int i=0; i<10*tc; i++){ 
   Move(); //Mi muovo per un multiplo del tempo di correlazione.
   Measure();
   //cout << "U/N: " << walker[iv]/(double)npart + vtail << endl;
   //cout << "Pressure: " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
   }
  } else {
   cerr << "************ PROBLEMA APERTURA FILE TEMPI DI CORRELAZIONE ***************" << endl;
 }
}

void Clean(){
 remove("output.epot.dat");
 remove("output.pres.dat");
 remove("output.gofr.dat");
 remove("output.gave.dat");
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
