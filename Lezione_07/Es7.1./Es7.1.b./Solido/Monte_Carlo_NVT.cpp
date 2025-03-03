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
   cout << "Passo " << delta << " acc " << a << endl;
   delta+=incr;
 }while(a > 0.5);

  delta=delta-2*incr; //Il ciclo si accorge di aver trovato il valore giusto dopo aver incrementato il passo due volte.

//Parto con la simulazione.
  Clean();
  int L=10;
do{ 
 Input(); 			//Inizialization
 //Equilibration(); 		//Porto il sistema all'equilibrio.
//Dopo che ho portato all'equilibrio il sistema, posso iniziare a misurare,
  nblk=(int)pow(10,5)/L;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= L; ++istep)
    {
      cout << nblk << " " << iblk << " " << istep << endl;
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk,L);   //Print results for current block
  }
  L=L+100;
} while( L < 5*pow(10,3)); 


  return 0;
}



void Input(void)
{
  ifstream ReadInput,ReadConf;

//  cout << "Classic Lennard-Jones fluid        " << endl;
//  cout << "Monte Carlo simulation             " << endl << endl;
//  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
//  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
//  cout << "The program uses Lennard-Jones units " << endl;

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
//  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
 // cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
//  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
//  cout << "Volume of the simulation box = " << vol << endl;
//  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
//  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
//  cout << "Tail correction for the potential energy = " << vtail << endl;
//  cout << "Tail correction for the virial           = " << ptail << endl; 
  double appo;
  ReadInput >> appo;

  ReadInput >> nblk;

  ReadInput >> nstep;

//  cout << "The program perform Metropolis moves with uniform translations" << endl;
//  cout << "Moves parameter = " << delta << endl;
//  cout << "Number of blocks = " << nblk << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
 
  n_props = 2; //Number of observables

//measurement of g(r)
  igofr = 2;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
//  cout << "Read initial configuration from file config.0 " << endl << endl;
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
//  cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
//  cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
//  cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
//Read time of correlation from files.
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
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {

// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)

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
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
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
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Equilibration(){
 ifstream filein;
 filein.open("TempiDiCorrelazione");

 if (filein.is_open()){
   filein >> tcPot; //Tempo di correlazione del potenziale.
   filein >> tcPress; //Tempo di correlazione della pressione.
   cout << tcPot << " " << tcPress << endl;
   filein.close();
   int tc;

  // cout << "Equilibro il sistema" << endl;

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
  remove("output.epotLAST.0");
  remove("output.presLAST.0");
}

void Averages(int iblk,int L) //Print results for current block
{
    
   double r, gdir;
   ofstream  EpotLAST, PresLAST;
   const int wd=12;
    
//    cout << "Block number " << iblk << endl;
//    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    EpotLAST.open("output.epotLAST.0",ios::app);
    PresLAST.open("output.presLAST.0",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

//Potential energy per particle
    if(iblk == nblk){ EpotLAST << L <<  " " << err_pot << endl; cout << L <<  " " << err_pot << endl; }
//Pressure
    if(iblk == nblk){  PresLAST << L <<  " " << err_press << endl; cout << L <<  " " << err_press << endl; }

//     cout << "L: " <<  L <<  "err_pot: " << err_pot << endl;
//         cout << "L: " <<  L <<  "err_press: " << err_press << endl;


//    cout << "----------------------------" << endl << endl;

    EpotLAST.close();
    PresLAST.close();
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
