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
#include <string>
#include "random.h"
#include "blocchi.h"
#include <cmath>

using namespace std;


int main (int argc, char *argv[]){

 Random rnd;
 double sum=0.;
 int M=10000000;
 int N=100;
 blocchi bl(M,N);

 for(int i=0; i<bl.N; i++)
 {
  sum=0.;
  for(int j=0; j<bl.L; j++) sum+=rnd.Rannyu();
  bl.MediaBlocco(sum,i);
 }
 
 bl.ErroreBlocco();
	
 for(int i=0; i<bl.N; i++) cout << i << " " << bl.sum_prog[i] << " " << bl.err_prog[i] << endl;
 bl.Print();
  
return 0;
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

