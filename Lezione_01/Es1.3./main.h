#include <iostream>
#include <fstream>

using namespace std;



 double L=2.;//lunghezza dell'ago
 double D=4.;//distanza tra due assi del pavimento.
 double lanci=10000; //mi serve per fare una stima di pigreco.
 double sum=0;

 int M=10000; //numero totale di lanci.
 int N=100; //numero totali di blocchi.
 double S=(double)M/N;//lanci all'interno di ogni blocco.
 
