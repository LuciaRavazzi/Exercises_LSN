#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;

 
int main (int argc, char *argv[]){

 Random rnd;	
 int esp=100;		//Numero di test del X.
 int M=100;		//Sottointervalli.
 double a=(double)1/M;
 int n=10000;		//Lanci per ogni test
 double evalue=n*((double)1/M);//Valore di apettazione in ogni sottointervallo.
 double b,X;
 double sum=0.;
 double test[M];
 int vett[M];		

 ofstream fileout;
 fileout.open("test.dat");

 for(int k=0; k<esp; k++) //ripeto 100 volte il test del X.
 {
  sum=0.;
  for(int i=0; i<M; i++) vett[i]=0;
   for(int i=0; i<n; i++) //devo lanciare n numeri random.
   {
    b=rnd.Rannyu();
    for(int j=0; j<M; j++) 
    { 
     if(b>=j*a && b<(j+1)*a) vett[j]+=1;
    }
   } 

   for(int i=0; i<M; i++) sum+=pow(vett[i]-evalue,2);
	
   test[k]=(double)sum/evalue;
   cout << k << " " << test[k] << endl;
   fileout << k+1 << " " << test[k] << endl;
  }

  fileout.close();
	
return 0;
}



