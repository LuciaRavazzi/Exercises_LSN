#include <iostream>
#include <cmath>
#include <fstream>
#include "blocchi.h"
#include "random.h"

using namespace std;

blocchi :: blocchi(int m,int n){
 M=m;
 N=n;
 L=(int)M/N;
}

blocchi :: ~blocchi(){
 delete [] media;
 delete [] sum_prog;
 delete [] su2_prog;
 delete [] err_prog;
}

void blocchi :: SetParametri(int a, int b){
 M=a;
 N=b;
 L=(int) M/N;
}

void blocchi :: MediaBlocco(double sum, int i){
 media[i]=sum/(double)L;
}


void blocchi :: ErroreBlocco(){
 for(int i=0; i<N; i++)
 {
  for(int j=0; j<i+1; j++)
  {
   sum_prog[i]+=media[j]; //sommo le medie
   su2_prog[i]+=pow(media[j],2.); //sommo le medie al quadrato
  }
  sum_prog[i]=sum_prog[i]/(i+1); //medie progressive	
  su2_prog[i]=su2_prog[i]/(i+1); //media progressive sui quadrati	
  err_prog[i]=Error(sum_prog,su2_prog,i); //incertezza statistica
 }
}


double blocchi :: Error(double* sum_prog,double* su2_prog,int i){
 if(i==0)
 {
  return 0;
 } else {
  return sqrt((su2_prog[i]-pow(sum_prog[i],2))/(double)i);
 }
}


void blocchi :: Print(){
 ofstream fileout;
 fileout.open("risultati.dat");
	
 for(int i=0; i<N; i++) fileout << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
 fileout.close();
}






