#include <iostream>
#include <fstream>
#include <cmath>
#include "randomWalk.h"
#include "random.h"
#include "main.h"

using namespace std;

int main(){

 ofstream fileout;
 randomWalk path(1.0);

 fileout.open("Statistica.dat");	

 for(int k=0; k<passi; k++){
   sum=0;
   ResetVettori();
   for(int j=0; j<N; j++){	
     sum=0;
       for(int i=0; i<L; i++) sum+=path.move(k);				
     media[j]=(double)sum/L;
   }

  Statistica(k);
  cout << k << " " << sqrt(sum_prog[N-1]) << " " << err_prog[N-1] << endl;
  fileout << k << " " << sqrt(sum_prog[N-1]) << " " << err_prog[N-1] << endl;
 }

 fileout.close();

return 0;
}


