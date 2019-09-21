#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"
#include "blocchi.h"
#include "ago.h"
#include "main.h"
using namespace std;

int main(){

 Blocchi bl;
 Ago needle;

 bl.setParametri(M,N);
 needle.setParametri(L,D,lanci);

 for(int j=0; j<N; j++){ //su tutti i blocchi.
   sum=0;
   for(int i=0; i<S; i++) sum+=needle.pigreco();
   bl.mediaBlocco(sum,j);
 }
	
 bl.erroreBlocco();
 bl.fileFinale();


return 0;
}
