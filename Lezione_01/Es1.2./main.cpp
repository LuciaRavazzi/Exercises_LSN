#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "main.h"
#include "random.h"

using namespace std;


int main (int argc, char *argv[]){
//Distribuzione uniforme.

 int myArray[] {1, 2, 10, 100};
 for(int i=0; i<4; i++){ 
  distribuzione(myArray[i],&rnd,0);
  distribuzione(myArray[i],&rnd,1);
  distribuzione(myArray[i],&rnd,2);
 }	

return 0;
}




