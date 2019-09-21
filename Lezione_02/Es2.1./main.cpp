#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include "main.h"
#include <iomanip>

using namespace std;


double functionI(double);


int main (){

//Primo punto

 for(int i=0; i<N; i++){ //ho calcolato l'integrale con il metodo della media
   sum=0;
   for(int j=0; j<L; j++) sum+=function(rnd.Rannyu()); 
   inte[i]=(double)sum/L;
 }
			
 PrintStatistica(0);
 ResetVettori();
//Secondo punto
//Importance sampling: devo trovare una buona funzione di probabilità
//Sembra che la funzione di probabilità più adatta a tale scopo sia proprio una retta.
//La retta che approssima meglio la funzione data è proprio -2x+2
//Devo costruire una funzione che mi restituisca x distribuite in tale modo
	
 for(int i=0; i<N; i++){ //ho calcolato l'integrale con il metodo della media
   sum=0;
   for(int j=0; j<L; j++) sum+=functionI(1+sqrt(1-rnd.Rannyu()));//variabile distribuita come una retta
   inte[i]=(double)sum/L;
 }
	
 PrintStatistica(1);	




return 0;
}



