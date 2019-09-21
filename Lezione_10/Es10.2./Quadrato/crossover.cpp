#include <iostream>
#include <fstream>
#include <cmath>
#include "crossover.h"
#include "percorso.h"

using namespace std; 

void crossover :: Input(percorso m, percorso d, int N, int a){
 N_city=N;
 quanti_scambi=a; 
//Costruisco la mamma.
 Mom=m;
//Costruisco il papà.
 Dad=d;
//Costruisco il primo figlio, è uguale alla mamma.
 FirstSon.Inizializza(N);
 FirstSon=Mom;
//Costruisco il secondo figlio, deve essere uguale al papà.
 SecondSon.Inizializza(N);
 SecondSon=Dad;
}

void crossover :: Accoppiamento(Random* rnd){

 int h=0;
 if(rnd->Rannyu() > 0.5){
//Generazione SecondSon.
 for(int i=0; i<N_city; i++){
  bool boolean=false;
  for(int j=0; j<N_city-quanti_scambi; j++){
   if(FirstSon.Get(i) == SecondSon.Get(j)) boolean=true;
  }
  if(boolean == false){ SecondSon.SetComponent(FirstSon.Get(i),N_city-quanti_scambi+h); h++;}
 }

 //Dovrà avere la sua distanza, non quella della madre! Cambio nel main.

//Generazione secondo figlio.
 h=0;
 for(int i=0; i<N_city; i++){
  bool boolean=false;
  for(int j=0; j<N_city-quanti_scambi; j++){
   if(FirstSon.Get(j) == SecondSon.Get(i)) boolean=true;
  }
  if(boolean == false){ FirstSon.SetComponent(SecondSon.Get(i),N_city-quanti_scambi+h); h++;}
 }
 }
 
  //Dovrà avere la sua distanza, non quella del papà! Cambio nel main.
}

percorso crossover :: GetFirtsSon(){
 return FirstSon;
}
 
percorso crossover :: GetSecondSon(){
 return SecondSon;
}

double crossover :: GetDistanzaFirtsSon(){
 return FirstSon.GetDistanza();
}
 
double crossover :: GetDistanzaSecondSon(){
 return SecondSon.GetDistanza();
}

void crossover :: Swap(int i){
 int a;

 a=FirstSon.Get(i);
 FirstSon.SetComponent(SecondSon.Get(i),i);
 SecondSon.SetComponent(a,i);
}

