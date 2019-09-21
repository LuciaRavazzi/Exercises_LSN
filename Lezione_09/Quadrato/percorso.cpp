#include <iostream>
#include <fstream>
#include <cmath>
#include "percorso.h"
#include <vector>
#include <algorithm>
using namespace std;


void percorso :: Inizializza(int M){
 N=M;	//Numero delle città.
 tragitto=new int[N]; //Percorso.
 for(int i=0; i<N; i++) tragitto[i]=i; 
 distanza=0.;
}

void percorso :: Swap(int i,int j){
 int a;
 a=tragitto[pbc(i)];
 tragitto[pbc(i)]=tragitto[pbc(j)];
 tragitto[pbc(j)]=a;
}


int percorso :: pbc(int i){
 if(i >= N) i=i-N;
 if(i < 0) i=i+N;
 return i;
}


void percorso ::  Verify(){
//Controllo che tutte le città siano diverse. 
 for(int i=0; i<N; i++){
  for(int j=0; j<N; j++){
   if(tragitto[i] == tragitto[j] && i!=j){cout << "Errore: visito due volte la stessa città" << endl;}
  }
 }
//Controllo se ho dimenticato di visitare qualche città.
 int sum=0.;
 int sum2=0.;
 for(int i=0; i<N; i++){
  sum += i;
  sum2 += tragitto[i];
 }
 if(sum != sum2){ cout << "Non sto visitando tutte le città " << endl;} //Questo vale perchè ho valori non negativi.
}

int percorso :: Get(int i){
  return tragitto[pbc(i)];
}

void percorso :: SetDistanza(double L){
 distanza=L;
}

double percorso :: GetDistanza(){
 return distanza;
}


void percorso :: Mutazione(Random* rnd){

//Prima mutazione: scambio di due elementi.
 if(rnd->Rannyu() < 0.1){
  int i=(int)(N*rnd->Rannyu());
  int j=(int)(N*rnd->Rannyu());
  Swap(pbc(i),pbc(j));
 }
 Verify();

//Seconda mutazione: shift.
 if(rnd->Rannyu() < 0.2){
  vector <int> vec1;
  for(int i=0; i<N; i++) vec1.push_back(tragitto[i]);  

  std :: rotate(vec1.begin(),vec1.begin()+2, vec1.end());

  for(int i=0; i<N; i++) tragitto[i]=vec1[i]; 
 }
 Verify();

}

void percorso :: Print(){
 cout << endl << "Stampo percorso " << endl;
 for(int i=0; i<N; i++) cout << tragitto[i] << " ";
 cout << endl;
}

void percorso :: SetComponent(int p, int i){
  if(p >= 0 && p < N){ 
   tragitto[pbc(i)]=p;
  } else {
   cout << "Stai inserendo una città che non esiste in un percorso" << endl;
  }
}

percorso& percorso::operator= (const percorso &per){
 N = per.N;
 tragitto = new int[N];
 for(int i=0; i<N; i++) tragitto[i] = per.tragitto[i]; 

return *this;
}



