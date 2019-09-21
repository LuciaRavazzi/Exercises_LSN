#include <iostream>
#include <fstream>
#include <cmath>
#include "city.h"
#include "percorso.h"
#include "costFunction.h"
#include <algorithm>
#include "crossover.h"

using namespace std;

bool Metropolis(percorso,percorso,Random*,double);

int main(){

  Random rnd;	
  ofstream fileout;
  fileout.open("Distanza-beta.dat");
//Creo le posizioni delle città.
  int N_city=30;
  double raggio=1.0;

  city* ct=new city[N_city];
  for(int i=0; i<N_city; i++) ct[i].Set(rnd.Rannyu(0,2*M_PI), raggio); 
 
//Metto in ordine.
 int i, j, min;
 double temp;

 for (i = 0; i < N_city - 1; i++) {
  min = i;
  for (j = i + 1; j < N_city; j++)
   if(ct[j].GetAngle() < ct[min].GetAngle()) min = j;
  
  temp=ct[i].GetAngle();
  ct[i].Set(ct[min].GetAngle(),raggio);
  ct[min].Set(temp,raggio);
 }

 for(int i=0; i<N_city; i++) cout << i << " " << ct[i].GetX() << " " << ct[i].GetY() << " " << ct[i].GetAngle() <<  endl;
  

//Simulated annealing. Immergo il sistema in un bagno termico a temperatura 1/beta. Ogni mossa consiste nella generazione di un nuovo percorso il quale ha subito delle mutazioni. L'algoritmo di Metropolis regola se la mossa può essere accettata o meno.

  double betaMax=80;
  double betaMin=0.15;
  int steps=1000;
  double beta_step=(betaMax-betaMin)/(double)steps;
  double beta=betaMin;

  //Genero il primo percorso.
  costFunction L;

  percorso old;
  old.Inizializza(N_city); 
  int N_scambi=30;
 
  for(int j=0; j<N_scambi; j++) old.Swap((int)(rnd.Rannyu()*N_city),(int)(rnd.Rannyu()*N_city));
  old.SetDistanza(L.cost(old,ct,N_city)); 
  old.Print();
  cout << "Distanza primo percorso " << old.GetDistanza() << endl;
  old.Verify();


  //Diminuisco dolcemente la temperatura.
  do{
    cout << "					Temperatura: " << 1./beta << endl;
    cout << "					Beta: " << beta  << endl;
    //Alla temperatura fissata deve muovermi di steps.
   for(int i=0; i<=steps; i++){
      //Genero una mossa casuale dal vecchio passo, creo il nuovo passo che potrebbe essere accettato.
      percorso New;
      New.Inizializza(N_city); 
      New=old;
      New.SetDistanza(L.cost(old,ct,N_city));
      //Mutazione: permutazione di alcuni elementi.
      New.MutazionePermutazione(&rnd);
      New.SetDistanza(L.cost(New,ct,N_city));
      New.Verify();
      //cout << "Prima" << endl;
      //New.Print();
     //cout << i << " Distanza percorso " << New.GetDistanza() << endl;
      //New.Print();
     
      if(Metropolis(New,old,&rnd,beta) == true){ 
        old=New; //La nuova mossa accettata diventa quella vecchia.
        old.SetDistanza(L.cost(New,ct,N_city));
        //cout << "Accettato" << endl;
        // cout << i << " Distanza percorso " << New.GetDistanza() << endl;
        //New.Print();
        if(i == steps){ fileout << beta << " " << New.GetDistanza() << endl; cout << beta << " " << New.GetDistanza() << endl; }
      }

      //Se non è stata accettata, rimango su quella vecchia.
      New=old; //La nuova mossa è sempre quella vecchia che verrà mutata.
      New.SetDistanza(L.cost(old,ct,N_city));
     //Mutazione: Shift.
      New.MutazioneShift(&rnd);
      New.SetDistanza(L.cost(New,ct,N_city)); 
      New.Verify();
      if(Metropolis(New,old,&rnd,beta) == true){ 
       old=New; //La nuova mossa accettata diventa quella vecchia.
       old.SetDistanza(L.cost(New,ct,N_city));
       //cout << i << " Distanza percorso " << New.GetDistanza() << endl;
       if(i == steps){ fileout << beta << " " << New.GetDistanza() << endl; cout << beta << " " << New.GetDistanza() << endl; }
       //New.Print();
      }

      New=old;
      New.SetDistanza(L.cost(old,ct,N_city));
   }
    beta+=beta_step; 
 }while(beta <= betaMax);
 
 fileout.close();
 
 ofstream out;
 out.open("distanza.dat");
 for(int i=0; i<N_city; i++)
  out << ct[old.Get(i)].GetX() << " " << ct[old.Get(i)].GetY() << endl;

 out << ct[old.Get(0)].GetX() << " " << ct[old.Get(0)].GetY() << endl;
 out.close();

return 0;
}


bool Metropolis(percorso New,percorso old, Random* rnd,double beta){
  //Implemento l'algoritmo di Metropolis per il passo proposto.
  double w=exp(-beta*New.GetDistanza())/exp(-beta*old.GetDistanza());
  bool boolean=false;

  if(w >= 1){
    boolean=true;
  } else {
    if(rnd->Rannyu() <= w) boolean=true;
  }
   
return boolean;
}



