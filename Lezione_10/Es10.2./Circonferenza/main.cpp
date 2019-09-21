#include <iostream>
#include <fstream>
#include <cmath>
#include "city.h"
#include "percorso.h"
#include "costFunction.h"
#include <algorithm>
#include "crossover.h"
#include "mpi.h"
#include "random.h"

using namespace std;

bool Metropolis(percorso,percorso,Random*,double);
int findSmallestElement(double arr[], int n);

int main(int argc, char* argv[]){

  //Calcolo parallelo.
  MPI::Init(argc,argv);
  int size=MPI::COMM_WORLD.Get_size();
  int rank=MPI::COMM_WORLD.Get_rank();
  cout << "*****	Sono	*****  " << rank << endl;

  //Ogni rank ha un generatore diverso. Altrimenti le simulazioni sarebbero tutte uguali.
  Random rnd(rank);	

  //Creo le posizioni delle città: le creo in rank=0 e poi le condivido.
  int N_city=30;
  double raggio=1.0;

  city* ct=new city[N_city];

  if(rank == 0){
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
  }

  //Condivido con gli altri nodi le città.
  for(int i=0; i<N_city; i++){
    MPI_Bcast(&(ct[i].x),1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);
    MPI_Bcast(&(ct[i].y),1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);
    MPI_Bcast(&(ct[i].theta),1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);
  }
 
  //Simulated annealing. Ogni nodo eseguirà tale tecnica.

  double betaMax=100;
  double betaMin=0.01;
  int steps=1000;
  double beta_step=(betaMax-betaMin)/(double)steps;
  double beta=betaMin;

  //Genero il primo percorso (diverso per ogni nodo).
  costFunction L;

  percorso old;
  old.Inizializza(N_city); 
  int N_scambi=30;
 
  for(int j=0; j<N_scambi; j++) old.Swap((int)(rnd.Rannyu()*N_city),(int)(rnd.Rannyu()*N_city));
  old.SetDistanza(L.cost(old,ct,N_city)); 
  old.Print();
  //cout << "Distanza primo percorso " << old.GetDistanza() << endl;
  old.Verify();

  ofstream out;
  out.open("Distanze.dat");

  //Diminuisco dolcemente la temperatura.
  do{
   // cout << "					Temperatura: " << 1./beta << endl;
   // cout << "					Beta: " << beta  << endl;
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
        //if(i == steps){ fileout << beta << " " << New.GetDistanza() << endl; cout << beta << " " << New.GetDistanza() << endl; }
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
       //if(i == steps){ fileout << beta << " " << New.GetDistanza() << endl; cout << beta << " " << New.GetDistanza() << endl; }
       //New.Print();
      }

      New=old;
      New.SetDistanza(L.cost(old,ct,N_city)); 
   }
    if(rank == 0) out << beta << " " << old.distanza << endl;

    beta+=beta_step; 
 }while(beta <= betaMax);

 out.close();

 //Aspetto che tutti abbiano finito.
  MPI_Barrier(MPI::COMM_WORLD);

 //Creo un vettore per salvare le distanze finali.
  double vett[size];

 //Raccolgo tutti i valori nel rank=0 per poi confrontarli.
 MPI_Gather(&(old.distanza),1,MPI_DOUBLE_PRECISION,vett,1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);
 
 int min;
 if(rank == 0){ 
   cout << "risultati" << endl; 
   for(int i=0; i<size; i++) cout << i << " " << vett[i] << endl;
   //Trovo il minimo della distanza nel vettore.
   min=findSmallestElement(vett,size);
   cout << "L'elemento più piccolo " << vett[min] << " del rank " << min << endl;
 }

 //Devo dire a tutti gli altri chi è il più piccolo perchè il calcolo è stato fatto solo nel rank=0.
 MPI_Bcast(&min,1,MPI_INTEGER,0,MPI::COMM_WORLD);

 //Una volta trovato l'elemento più piccolo, devo recuperare i dati di quel rank.
 if(rank == min){ //È fonadamentale che tutti sappiano chi è min.
   cout << "Adesso stampo il miglior percorso " << endl;
   ofstream fileout;
   old.Print();
   fileout.open("BestPath.dat");
   for(int i=0; i<N_city; i++){
     fileout << ct[old.Get(i)].x << " " << ct[old.Get(i)].y << endl;
     cout << ct[old.Get(i)].x << " " << ct[old.Get(i)].y << endl;
   }
   fileout << ct[old.Get(0)].x << " " << ct[old.Get(0)].y << endl;
   fileout.close();
 }  

 MPI::Finalize();

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

int findSmallestElement(double arr[], int n){
 double temp=arr[0];
 int k=0;

 for(int i=0; i<n; i++){
  if(temp >= arr[i]){ 
   temp=arr[i];
   k=i;
  }
 }
  return k;
}


