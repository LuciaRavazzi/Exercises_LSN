#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "city.h"
#include "percorso.h"
#include "costFunction.h"
#include <algorithm>
#include "crossover.h"

using namespace std;

percorso* Ordina(percorso*,int,int);

int main(){

 Random rnd;	

//Creo le posizioni delle città.
 int N_city=30;
 double raggio=1.0;

 cout << "Città " << endl;
 city* ct=new city[N_city];
 for(int i=0; i<N_city; i++){ 
   ct[i].Set(&rnd);
   cout << ct[i].GetX() << " " << ct[i].GetY() << endl;
 }


//Metto in ordine in base alla distanza dall'origine.
 int i, j, min;
 double tempx, tempy;

 for (i = 0; i < N_city - 1; i++) {
  min = i;
  for (j = i + 1; j < N_city; j++){
   double a=sqrt(ct[j].GetX()*ct[j].GetX()+ct[j].GetY()*ct[j].GetY());
   double b=sqrt(ct[min].GetX()*ct[min].GetX()+ct[min].GetY()*ct[min].GetY());
   if(a < b) min = j;
  }
  tempx=ct[i].GetX();
  tempy=ct[i].GetY();
  ct[i].SetX(ct[min].GetX());
  ct[i].SetY(ct[min].GetY());
  ct[min].SetX(tempx);
  ct[min].SetY(tempy);
 }
  cout << "Città ordinate" << endl;
  for(int i=0; i<N_city; i++) cout << sqrt(ct[i].GetX()*ct[i].GetX()+ct[i].GetY()*ct[i].GetY()) << " " << ct[i].GetX() << " " << ct[i].GetY() << endl;

//Creo la popolazione.
 int N_per=900;
 int N_scambi=30;
 
 percorso* popolazione=new percorso[N_per];
 for(int i=0; i<N_per; i++){
  popolazione[i].Inizializza(N_city);
  for(int j=0; j<N_scambi; j++) popolazione[i].Swap((int)(rnd.Rannyu()*N_city),(int)(rnd.Rannyu()*N_city));
  popolazione[i].Verify();
 }

//Calcolo la distanza di ogni percorso.
 costFunction L;
 for(int i=0; i<N_per; i++)  popolazione[i].SetDistanza(L.cost(popolazione[i],ct,N_city));

//Ordino i percorsi in base alla distanza.
  Ordina(popolazione,N_per,N_city);
  cout << "Distanze iniziali ordinate in base alla distanza" << endl;
  for(int i=0; i<N_per; i++){ 
    //popolazione[i].Print(); 
    //cout << endl << i << " " << popolazione[i].GetDistanza() << endl; 
  }


 crossover acc;
//Creo una nuova generazione.
 int N_gen=100;
 percorso* next_generation;
 ofstream fileout, fileout2, out_conf;
 fileout.open("LBestPath.dat");
 fileout2.open("L-medie.dat");

 for(int k=0; k < N_gen; k++){
   next_generation=new percorso[N_per];
   for(int i=0; i<N_per; i++) next_generation[i].Inizializza(N_city);

   for(int i=0; i<N_per-1; i=i+2){
    //Seleziono due genitori a caso. La quarta potenza aiuta a selezionare i migliori genitori.
    int k=(int)(N_per*pow(rnd.Rannyu(),3.)); 
    int n=(int)(N_per*pow(rnd.Rannyu(),3.));
    acc.Input(popolazione[k],popolazione[n],N_city,5);
    acc.Accoppiamento(&rnd);
    //Copio il vettore del percorso (ma non la distanza!).
    next_generation[i+1]=acc.GetFirtsSon();
    //Calcolo ed inserisco la distanza del primo figlio nel primo figlio.
    next_generation[i+1].SetDistanza(L.cost(next_generation[i+1],ct,N_city));
    next_generation[i+1].Verify();
    //Lo stesso per il secondo figlio.
    next_generation[i]=acc.GetSecondSon();
    next_generation[i].SetDistanza(L.cost(next_generation[i],ct,N_city));
    next_generation[i].Verify();
    //Mutazione.
    next_generation[i+1].Mutazione(&rnd);
    next_generation[i+1].SetDistanza(L.cost(next_generation[i+1],ct,N_city));
    next_generation[i+1].Verify();
    next_generation[i].Mutazione(&rnd);
    next_generation[i].SetDistanza(L.cost(next_generation[i],ct,N_city));
    next_generation[i].Verify();
    //cout << "    next_generation[i+1].Print() " << endl; 
    //next_generation[i+1].Print();
    //cout << "    next_generation[i].Print() " << endl; 
    //next_generation[i].Print();
   }

  next_generation[0]=popolazione[0]; //La migliore mamma e il migliore papà di ogni generazione sopravvivono.
  next_generation[0].SetDistanza(popolazione[0].GetDistanza());

  //Salvo i dati su file.
  fileout << k << " " << next_generation[0].GetDistanza() << endl;
  double media=0.;
  for(int i=0; i<N_gen/2; i++) media+=next_generation[i].GetDistanza();
  fileout2 << k << " " << media/((double)N_gen/2.) << endl;

  //Salvo il path per diverse generazioni.
  if(k%9 == 0){
    //cout << k << endl;
    out_conf.open("./Configurazioni/configurazione" + to_string(k) + ".dat");
    for(int j=0; j<N_city; j++){
      int a=next_generation[0].Get(j);
      out_conf << ct[a].GetX() << " " <<  ct[a].GetY() << endl;
    }
    out_conf << ct[next_generation[0].Get(0)].GetX() << " " <<  ct[next_generation[0].Get(0)].GetY() << endl;
    out_conf.close();
  }

  Ordina(next_generation,N_per,N_city); 

  //Distruggo la vecchia popolazione.
  for(int i=0; i<N_per; i++){
    popolazione[i]=next_generation[i];
    popolazione[i].SetDistanza(next_generation[i].GetDistanza());
  }
}

  cout << "Configurazione finale " << endl;

  for(int i=0; i<N_per; i++){ 
  //popolazione[i].Print();
   cout << i << " " << popolazione[i].GetDistanza() << endl; 
 }

 fileout.close();
 fileout2.close();
 out_conf.close();

return 0;
}


percorso* Ordina(percorso* popolazione, int N_per,int N_city){
//Ordino i percorsi in base alla distanza.
 int i, j, min;
 int* temp_per=new int[N_city];
 double temp_lun;

 for (i = 0; i < N_per - 1; i++) {
  min = i;
  for (j = i + 1; j < N_per; j++)
   if(popolazione[j].GetDistanza() < popolazione[min].GetDistanza()) min = j;
  
  for(int k=0; k<N_city; k++) temp_per[k]=popolazione[i].Get(k);
  temp_lun=popolazione[i].GetDistanza(); 

  for(int k=0; k<N_city; k++) popolazione[i].SetComponent(popolazione[min].Get(k),k);
  popolazione[i].SetDistanza(popolazione[min].GetDistanza());

  for(int k=0; k<N_city; k++) popolazione[min].SetComponent(temp_per[k],k); 
  popolazione[min].SetDistanza(temp_lun);
 }

return popolazione; 
}




