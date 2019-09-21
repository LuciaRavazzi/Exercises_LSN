#include <iostream>
#include <fstream>
#include "random.h"
#include <cmath>
#include "main.h"


using namespace std;

int main(){
  mu=0.8;
  sigma=0.6;

  Input();
  Equilibrio();

  cout <<"Passo: " <<  delta << endl;

  for(int j=0; j<N; j++){
    Reset(); 
    for(int i=0; i<L; i++){
      Move();	
      Measure();
      Accumulate();  
   }
   Averages(j);
  }
 
 cout << endl << endl;
 Statistica();
 cout << "Media : " << mu << endl;

 Print();



return 0;
}

void Input(void){
  cout << endl << endl;
  cout << "	Ricerca dell'energia di ground state di una particella quantistica in una dimensione	" << endl << endl;
  x=0.; //Posizione iniziale.
  cout << "	Posizione iniziale x=" << x << endl;
}

void Equilibrio(){
  int B=1000;
  //cout << "		Sto cercando il passo per avere un'accettanza del 50%		" << endl;
  delta=0.;
  do{
    delta=delta+0.01;
    counter=0;
    for (int i=0; i<B; i++) Move();
    accettanza=(double)counter/B;
    //cout << "Passo: " << delta << " accettanza: " << accettanza << endl;
  } while (accettanza > 0.52 || accettanza < 0.47);
  //cout << "Il passo che mi dà un'accettanza del 50% nel Metropolis vale " << delta << endl;
  cout << endl << endl << endl << endl;

}

void Move(){     
 xnew=rnd.Rannyu(x-delta, x+delta);
 w=(double)Density(xnew)/Density(x);

 if(w >= 1)
 {
  x=xnew;
  counter++;
 } else {
  if(rnd.Rannyu() <= w)
  { 
   x=xnew;
   counter++;
   } 
 }
}


double Density(double x){
  return exp(-pow(x-mu,2)/pow(sigma,2))+exp(-pow(x+mu,2)/pow(sigma,2))+2.*exp(-(pow(x,2)+pow(mu,2))/pow(sigma,2));
}

double Wave_function(double x){
  return exp(-pow(x-mu,2)/(2.*pow(sigma,2)))+exp(-pow(x+mu,2)/(2.*pow(sigma,2)));
}

double Ekin(double x){
  pezzo1=exp(-pow(x-mu,2)/(2.0*pow(sigma,2)));
  pezzo2=exp(-pow(x+mu,2)/(2.0*pow(sigma,2)));
  pezzo3=1.0-(pow(x-mu,2)/pow(sigma,2));
  pezzo4=1.0-(pow(x+mu,2)/pow(sigma,2));

return ((1./(2.*pow(sigma,2)))*(pezzo1*pezzo3+pezzo2*pezzo4))/(double)Wave_function(x);
}

double Pot(double x){
  return ((pow(x,4)-2.5*pow(x,2))*Wave_function(x))/Wave_function(x);
}

void Reset(){
 sum=0.;
 counter=0;
 for(int i=0; i<N; i++) {
  media[N]=0.;
  media_prog[N]=0.;
  media2_prog[N]=0.;
  err_prog[N]=0.;
 } 
}

void Measure(){
  etot=Ekin(x)+Pot(x);
}

void Accumulate(void){
 sum += etot;
}

void Averages(int j){
 media[j]=sum/(double)L;	
}

void Statistica(void){
  for(int i=0; i<N; i++){
    media_prog[i]=0.;
    media2_prog[i]=0.;
    for(int j=0; j<i+1; j++){
      media_prog[i]+=media[j]; //sommo le medie
      media2_prog[i]+=media[j]*media[j]; //sommo le medie al quadrato
    }

    media_prog[i]=media_prog[i]/(float)(i+1); //medie progressive	
    media2_prog[i]=media2_prog[i]/(float)(i+1); //media progressive sui quadrati	
    err_prog[i]=error(i); //incertezza statistica
     
   // cout << i << "media[i] " << media[i] << " media_prog[i] " << media_prog[i]<< " " << endl;
   // cout << i << "pow(media_prog[i],2) " << pow(media_prog[i],2) << " media2_prog[i] " << media2_prog[i]<< " " << endl;
    cout << i << " " << media_prog[i] << " " << err_prog[i] << endl;
   // cout << endl;
  }
}

double error(int i){
  if(i==0){
    return 0;
  } else {
    return sqrt((media2_prog[i]-media_prog[i]*media_prog[i])/i);
  }
}

void Print(){
 ofstream fileout;
 fileout.open("Statistica");

 for(int i=0; i<N; i++)
    fileout << i << " " << media_prog[i] << " " << err_prog[i] << endl;
   
 
 fileout.close();
}

