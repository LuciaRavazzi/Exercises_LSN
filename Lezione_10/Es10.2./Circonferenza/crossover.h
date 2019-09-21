#ifndef __crossover__
#define __crossover__

#include "percorso.h"
#include "random.h"

class crossover{	
 public:
 void Input(percorso, percorso,int,int);
 void Accoppiamento(Random*);
 percorso GetFirtsSon();
 percorso GetSecondSon();
 double GetDistanzaFirtsSon();
 double GetDistanzaSecondSon();
 void Swap(int); 


 private:
 int N_city;
 int quanti_scambi;
 percorso Mom;
 percorso Dad;
 percorso FirstSon;
 percorso SecondSon;
};

#endif
