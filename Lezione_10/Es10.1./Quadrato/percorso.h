#ifndef __percorso__
#define __percorso__

#include "random.h"

class percorso{	
 public:
  percorso& operator= (const percorso &popolazione);
  void Inizializza(int);
  void Swap(int,int); 
  void Print();
  void Verify();
  int Get(int);
  void SetDistanza(double);
  double GetDistanza();
  void MutazioneShift(Random*);
  void MutazionePermutazione(Random*);
  void SetComponent(int,int);
  int pbc(int i);

  int N;
  private:
  int* tragitto;
  double distanza;
};

#endif
