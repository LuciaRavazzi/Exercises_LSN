#ifndef __randomWalk__
#define __randomWalk__

#include "random.h"  //necessario: altrimenti mi dà errore.

class randomWalk{
 public:
   randomWalk(double);
   void azzero();//azzero le somme in modo tale da poter fare un nuovo RW.
   double move(int);//faccio un random Walk.
   void setPasso(double);//definisco il passo.	

 private:
   double passo;//lunghezza di ogni passo del random Walk.
   double sumx,sumy,sumz;//somme progressive dei passi lungo una certa direzione.
   Random rnd;
};

#endif 
