#ifndef __Ago__
#define __Ago__

#include "random.h"

class Ago {
	public:
	Ago(); //costruttore
	~Ago(); //distruttore
	void setParametri(int, int, int);
	double pigreco();

	int L; //lunghezza ago.
	int D; //distanza tra le due assi del pavimento.
	int lanci; //numero dilanci per ottenere una stima di pigreco.
	Random *rnd;
};
#endif
