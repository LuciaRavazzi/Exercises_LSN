#ifndef __city__
#define __city__

#include "random.h"

//Una sola città.

class city{	
 public:
 void Set(Random*); 
 double GetX();
 double GetY();
 double SetX(double);
 double SetY(double);

 private:
 double theta;
 double x,y;  
};

#endif
