#ifndef __city__
#define __city__

#include "random.h"

//Una sola città.

class city{	
 public:
 void Set(double,double); 
 double GetX();
 double GetY();
 double GetAngle();
 void SetAngle(double);

 private:
 double theta;
 double x,y;  
};

#endif
