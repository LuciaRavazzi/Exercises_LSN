#include <fstream>

using namespace std;
 
 const double siniz=100;
 const double T=1;
 const double K=100;
 const double r=0.1;
 const double vol=0.25;
 int appo=100;

 double sum=0;
 double c=0;	
 int N=100;//numero blocchi
 int M=100000;//numero lanci totali
 double L=(double) M/N;
 double* call=new double[N];
 double* sum_prog=new double[N];
 double* su2_prog=new double[N];
 double* err_prog=new double[N];
 double passo=(double)T/appo;
 double s,p;
 double t=0;

void Statistica();
double Sdisc(double ,double , double ,double );	
double error(double* ,double* ,int );
double Scont(double t, double c);
double fmax(double a);





