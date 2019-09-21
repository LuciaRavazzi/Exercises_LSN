#ifndef __blocchi__
#define __blocchi__


class blocchi {
 public:
  blocchi(int,int); //costruttore
  ~blocchi(); //distruttore

  void SetParametri(int, int);
  void MediaBlocco(double, int);
  void ErroreBlocco();
  double Error(double*,double*,int );
  void Print();

	
  int M; //numero di tentativi totali.
  int N; //numero di blocchi totali.
  int L; //quanti tentativi nel blocco.
  double* media=new double[N];
  double* sum_prog=new double[N];
  double* su2_prog=new double[N];
  double* err_prog=new double[N];
};


#endif
