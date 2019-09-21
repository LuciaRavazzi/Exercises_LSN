#ifndef __blocchi__
#define __blocchi__


class Blocchi {
	public:
	Blocchi(); //costruttore
	~Blocchi(); //distruttore

	void setParametri(int, int);
	void mediaBlocco(double, int);
	void erroreBlocco();
	double error(double*,double*,int );
	void fileFinale();
	
	int M; //numero di tentativi totali.
	int N; //numero di blocchi totali.
	int L; //quanti tentativi nel blocco.
	double* media=new double[N];
	double* media2=new double[N];
	double* sum_prog=new double[N];
	double* su2_prog=new double[N];
	double* err_prog=new double[N];
};


#endif
