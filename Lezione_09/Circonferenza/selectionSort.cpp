#include "selectionSort.h"

// (10) ORDINA ARRAY DI double (M.1): ordina un array di double di dimensione passata in ordine crescente. Algoritmo SelectionSort
void ordinaArrayCrescente(double array[],int dim){
	// Utilizza un ciclo for. Il ciclo controlla il vettore per individuarne il minimo, ma trascura le parti che precedono il
	// proprio indice. Sposta poi il valore all'inizio dell'array.
	int a=0;
	for (int i=0; i<dim; i++){
		a=trovaMinimoDa(array,dim,i);
		scambiaValori(array,a,i);
	}
}

// (09) TROVA MINIMO A PARTIRE DA IN UN ARRAY DI double: restituisce l'indice di una delle occorrenze del minimo in un array di double a partire dall'indice passato.
int trovaMinimoDa(double array[], int dim, int da){
	double valLoser=array[da];
	int idLoser=da;
	
	for (int i=da; i<dim; i++){
		if (array[i]<valLoser){
			valLoser=array[i];
			idLoser=i;
		}
	}
	
	return idLoser;
}

// (07) SCAMBIA VALORI IN ARRAY DI double: scambia tra loro due valori di indice passato in un array di double.
void scambiaValori(double array[],int a,int b){
	double valA,valB;
	valA=array[a];
	valB=array[b];
	array[a]=valB;
	array[b]=valA;
}
