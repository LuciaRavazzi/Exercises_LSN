
MolDyn_NVE.exe: MolDyn_NVE.o random.o
	g++ MolDyn_NVE.o random.o -o MolDyn_NVE.exe

MolDyn_NVE.o: MolDyn_NVE.cpp MolDyn_NVE.h random.h 
	g++ -c MolDyn_NVE.cpp -o MolDyn_NVE.o 

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
