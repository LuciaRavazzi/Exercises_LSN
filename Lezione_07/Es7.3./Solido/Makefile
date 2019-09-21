
MolDyn_NVE.exe : MolDyn_NVE.o random.o
	g++ random.o MolDyn_NVE.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp
	g++ -c MolDyn_NVE.cpp -o MolDyn_NVE.o 
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o MolDyn_NVE.exe
