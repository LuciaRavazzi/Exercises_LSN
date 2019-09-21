
qmc1d.exe: qmc1d.o random.o 
	g++ qmc1d.o random.o -o qmc1d.exe

qmc1d.o: qmc1d.cpp functions.h constants.h random.h
	g++ -c qmc1d.cpp -o qmc1d.o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
