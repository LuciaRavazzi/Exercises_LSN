
main.exe : main.o random.o blocchi.o
	g++ random.o main.o blocchi.o -o main.exe
main.o : main.cpp random.h blocchi.h
	g++ -c main.cpp -o main.o 
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
blocchi.o: blocchi.cpp blocchi.h random.h
	g++ -c blocchi.cpp -o blocchi.o
clean :
	rm *.o main.exe seed.out
