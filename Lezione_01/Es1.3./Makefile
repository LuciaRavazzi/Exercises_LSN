
main.exe : main.o random.o blocchi.o ago.o
	g++ random.o main.o blocchi.o ago.o -o main.exe

main.o : main.cpp random.h blocchi.h main.h ago.h
	g++ -c main.cpp -o main.o 

random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 

blocchi.o: blocchi.cpp blocchi.h random.h
	g++ -c blocchi.cpp -o blocchi.o

ago.o: ago.cpp ago.h random.h
	g++ -c ago.cpp -o ago.o
clean :
	rm *.o main.exe seed.out
