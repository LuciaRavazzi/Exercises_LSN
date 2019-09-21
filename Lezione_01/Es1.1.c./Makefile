
main.exe : main.o random.o
	g++ random.o main.o -o main.exe
main.o : main.cpp
	g++ -c main.cpp -o main.o 
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o main.exe seed.out
