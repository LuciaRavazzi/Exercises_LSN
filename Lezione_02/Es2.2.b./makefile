
main.exe:main.o random.o randomWalk.o 
	g++ main.o random.o randomWalk.o -o main.exe

main.o: main.cpp randomWalk.h random.h main.h
	g++ -c main.cpp -o main.o

randomWalk.o: randomWalk.cpp randomWalk.h random.h
	g++ -c randomWalk.cpp -o randomWalk.o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o





