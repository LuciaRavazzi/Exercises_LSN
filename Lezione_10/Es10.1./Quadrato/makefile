main.exe: main.o city.o random.o percorso.o costFunction.o crossover.o
	g++ main.o city.o random.o percorso.o costFunction.o crossover.o -o main.exe

main.o: main.cpp percorso.h city.h costFunction.h crossover.h
	g++ -c main.cpp -o main.o

city.o: city.cpp city.h random.h 
	g++ -c city.cpp -o city.o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

percorso.o: percorso.cpp percorso.h random.h
	g++ -c percorso.cpp -o percorso.o

costFunction.o: costFunction.cpp costFunction.h city.h percorso.h
	g++ -c costFunction.cpp -o costFunction.o

crossover.o:crossover.cpp crossover.h percorso.h random.h
	g++ -c crossover.cpp -o crossover.o
