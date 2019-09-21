
Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	g++ random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	g++ -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o 
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o 
clean :
	rm *.o Monte_Carlo_NVT.exe
