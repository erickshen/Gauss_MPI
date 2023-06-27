all: Gauss

Gauss: main.o ErrorNorm.o MatrixNorm.o MatrixInput.o MatrixOutput.o Function.o Gauss.o
	mpic++ main.o ErrorNorm.o MatrixNorm.o MatrixInput.o MatrixOutput.o Function.o Gauss.o -lpthread -o Gauss

main.o: main.cpp
	mpic++ -c main.cpp -lpthread -O3

ErrorNorm.o: ErrorNorm.cpp
	mpic++ -c ErrorNorm.cpp -lpthread -O3

MatrixNorm.o: MatrixNorm.cpp
	mpic++ -c MatrixNorm.cpp -lpthread -O3

MatrixInput.o: MatrixInput.cpp
	mpic++ -c MatrixInput.cpp -lpthread -O3

MatrixOutput.o: MatrixOutput.cpp
	mpic++ -c MatrixOutput.cpp -lpthread -O3

Function.o: Function.cpp
	mpic++ -c Function.cpp -lpthread -O3

Gauss.o: Gauss.cpp
	mpic++ -c Gauss.cpp -lpthread-O3

clean: 
	rm -rf *.o Gauss
