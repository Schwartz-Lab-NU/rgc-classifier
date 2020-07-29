CBLAS = "c:/msys64/mingw64/include/OpenBLAS"
FLAGS =

all: test train print ecoc
debug: FLAGS += -g
debug: test train print ecoc

train: treeObjects.o train.cpp
	g++ -m64 -std=c++17 -static-libgcc -static-libstdc++ -fopenmp -O3 -I$(CBLAS) train.cpp treeObjects.o -lblas -lm -o train $(FLAGS)

test: treeObjects.o test.cpp
	g++ -m64 -std=c++17 -static -O3 -I$(CBLAS) test.cpp treeObjects.o -lblas -lm -o test $(FLAGS)

print: treeObjects.o print.cpp
	g++ -m64 -std=c++17 -static -I$(CBLAS) print.cpp treeObjects.o -lblas -lm -o print $(FLAGS)

ecoc: treeObjects.o ecoc.cpp
	g++ -m64 -std=c++17 -static -I$(CBLAS) ecoc.cpp treeObjects.o -lblas -lm -o ecoc $(FLAGS)	

treeObjects.o: treeObjects.cpp
	g++ -m64 -std=c++17 -static -O3 -c treeObjects.cpp $(FLAGS)

clean:
	-rm -f train test print ecoc treeObjects.o