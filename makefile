CBLAS = "c:/msys64/mingw64/include/OpenBLAS"
FLAGS =

all: FLAGS += -O3
all: test transform train print ecoc
debug: FLAGS += -g
debug: test transform train print ecoc

train: treeObjects.o train.cpp
	g++ -m64 -std=c++17 -static-libgcc -static-libstdc++ -fopenmp -I$(CBLAS) train.cpp treeObjects.o -lblas -lm -o train $(FLAGS)

transform: treeObjects.o transform.cpp
	g++ -m64 -std=c++17 -static -I$(CBLAS) transform.cpp treeObjects.o -lblas -lm -o transform $(FLAGS)

test: treeObjects.o test.cpp
	g++ -m64 -std=c++17 -static -I$(CBLAS) test.cpp treeObjects.o -lblas -lm -o test $(FLAGS)

print: treeObjects.o print.cpp
	g++ -m64 -std=c++17 -static -I$(CBLAS) print.cpp treeObjects.o -lblas -lm -o print $(FLAGS)

ecoc: treeObjects.o ecoc.cpp
	g++ -m64 -std=c++17 -static -I$(CBLAS) ecoc.cpp treeObjects.o -lblas -lm -o ecoc $(FLAGS)	

treeObjects.o: treeObjects.cpp
	g++ -m64 -std=c++17 -static -c treeObjects.cpp $(FLAGS)

clean:
	-rm -f train transform test print ecoc treeObjects.o