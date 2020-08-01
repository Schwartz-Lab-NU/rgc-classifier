#include <iostream>
#include "treeObjects.h"
#include <fstream>
#include <string>

//usage: ./transform PARAMS_ID
int main(int argc, char* argv[]) {

	//load the params
	std::cout << "Using parameters from " << (std::string)argv[1] + "params.in" << std::endl;
	std::fstream paramfile((std::string)argv[1] + "params.in", std::ios::in);
	paramSet p(paramfile);
	// p.print();

	//load the error correction matrix
	int* ecoc = new int[p.nLabels*p.ensembleSize];
	for (int i=0;i<p.nLabels*p.ensembleSize;i++) {
		paramfile >> ecoc[i];
	}

	int f = 1;
	// double ff;

	int c = 0;
	c = paramfile.peek();
	while (c!=EOF) {
		paramfile >> c;
		f++;
		if (f>12) {
			std::cout << "All folds already transformed. Exiting." << std::endl;
			delete[] ecoc;
			return 0;
		}
		c = paramfile.peek();
	}
	paramfile.close();

	if (f<4) {
		std::cout << "Folds are not done training! Exiting." << std::endl;
		delete[] ecoc;
		return 0;
	}

	// int trainFold = ((f-4)/2) + 1; //fold to train transform on
	// int transformFold = ((f % 2) + trainFold + 1) % 3; //fold to test transform against
	// if (transformFold==0) {
	// 	transformFold = 3;
	// }
	int trainFold = ((f-4)/3) + 1;
	int transformFold = (f-4)%3 + 1;

	// std::string foldname = "fold" + std::to_string(f) + ".in";
	//load test fold data.in
	std::ifstream infile("fold" + std::to_string(transformFold) + ".in",std::ios::binary);

	psthSet* data;
	data = new psthSet(infile);
	infile.close();


	int* sampleCounts = new int[p.nLabels]();
	for (int i=0; i<data->N; i++) {
		sampleCounts[data->data[i]->label -1] +=1;
	}


	Ensemble* ens = new Ensemble;
	ens->params = &p;
	ens->N = data->N;
	ens->ecoc=ecoc;
	ens->rootDir = (std::string)argv[1] + "/fold" + std::to_string(trainFold) + "/";

	// char* pred = new char[data->N];

	std::cout << "Normalizing output from fold "<< trainFold<< " using labels from fold " << transformFold << std::endl ;

	ens->normalize(data, sampleCounts, transformFold);

	delete ens;
	delete data;
	delete[] sampleCounts;
	delete[] ecoc;

	paramfile.open((std::string)argv[1] + "params.in", std::ios::app);
	paramfile << std::endl << (int) 1; //time * thread
	paramfile.close();

	std::cout << "Transform successful" << std::endl;
	

	return 0;
}