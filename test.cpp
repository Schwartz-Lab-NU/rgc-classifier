#include <iostream>
#include "treeObjects.h"
#include <fstream>
#include <string>

//usage: ./test PARAMS_ID train_fold test_fold
int main(int argc, char* argv[]) {

	//load the params
	std::cout << "Using parameters from " << (std::string)argv[1] + "params.in" << std::endl;
	std::fstream paramfile((std::string)argv[1] + "params.in", std::ios::in);
	paramSet p(paramfile);
	p.print();

	//load the error correction matrix
	int* ecoc = new int[p.nLabels*p.ensembleSize];
	for (int i=0;i<p.nLabels*p.ensembleSize;i++) {
		paramfile >> ecoc[i];
	}
	paramfile.close();

	// for (int i=0;i<p.ensembleSize; i++) {
	// 	for (int j=0;j<p.nLabels; j++) {
	// 		std::cout << ecoc[i*p.nLabels+j] << ",";
	// 	}
	// 	std::cout << std::endl;
	// }
	// std::cout << std::endl;

	//load test fold data.in
	std::ifstream infile("fold" + (std::string)argv[3] + ".in",std::ios::binary);

	psthSet* data;
	data = new psthSet(infile);
	infile.close();

	Ensemble* ens = new Ensemble;
	ens->params = &p;
	ens->N = data->N;
	ens->ecoc=ecoc;
	ens->rootDir = (std::string)argv[1] + "/fold" + (std::string) argv[2] + "/";


	char* pred = new char[data->N];
	ens->test(data, pred);

	double acc = 0.;
	for (int i = 0; i<data->N; i++) {
		std::cout <<"("<<(int) data->data[i]->label<<","<<(int) pred[i]+1 << "),";
		if (data->data[i]->label == pred[i]+1) {
			acc += 1./data->N;
		}
	} 
	std::cout << std::endl << acc << std::endl;

	return 0;
}