#include <iostream>
#include "treeObjects.h"
#include <fstream>
#include <string>

//usage: ./test PARAMS_ID train_fold test_fold transform_fold
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

	char transform;
	if (argc>4) {
		transform = (char) atoi(argv[4]); //lol
	} else {
		transform = 0;
	}

	char* pred = new char[data->N];

	// for (int i=0; i<data->N; i++) {
	// 	std::cout << (int) data->data[i]->label << ",";
	// }
	// std::cout << std::endl;
	int* sampleCounts = new int[p.nLabels]();
	
	for (int i = 0; i<data->N; i++) {
		sampleCounts[data->data[i]->label - 1]+=1;
	}	

	std::cout << "Testing model trained from fold " + (std::string) argv[2] + " on fold " + (std::string) argv[3];
	if (transform) {
		std::cout << " using transform from fold " + (std::string) argv[4];
	}
	std::cout << std::endl;

	double loss = ens->test(data, pred, sampleCounts, transform);

	// for (int i=0; i<data->N; i++) {
	// 	sampleCounts[data->data[i]->label -1] +=1;
	// }

	double acc = 0.;
	int* classAcc = new int[p.nLabels]();
	for (int i = 0; i<data->N; i++) {
		// std::cout <<"("<<(int) data->data[i]->label<<","<<(int) pred[i]+1 << "),";
		if (data->data[i]->label - 1 == pred[i]) {
			acc += 1.;//./data->N;
			classAcc[data->data[i]->label - 1] +=1;
		}
	} 
	std::cout << "Total accuracy: " << 100.*acc / data->N<<"%" << std::endl;

	acc=0.;
	for (int i =0; i<p.nLabels; i++) {
		acc+= (double)classAcc[i]/sampleCounts[i];
	}

	std::cout << "Class-weighted accuracy: " << 100*acc / p.nLabels <<"%" <<std::endl;


	paramfile.open((std::string)argv[1] + "params.in", std::ios::app);
	paramfile << std::endl << argv[2] << ","<< argv[3];
	if (transform) {
		paramfile << "," << argv[4] << "," << loss;
	} else {
		paramfile << "," << 100*acc / p.nLabels;	
	}
	paramfile.close();

	std::cout << "Successfully wrote loss to file" << std::endl;

	delete[] sampleCounts;	
	delete[] classAcc;
	delete[] ecoc;
	delete data;
	delete ens;
	delete[] pred;

	return 0;
}