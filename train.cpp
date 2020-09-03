#include <iostream>
#include "treeObjects.h"
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <omp.h>
#include <stdint.h>
#include <filesystem>
#include <ctime>
#include <sys/time.h>
#include <sstream>

//usage: ./train PARAMS_ID num_threads
int main(int argc, char* argv[]) {
	// test compiler version
	// #if UINTPTR_MAX == 0xffffffff
	// std::cout << "32bit compile" <<std::endl;
	// #elif UINTPTR_MAX == 0xffffffffffffffff
	// std::cout << "64bit compile" <<std::endl;
	// #else
	// std::cout << "unknown compile" <<std::endl;
	// #endif

	//load the parameters from the file
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
		if (f>3) {
			std::cout << "All folds already trained. Exiting." << std::endl;
			delete[] ecoc;
			return 0;
		}
		c = paramfile.peek();
	}
	paramfile.close();

	std::string foldname = "fold" + std::to_string(f) + ".in";
	// std::string foldname = "fold" + std::to_string(f==1?2:1) + ".in";
	//load the PSTH data
	std::ifstream infile(foldname,std::ios::binary);

	psthSet* data;
	data = new psthSet(infile);
	infile.close();

	//load the second dataset

	// foldname = "fold" + std::to_string(f==3?2:3) + ".in";
	// infile.open(foldname, std::ios::binary);
	// psthSet* data2;
	// data2 = new psthSet(infile);
	// infile.close();

	// data->append(data2);

	int* sampleCounts = new int[p.nLabels]();
	for (int i=0; i<data->N; i++) {
		sampleCounts[data->data[i]->label -1] +=1;
	}

	std::filesystem::path rootDir = argv[1];
	rootDir /= ("fold" + std::to_string(f));

	int NT = atoi(argv[2]); //number of threads
	
	// std::filesystem::remove_all(rootDir); //we want to clear the prior contents of this directory if any
	std::filesystem::create_directories(rootDir); //re-create the root dir

	//begin timer
	// std::clock_t start = std::clock();
	struct timeval start, end;
	gettimeofday(&start, NULL);

	std::cout << "Training against fold " << f << std::endl;


	#pragma omp parallel num_threads(NT) //private(sampleCounts)//private(ensl, datal)
	{
		int t = omp_get_thread_num();
		// int err = 0;

		Ensemble ensl(&p, data->N, t, rootDir, ecoc);
		psthSet datal(data);

		#pragma omp for schedule(dynamic, 1)
			for (int i=0; i<p.ensembleSize; i++) {
				ensl.train(&datal, sampleCounts, i);

				#pragma omp critical
				std::cout << ensl.buffer.rdbuf();
			}
		// std::cout << "exiting on thread " << t << std::endl;
	}
	gettimeofday(&end, NULL);

	std::cout << "Done training against fold " << f << "!" << std::endl;

	paramfile.open((std::string)argv[1] + "params.in", std::ios::app);
	paramfile << std::endl << (int) (end.tv_sec - start.tv_sec)*NT; //time * thread
	paramfile.close();

	delete[] ecoc;
	delete[] sampleCounts;
	delete data;
	// delete data2;


	// std::cout << "Line " << __LINE__  << ", File " << __FILE__ <<std::endl; //kept here for reuse
	return 0;
}