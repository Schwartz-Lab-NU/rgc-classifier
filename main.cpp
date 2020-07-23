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
	std::cout << "Using parameters from params" + (std::string)argv[1] + ".in" << std::endl;
	std::fstream paramfile((std::string)argv[1] + "params.in", std::ios::in);
	paramSet p(paramfile);
	p.print();

	//load the error correction matrix
	int* ecoc = new int[p.nLabels*p.ensembleSize];
	for (int i=0;i<p.nLabels*p.ensembleSize;i++) {
		paramfile >> ecoc[i];
	}


	int f;
	double ff;
	std::string foldname;
	if (paramfile.eof()) {
		f = 1;
		foldname = "fold1.in";
	} else {
		paramfile >> ff;
		// std::cout << ff << std::endl;
		if (paramfile.eof()) {
			f = 2;
			foldname = "fold2.in";
		} else {
			paramfile >> ff;
			// std::cout << ff << std::endl;
			if (paramfile.eof()) {
				f = 3;
				foldname = "fold3.in";
			} else {
				std::cout << "All folds already trained. Exiting." << std::endl;
				delete[] ecoc;
				return 0;
			}
		}
	}

	paramfile.close();


	//load the PSTH data
	std::ifstream infile(foldname,std::ios::binary);

	psthSet* data;
	data = new psthSet(infile);

	infile.close();


	std::filesystem::path rootDir = argv[1];
	rootDir /= ("fold" + std::to_string(f));

	int NT = atoi(argv[2]); //number of threads
	
	std::filesystem::create_directories(rootDir); //create the root dir if it doesn't exist yet

	//begin timer
	// std::clock_t start = std::clock();
	struct timeval start, end;
	gettimeofday(&start, NULL);

	std::cout << "Training on fold " << f << std::endl;


	#pragma omp parallel num_threads(NT) //private(ensl, datal)
	{
		int t = omp_get_thread_num();
		// int err = 0;

		Ensemble ensl(&p, data->N, t, rootDir, ecoc);
		psthSet datal(data);

		#pragma omp for schedule(dynamic, 1)
			for (int i=0; i<p.ensembleSize; i++) {
				ensl.train(&datal, i);
			}
		// std::cout << "exiting on thread " << t << std::endl;
	}
	gettimeofday(&end, NULL);

	std::cout << "Done training on fold " << f << "!" << std::endl;

	paramfile.open((std::string)argv[1] + "params.in", std::ios::app);
	paramfile << std::endl << (end.tv_sec - start.tv_sec)*NT; //time * thread
	paramfile.close();

	delete[] ecoc;
	delete data;


	// std::cout << "Line " << __LINE__  << ", File " << __FILE__ <<std::endl; //kept here for reuse
	return 0;
}