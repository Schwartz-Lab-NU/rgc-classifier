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
#include <iomanip>

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
	int nFolds = 3;

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
	paramfile.close();

	psthSet** data = new psthSet*[nFolds];
	int** sampleCounts = new int*[nFolds];
	std::filesystem::path* rootDir = new std::filesystem::path[3];
	// int N;

	for (int f=0;f<nFolds;f++) {
		std::string foldname = "fold" + std::to_string(f+1) + ".in";
		std::ifstream infile(foldname,std::ios::binary);

		//load the PSTH data
		data[f] = new psthSet(infile);
		infile.close();

		// N = std::max(N, data[f]->N);

		sampleCounts[f] = new int[p.nLabels]();
		for (int i=0; i<data[f]->N; i++) {
			sampleCounts[f][data[f]->data[i]->label -1] +=1;
		}

		rootDir[f] = argv[1];
		rootDir[f] /= ("fold" + std::to_string(f+1));
		
		std::filesystem::create_directories(rootDir[f]); //create the root dir if necessary
	}

	int NT = atoi(argv[2]); //number of threads
	

	//begin timer
	// std::clock_t start = std::clock();
	struct timeval start, end;
	gettimeofday(&start, NULL);

	std::cout << "Beginning training... " << std::endl;

	#pragma omp parallel num_threads(NT) //private(sampleCounts)//private(ensl, datal)
	{
		int t = omp_get_thread_num();
		// int err = 0;

		for (size_t f=0; f<nFolds; f++) {
			Ensemble ensl(&p, data[f]->N, t, rootDir[f], ecoc);
			psthSet datal(data[f]);

			#pragma omp for schedule(dynamic, 1) nowait
				for (int i=0; i<p.ensembleSize; i++) {
					ensl.train(&datal, sampleCounts[f], i);

					#pragma omp critical
					std::cout << "Fold[" << f+1 << "]:" << std::endl << ensl.buffer.rdbuf();
					//only print from 1 thread at a time, to avoid illegible results
				}
		}
		// std::cout << "exiting on thread " << t << std::endl;
	}
	gettimeofday(&end, NULL);

	std::cout << "Done training!" << std::endl;

	paramfile.open((std::string)argv[1] + "params.in", std::ios::out);

	p.print(paramfile); //re-write the params...

	paramfile << std::endl;
	for (int i=0;i<p.ensembleSize;i++) {
		int stride = i*p.nLabels;
		for (int j=0;j<p.nLabels;j++) {
			paramfile << std::setw(3)<< ecoc[stride+j];
		}
		paramfile << std::endl;
	} //re-write the ecoc matrix...

	paramfile << (int) (end.tv_sec - start.tv_sec)*NT; //time * thread
	paramfile.close();

	delete[] ecoc;
	delete[] sampleCounts;
	delete data;
	// delete data2;


	// std::cout << "Line " << __LINE__  << ", File " << __FILE__ <<std::endl; //kept here for reuse
	return 0;
}