#include <iostream>
#include "treeObjects.h"
#include <fstream>
#include <string>
#include <cstring>
#include <climits>
#include <iomanip>

int do_ecoc(int*, int, int, double, double);
bool unique_rows(int*, int, int);
bool traversible(int*, int, int);
bool visit(int, int, int, int*, bool*, bool*, int, int, int);

#define MAX_ITER 10000

//usage: ./ecoc PARAMS_ID 
int main(int argc, char* argv[]) {

	//load the params
	std::cout << "Using parameters from " << (std::string)argv[1] + "params.in" << std::endl;
	std::fstream paramfile((std::string)argv[1] + "params.in", std::ios::in);
	paramSet p(paramfile);
	p.print();

	if (!paramfile.eof()) {
		std::cout << "ECOC matrix already exists. Previous results will be overwritten." << std::endl;
		// paramfile.close();
		// return 0;
	}
	paramfile.close();

	//load the error correction matrix
	int* ecoc = new int[p.nLabels*p.ensembleSize];
	
	int error = do_ecoc(ecoc, p.nLabels, p.ensembleSize, p.probNeg, p.probPos);

	if (error) {
		delete[] ecoc;
		std::cout << "Error: unable to find valid ECOC matrix" << std::endl;
		return 0;
	}


	std::ofstream paramsout;
	paramsout.open((std::string)argv[1] + "params.in");
	p.print(paramsout); //re-write the params...

	paramsout << std::endl;
	for (int i=0;i<p.ensembleSize;i++) {
		int stride = i*p.nLabels;
		for (int j=0;j<p.nLabels;j++) {
			paramsout << std::setw(3)<< ecoc[stride+j];
		}
		if (i<p.ensembleSize-1) {
			paramsout << std::endl;
		}	
	}

	paramsout.close();
	delete[] ecoc;
	std::cout << "Successfully wrote ECOC matrix for " << (std::string) argv[1] << std::endl;
	return 0;
}

int do_ecoc(int* ecocOut, int K, int L, double neg, double pos) {
	randomDraws rd(0);

	bool* diff = new bool[L];
	bool* negDiff = new bool[L];
	int* ecoc = new int[K*L];

	int bestHamming = -1;
	bool updated = false;

	int i=0;
	while (i<MAX_ITER) {
		//make a random matrix
		int l=0;
		int stride = 0;
		std::memset(diff, 0, L*sizeof(bool));
		std::memset(negDiff, 0, L*sizeof(bool));

		while (l<L) {
			bool anyNeg = false;
			bool anyPos = false;
			bool notDiff = true;
			bool notNegDiff = true;

			for (int k=0;k<K;k++) {
				double r = rd.getR();
				if (r<neg) {
					anyNeg = true;
					ecoc[stride+k] = -1;
				} else if (1-r<pos) {
					anyPos = true;
					ecoc[stride+k] = 1;
				} else {
					ecoc[stride+k] = 0;
				}
				//validation against previous columns
				if (notDiff) {
					notDiff = false; // this will stick if this vector become different from all the others
					for (int ll=0;ll<l;ll++) { //check for uniqueness of column
						if (!diff[ll]) {
							if (ecoc[stride+k] == ecoc[ll*K+k]) { //we still haven't differentiated from this earlier learner
								notDiff = true;
							} else {
								diff[ll] = true;
							}
						}
					}
				}
				if (notNegDiff) { //we also care if the learners are just negative copies of another
					notNegDiff = false;
					for (int ll=0;ll<l;ll++) { //check for uniqueness of column
						if (!negDiff[ll]) {
							if (ecoc[stride+k] == -ecoc[ll*K+k]) { //we still haven't differentiated from this earlier learner
								notNegDiff = true;
							} else {
								negDiff[ll] = true;
							}
						}
					}	
				}
			}
			//advance if all validation passes
			if (anyPos && anyNeg && !notDiff) {
				l++;
				stride = l*K;	
			} else {
				// std::cout << "bad column" << std::endl;
			}
		}
		// std::cout << "built matrix" << std::endl;

		//further validation now that we have all the columns
		bool bad = false;

		//all rows must be unique
		bad = unique_rows(ecoc, K, L);
		if (bad) {
			// std::cout << "encountered matrix with duplicate rows" <<std::endl;
			continue;
		}
		bad = traversible(ecoc, K, L);
		if (bad) {
			// std::cout << "encountered non-traversible matrix" << std::endl;
			continue;
		}

		// int hamming = INT_MAX;
		int allHamming = INT_MAX;
		for (int k=0;k<K;k++) {
			for (int kk=0; kk<k; kk++) {
				int hamming=0;
				for (int l=0; l<L; l++) {
					hamming += (ecoc[l*K + k]!=0)*(ecoc[l*K + kk]!=0)*(ecoc[l*K + k]!=ecoc[l*K + kk]);
				}
				allHamming = std::min(allHamming, hamming);
			}
		}

		if (allHamming>bestHamming) {
			//save the results
			updated = true;
			bestHamming = allHamming;
			std::memcpy(ecocOut, ecoc, K*L*sizeof(int));
		}

		i++;
	}
	delete[] diff;
	delete[] negDiff;
	delete[] ecoc;

	if (updated) {
		std::cout << "Best Hamming Distance (higher is better): " << bestHamming << std::endl;
		return 0;
	} else {
		return 1;
	}

}

bool unique_rows(int* ecoc, int K, int L) {
	for (int k=0; k<K; k++) {
		for (int kk=0; kk<k; kk++) {
			bool diff = false;
			for (int l=0; l<L; l++) {
				if (ecoc[l*K + k] != ecoc[l*K + kk]) {
					diff = true;
					break;
				}
			}
			if (!diff) {
				return true;
			}
		}
	}

	//if we make it to the end then each row is different from each other row
	return false;
}

bool traversible(int* ecoc, int K, int L) {
	//test whether we can classify results by traversing the matrix

	//rules:
		//can move along K from (-1 to 1) or (1 to -1), only once per L
		//can move along L from (-1 or 1) to (-1 or 1)
		//must be able to move from any one row to any other
	bool* visited = new bool[L];
	bool* traversed = new bool[K];

	for (int k=0; k<K; k++) {
		for (int kk=0;kk<k;kk++) {
			bool sep = false;

			//in the trivial case, there exists a single column where the two classes are nonzero and different
			for (int l=0;l<L;l++) {
				if (ecoc[l*K + k]!=0 && ecoc[l*K+kk]!=0 && ecoc[l*K + k]!=ecoc[l*K+kk]) {
					sep = true;
					break;
				}
			}

			//otherwise it's a maze traversal problem
			std::memset(visited, 0, L*sizeof(bool));
			std::memset(traversed, 0, K*sizeof(bool));
			int unvisited = L;
			int l=0;
			traversed[k]=true;
			traversed[kk]=true;
			while (!sep && l<L) {
				sep = visit(l, k, kk, ecoc, visited, traversed, K, L, unvisited);
				if (sep) {
					break;
				}

				l++;
			}

			if (!sep) { //we couldn't separate these rows, so the matrix is bad
				delete[] visited;
				delete[] traversed;
				return true;
			}
		}
	}
	delete[] visited;
	delete[] traversed;
	return false;
}

bool visit(int l, int k1, int k2, int* ecoc, bool* visited, bool* traversed, int K, int L, int unvisited) {
	//k1 and k2 are the target rows
	//l is the current column
	//visited indicates the legal columns

	// bool success = false;
	if (ecoc[l*K+k1]) {
		if (ecoc[l*K+k2] && ecoc[l*K+k1] != ecoc[l*K+k2]) {
			return true; //we've reached the end
		} else {
			if (unvisited<=1) {
				return false; //there are no more columns left to visit
			}
			unvisited--;
			visited[l] = true;
			
			//try all the vertical moves at this column
			for (int k=0;k<K;k++) {
				if (!traversed[k] && ecoc[l*K+k]!=0 && ecoc[l*K+k]!=ecoc[l*K+k1]) { //we can make this vertical move
					//we can solve the matrix if we can migrate between k2 and k, via traversing this column
					//let's try all the other columns
					traversed[k] = true; //we don't need to revisit this row on this route... otherwise would end up in inifinite loop
					for (int ll=0;ll<L;ll++) {
						if (!visited[ll] && ecoc[ll*K+k]) {
							bool success = visit(ll, k, k2, ecoc, visited, traversed, K, L, unvisited);
							if (success) { //we managed to link these rows
								return true;
							}
						} 
					}
					traversed[k] = false;
				}
			}
			// we failed to return, indicating that there is no viable path
			visited[l] = false;
			return false; //we couldn't find a route
		}
	} else if (ecoc[l*K+k2]) {
		//we can just swap k1 and k2 and try again
		return visit(l, k2, k1, ecoc, visited, traversed, K, L, unvisited);
	} else {
		return false; //we can't move vertically in this column 
	}
}