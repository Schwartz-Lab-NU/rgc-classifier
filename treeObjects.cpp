#include "treeObjects.h"
#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm> //upper_bound, sort...
#include <fstream>
#include <string>
#include <iterator>
#include <cblas.h>
#include <math.h>
#include <iomanip>
#include <filesystem>

//loading psth files
psthSet::psthSet(int N):N(N) {
	data = new psth*[N]; //need to alloc the individual psth's...
	light = true;//the child data are just pointers, we malloced on load
}

psthSet::psthSet(std::ifstream& datafile) {

	//read the number of cells
	datafile.read(reinterpret_cast<char*>(&N), sizeof(unsigned short int));

	// N = 100; //just to speed things up a bit

	datafile.read(reinterpret_cast<char*>(&nt), sizeof(unsigned short int));
	datafile.read(reinterpret_cast<char*>(&ns), sizeof(unsigned short int));
	
	//read the time-dimension metadata
	datafile.read(reinterpret_cast<char*>(&t0), sizeof(double));
	datafile.read(reinterpret_cast<char*>(&tf), sizeof(double));
	datafile.read(reinterpret_cast<char*>(&dt), sizeof(double));

	//read the size-dimension metadata
	datafile.read(reinterpret_cast<char*>(&s0), sizeof(unsigned short int));
	datafile.read(reinterpret_cast<char*>(&sf), sizeof(unsigned short int));
	
	data = new psth*[N];
	// weights = new double[N];
	
	ne = nt * ns;
	double a = 1./N;

	unsigned char cellNameLength;
	for (int i=0; i<N; i++) {
		data[i] = new psth;
		
		//read the length of the cell name
		datafile.read(reinterpret_cast<char*>(&(cellNameLength)), sizeof(unsigned char)); //get the size of the cell name
		data[i]->cellName.resize(cellNameLength); //resize the cell name to match
		
		//read the cell name
		datafile.read(&(data[i]->cellName[0]), cellNameLength*sizeof(unsigned char));

		//read the cell label (on [1,255])
		datafile.read(reinterpret_cast<char*>(&(data[i]->label)), sizeof(char));

		//allocate memory for data
		data[i]->data = new double[ne];
		data[i]->qualityT = new char[nt];
		data[i]->qualityS = new unsigned int[ns];

		//read the quality arrays
		datafile.read(&(data[i]->qualityT[0]), nt*sizeof(char));
		datafile.read(reinterpret_cast<char*>(&(data[i]->qualityS[0])), ns*sizeof(unsigned int));

		//read the data image
		datafile.read(reinterpret_cast<char*>(&(data[i]->data[0])), ne*sizeof(double));

		// weights[i] = a;
		
	}

	qualityT = new unsigned int[nt]();
	qualityS = new unsigned long long int[ns]();

	pdfT = new double[nt];
	pdfS = new double[ns];
	cdfT = new double[nt];
	cdfS = new double[ns];


	stats();
	light = false;
}

psthSet::psthSet(const psthSet* indata) {
	//deep copy of another psthSet
	N = indata->N;
	npos = indata->npos;
	dt = indata->dt;
	t0 = indata->t0;
	tf = indata->tf;
	nt = indata->nt;
	ns = indata->ns;
	s0 = indata->s0;
	sf = indata->sf;
	ne = indata->ne;
	// light = indata->light;
	light = true;//the child data are just pointers, we malloced on load

	data = new psth*[N];
	// weights = new double[N];

	double a = 1./N;

	for (int i=0; i<N; i++) {
		data[i] = indata->data[i]; //just copy the pointer... may slow down operation a bit but reduces memory requirement
		// data[i] = new psth(*indata->data[i], nt, ns);
		// weights[i] = a;
	}

	qualityT = new unsigned int[nt]();
	qualityS = new unsigned long long int[ns]();

	pdfT = new double[nt];
	pdfS = new double[ns];
	cdfT = new double[nt];
	cdfS = new double[ns];

	stats();
}

void psthSet::stats() {
	unsigned int tsum = 0;
	
	for (int i=0;i<N;i++) {
		//update set-level quality
		for (int j=0;j<nt;j++) {
			qualityT[j] += (unsigned int) data[i]->qualityT[j];
			tsum += (unsigned int) data[i]->qualityT[j];
		}
		//update set-level quality
		for (int j=0;j<ns;j++) {
			qualityS[j] += (unsigned long long int) data[i]->qualityS[j];
			// ssum += (unsigned long long int) data->data[i]->qualityS[j];
		}
	}

	//update the pdfs

	double a = 1./N;
	double a2 = 1./tsum;

	pdfT[0] = qualityT[0] * a2; //for t we are just collecting a simple pdf
	cdfT[0] = pdfT[0];
	for (int j=1;j<nt;j++) {
		pdfT[j] = qualityT[j] * a2;
		cdfT[j] = cdfT[j-1] + pdfT[j];
	}
	cdfT[nt-1] = 1.;

	//we need to make multiple passes through s
	double maxS = 0;
	double ssum = 0;

	pdfS[0] = std::max(log(qualityS[0] * a),0.);
	maxS = std::max(pdfS[0], maxS);
	for (int j=1; j<ns; j++) {
		pdfS[j] = std::max(log(qualityS[j] * a),0.);
		maxS = std::max(pdfS[j], maxS);
		ssum += pdfS[j];
	}
	ssum = 1./(maxS*ns - ssum); //this will be the final denominator


	pdfS[0] = (maxS - pdfS[0]) * ssum;
	cdfS[0] = pdfS[0];
	for (int j=1;j<ns; j++) {
		pdfS[j] = (maxS - pdfS[j]) * ssum;
		cdfS[j] = cdfS[j-1] + pdfS[j];
	}
	cdfS[ns-1] = 1.;
}

void psthSet::copy(psthSet* data) {
	ns = data->ns;
	nt = data->nt;
	ne = data->ne;

	qualityT = new unsigned int[nt]();
	qualityS = new unsigned long long int[ns]();
	
	pdfT = new double[nt];
	pdfS = new double[ns];
	cdfT = new double[nt];
	cdfS = new double[ns];
}

psthSet::~psthSet() {
	if (!light) {
		// delete[] weights;
		for (int i=0; i<N; i++) {
			delete data[i];
		}
	} 
	delete[] data;
	delete[] qualityT;
	delete[] qualityS;
	delete[] pdfT;
	delete[] pdfS;
	delete[] cdfT;
	delete[] cdfS;
}

// psth::psth(const psth& indata, int nt, int ns) {
// 	cellName = indata.cellName;
// 	data = new double[nt*ns];
// 	qualityT = new char[nt];
// 	qualityS = new unsigned int[ns];

// 	int stride;
// 	for (int i=0; i<nt; i++) {
// 		qualityT[i] = indata.qualityT[i];
// 		stride = ns*i;
// 		for (int j=0; j<ns; j++) {
// 			data[stride+j] = indata.data[stride+j];
// 		}
// 	}

// 	for (int j=0; j<ns; j++) {
// 		qualityS[j] = indata.qualityS[j];
// 	}
// }

psth::~psth() {
	delete[] data;
	delete[] qualityS;
	delete[] qualityT;
}

//training dependencies
randomDraws::randomDraws(int thread):mt{std::random_device{}() + thread} { //initialize the twister with a random seed
}

double randomDraws::getR() { //draw from the generator and update the seed
	return dist(mt);
}

void randomDraws::getP(int* ind, int N) {
	std::shuffle(ind, ind+N, mt); //shuffles N-many elements at pointer ind
}

featureList::featureList(int F_):F(F_) {
	feats = new featureItem[F];
}

void featureList::get(int N, psthSet* data, randomDraws* r, double* sampleWeights) {

	for (int i=0;i<N;i++) {
		
		while (true) { //loop until nonzero std
			j = r->getR();
			indT = 0;
			while (data->cdfT[indT]<j) {
				indT++;
			}

			j = r->getR();
			indS = 0;
			while (data->cdfS[indS]<j) {
				indS++;
			}
			

			ind = indS*data->nt + indT;
			
			//now we want the mean and std for this feature...

			for (int k=0;k<data->N;k++) {
				feats[i].mean += data->data[k]->data[ind] * sampleWeights[k];
			}
			
			for (int k=0;k<data->N;k++) {
				delta = data->data[k]->data[ind] - feats[i].mean;
				feats[i].istd += delta * delta * sampleWeights[k];
			}
			
			if (feats[i].istd> 0) { //we will keep these indices and populate the featureItem
				
				feats[i].t = indT;
				feats[i].qualityT = data->pdfT[indT] * data->nt;


				feats[i].s = indS;
				feats[i].qualityS = data->pdfS[indS] * data->ns;

				feats[i].quality = feats[i].qualityS * feats[i].qualityS;
				feats[i].istd = feats[i].quality / std::sqrt(feats[i].istd);

				break;
			} else {
				feats[i].mean = 0.;
				feats[i].istd = 0.;
			}

		}
	}
}

featureList::~featureList() {
	delete[] feats;
}

Tree::Tree() {
}

void Tree::init() {
	light = true;

	//inherit properties from ensemble
	inds = ensemble->inds;
	feat = ensemble->feat;
	pred = ensemble->pred;
	bestPred = ensemble->bestPred;
	error = ensemble->error;

	scoreP = ensemble->scoreP;
	scoreN = ensemble->scoreN;
	allPos = ensemble->allPos;
	allNeg = ensemble->allNeg;
	
	labels = ensemble->labels;
	sampleWeights = ensemble->sampleWeights;
	featureWeights = ensemble->featureWeights;

	reg = ensemble->reg;
	rd = ensemble->rd;
	features = ensemble->features;

	treeWeight = 0.;

	minSize = params->minSize;
}
	
int Tree::train(psthSet* data) {

	N = data->N;
	
	//make a copy of the current RNG state
	rd_original = new randomDraws;
	*rd_original = *rd;


	decisionNode* n; //we will train as decision nodes
	root = new decisionNode(N, this, 0); //don't want this going out of scope!!
	alloced=true;

	// n->start = 0; //the root node will see all the tree elements
	// n->tree = this;
	root->depth = 0;

	int err=1;
	int nerr = 0;
	while (err) {
		err = root->train(data, 1);
		if (err) {
			if (nerr<5) {
			std::cout << "Error training root node (bad feature pick, learner "<< ensemble->current<<"). Trying again." << std::endl;
			nerr++;
		} else {
			std::cout << "Error training root node (bad feature pick, learner "<< ensemble->current<<"). Lowering minSize and trying again." << std::endl;
			minSize = params->minSize - nerr + 4;
			nerr++;
			if (minSize<1) {
				delete root;
				root = new Leaf(0);
				root->train(data, 1);
				// return 1;
			}
		}
		}
	}
	// n->complete = true;
	accuracy = root->accuracy;

	// root = n; //this construct (virtualization) allows us to store a pointer to child nodes regardless of what node type they are
	return 0;
}

void Tree::save(std::ofstream& datafile) {
	datafile.write(reinterpret_cast<char*>(&N), sizeof(int));
	datafile.write(reinterpret_cast<char*>(params), sizeof(paramSet));
	datafile.write(reinterpret_cast<char*>(rd_original), sizeof(randomDraws));
	datafile.write(reinterpret_cast<char*>(rd), sizeof(randomDraws));
	datafile.write(reinterpret_cast<char*>(sampleWeights), N*sizeof(double));
	datafile.write(reinterpret_cast<char*>(&treeWeight), sizeof(double));
	root->save(datafile);
}

Tree::Tree(std::ifstream& datafile) {
	light = false;
	//load params

	datafile.read(reinterpret_cast<char*>(&N), sizeof(int));

	params = new paramSet;
	datafile.read(reinterpret_cast<char*>(params), sizeof(paramSet));

	rd_original = new randomDraws;
	datafile.read(reinterpret_cast<char*>(rd_original), sizeof(randomDraws));

	rd = new randomDraws;
	datafile.read(reinterpret_cast<char*>(rd), sizeof(randomDraws));

	sampleWeights = new double[N];
	datafile.read(reinterpret_cast<char*>(sampleWeights), N*sizeof(double));

	datafile.read(reinterpret_cast<char*>(&treeWeight), sizeof(double));

	//load tree
	if (datafile.get()==1) { //we have a decision node
		decisionNode* n;
		n = new decisionNode;
		root = n;
		root-> tree = this;
		root->load(datafile);
		N = root->N;
	} //could be a compact node...
}

void Tree::print(bool o, bool v) {
	if (v) {
		params->print();
	}

	if (o) { // print bottom up
		std::cout << "Printing decision tree (bottom-up, left-to-right): " << std::endl;
		for (int i = params->maxDepth-1; i>0; i--) {
			root->print(i, 0, v);	
		}
	} else { //print top down
		std::cout << "Printing decision tree (top-down, left-to-right): " << std::endl;
		for (int i = 0; i<params->maxDepth; i++) {
			root->print(i, 0,v);	
		}
	}
	std::cout << std::endl;
}

Tree::~Tree() {
	if (alloced) {
		delete root;
		delete rd_original;
	}
	if (!light) { //the following were allocated within the tree so should be deleted here
		delete params;
		delete rd;
	}
}

// void Node::train(psthSet* data, paramSet* params, randomDraws* rd, regression* reg) {
//}

int Node::train(psthSet* data, double sumWeights) {
	return 0;	
}

void Node::save(std::ofstream& datafile) {
}

void Node::load(std::ifstream& datafile) {
}

void Node::print(int d, int o, bool v) {
}

int decisionNode::train(psthSet* data, double sumWeights) {

	featureList* allFeat = tree->features;
	paramSet* params = tree->params;
	randomDraws* rd = tree->rd;
	regression* reg = tree->reg;
	

	int* inds = &(tree->inds[start]);
	double* feat = &(tree->feat[start*params->nFeat]); //do we really need to refer to start here?
	double* pred = &(tree->pred[start]);
	double* bestPred = &(tree->bestPred[start]);
	bool* error = &(tree->error[start]);

	double* scoreP = &(tree->scoreP[start]);
	double* scoreN = &(tree->scoreN[start]);
	double* allPos = &(tree->allPos[start]); //lol...
	double* allNeg = &(tree->allNeg[start]);

	char* labels = &(tree->labels[start]);
	double* sampleWeights = &(tree->sampleWeights[start]);

	int nleft, nright, pleft, pright;
	double wleft, wright;
	
	uncertainty u;
	u.set_params(labels, sampleWeights, N, tree->minSize);
	classU = INFINITY;

	//double posFrac = 0;
	double incr = 1./N;
	int ind;
	int strideI, strideJ, stride;
	bool cv;
	// int nnz;

	if (data->npos > 1 && data->npos < N-1) {
		cv = true;
		allFeat->get(params->nFeat*params->nRepeats, data, rd, sampleWeights);		
	} else { //we don't have enough samples to regress properly, so instead just look at single features
		cv=false;
		allFeat->get(params->nRepeats, data, rd, sampleWeights);
		//this should be easy-ish since there's only 1 sample...
		//just find the feature that puts it the farthest to the end
		lambda = NAN;
		deviance = NAN;
		nnz=1;
		beta1[0] = 1.;
	}
	for (int i=0;i<params->nRepeats;i++) {

		if (cv) {
			strideI = params->nRepeats * i;

			//pull the data up into a local array
			for (int j=0;j<N;j++) {
				strideJ = params->nFeat * j;
				for (int k=0;k<params->nFeat;k++) {
					//normalize to 0 mean, unit sdev
					ind = allFeat->feats[strideI + k].s * data->nt + allFeat->feats[strideI + k].t;
					feat[k + strideJ] =  (data->data[j]->data[ind] - allFeat->feats[strideI+k].mean) * allFeat->feats[strideI+k].istd;
				}
			}
			//
			for (int j=0; j<params->nFeat; j++) {
				tree->featureWeights[j] = allFeat->feats[strideI + j].quality;
			}
			//do logistic regression
			reg->do_fit(feat, labels, sampleWeights, tree->featureWeights, N, data->npos, rd);

			//the regression fit contains the beta parameters as well as the mean and standard deviation (weighted)
			//we will apply these to the data samples to combine the predictors
			for (int j=0;j<N;j++) {
				pred[j] = 0.;
				strideJ = params->nFeat * j;
				//stride = strideI + j;
				for (int k=0;k<params->nFeat;k++) {
					ind = k+strideJ;
					pred[j] += feat[ind] * reg->beta1[k]; //note that we normalized using the same pointers!!!
				}
			}
		} else {
			for (int j=0;j<N;j++) {
				ind = allFeat->feats[i].s * data->nt + allFeat->feats[i].t;
				pred[j] =  (data->data[j]->data[ind] - allFeat->feats[i].mean) * allFeat->feats[i].istd;
			}
		}

		//with the combined predictors we will pick the threshold using weighted entropy
		u.get_uncertainty(pred, labels, sampleWeights, inds);

		if (u.r == i && !isinf(u.minU)) { //this iteration is the best so far, so store the results
			nleft = u.nleft;
			pleft = u.pleft;
			wleft = u.wleft;
			nright = u.nright;
			pright = u.pright;
			wright = u.wright;

			if (cv) {
				lambda = reg->lambda;
				deviance = reg->deviance;
				classU =u.minU;
				support = u.support;
				beta0 = u.bias;

				ind=0;
				for (int j=0;j<params->nFeat;j++) {
					if (reg->beta1[j] != 0) {
						t[ind] = allFeat->feats[strideI + j].t;
						s[ind] = allFeat->feats[strideI + j].s;
						qualityT[ind] = allFeat->feats[strideI + j].qualityT;
						qualityS[ind] = allFeat->feats[strideI + j].qualityS;
						quality[ind] =  allFeat->feats[strideI + j].quality;

						beta1[ind] = reg->beta1[j];
						mean[ind] = allFeat->feats[strideI+j].mean;
						istd[ind] = allFeat->feats[strideI+j].istd;

						ind++;
					}
				}
				nnz = ind;
				for (int j=ind;j<params->nFeat;j++) {
					t[j] = NAN;
					s[j] = NAN;
					qualityT[j] = NAN;
					qualityS[j] = NAN;
					quality[j] =  NAN;

					beta1[j] = NAN;
					mean[j] = NAN;
					istd[j] = NAN;

				}
			} else {
				t[0] = allFeat->feats[i].t;
				s[0] = allFeat->feats[i].s;
				qualityT[0] = allFeat->feats[i].qualityT;
				qualityS[0] = allFeat->feats[i].qualityS;
				quality[0] =  allFeat->feats[i].quality;
				mean[0] = allFeat->feats[i].mean;
				istd[0] = allFeat->feats[i].istd;

				ind++;
			}
			for (int j=0;j<N; j++) {
				// std::cout << pred[j] << ",";
				bestPred[j] = pred[j];
			}
			// std::cout << std::endl;
		}
	}
	if (isinf(classU)) {
		// std::cout <<"error (depth, N, npos, isRight): "<< depth << "," << N <<"," << data->npos << "," << leftRight <<std::endl;
		// std::cout << u.minU << std::endl;
		return 1; //error, convert to leaf
	}	

	//sort the indices by the regression result for the desired r,c...
	for (int i=0;i<N;i++) {
		inds[i] = i; //initialize sorting...
	}
	// stride = u.r*N;
	std::stable_sort(inds, &(inds[N]), [&bestPred](size_t a, size_t b) {return bestPred[a] < bestPred[b];});
	
	//int last = data->N;
	// ind = INT_MAX;
	psthSet subData(std::max(nleft, nright));
	subData.copy(data);
	subData.N = nleft;
	subData.npos = pleft;

	// sort the data so that the first nleft-many elements are at the beginning
	for (int i=0;i<nleft;i++) {
		// std::cout << bestPred[i] << ",";
		while (inds[i] != i) {
			//need to swap data[inds[i]] with data[ind[inds[i]]]!!!
			std::swap(data->data[inds[i]], data->data[inds[inds[i]]]); //just swaps the pointers to each cell
			std::swap(labels[inds[i]],labels[inds[inds[i]]]);
			std::swap(sampleWeights[inds[i]],sampleWeights[inds[inds[i]]]);
			
			std::swap(error[inds[i]], error[inds[inds[i]]]);
			std::swap(scoreP[inds[i]], scoreP[inds[inds[i]]]);
			std::swap(scoreN[inds[i]], scoreN[inds[inds[i]]]);	

			std::swap(allPos[inds[i]], allPos[inds[inds[i]]]);	
			std::swap(allNeg[inds[i]], allNeg[inds[inds[i]]]);	

			std::swap(inds[i], inds[inds[i]]);
		}
	} //once the left part of the tree is done the order of the right tree doesn't matter

	for (int i=0; i<nleft;i++) {
		subData.data[i] = data->data[i];
		sampleWeights[i] *= wleft; //renormalize the weights
	}

	//calculate the pdf, etc. for the subset...
	//train the left half of the subtree
	Leaf *ll, *rl;
	decisionNode *ln, *rn;
	int err;

	// subData.weights = data->weights;
	// std::cout << "Left: "<< nleft<<","<< depth << "," << pleft << "," << std::endl;
	if (depth < params->maxDepth-2 && 2*pleft-1>=tree->minSize && 2*(nleft-pleft)-1>=tree->minSize) { //make a decision node

		subData.stats();
		
		left = new decisionNode(nleft, tree, start);
		// left->start = start;
		
		left->parent = this;
		// ln->tree = tree;
		left->leftRight = false;
		left->depth = depth+1;

		err = left->train(&subData, sumWeights/wleft);

		if (err) {
			delete left;
		}
	} else {
		err = true;
	}

	if (err) {//make a new leaf
		left = new Leaf(start);
		// left->start = start;
		left->parent = this;
		left->tree = tree;
		left->leftRight = false;
		left->depth = depth+1;

		left->train(&subData, sumWeights/wleft);

		// left = ll;
		// left->print(depth+1, 0, false);
	}

	//unwind the rescaling
	wleft = 1./wleft;
	for (int i=0; i<nleft;i++) {
		sampleWeights[i] *= wleft; 
	}

	//train the right half of the subtree
	//note: we can reliably pass down the data here because the left half ignores the right
	subData.N = nright;
	subData.npos = pright;

	for (int i = nleft;i<N;i++) {
		subData.data[i-nleft] = data->data[i];
		sampleWeights[i] *= wright; //renormalize the weights
		
	}
	// subData.weights = &(data->weights[nleft]);

	// std::cout << "Right: "<< nright<<","<< depth << "," << pright << "," << std::endl;
	if (depth < params->maxDepth-2 && 2*pright-1>=tree->minSize && 2*(nright-pright)-1>=tree->minSize) { //make a decision node
		
		subData.stats();
		
		right = new decisionNode(nright, tree, start+nleft);
		// rn->start = start + nleft;
		
		right->parent = this;
		// rn->tree = tree;
		right->leftRight = true;
		right->depth = depth+1;

		err = right->train(&subData, sumWeights/wright);
		
		if (err) {
			delete right;
		}
	} else {
		err = true;
	}

	if (err) {//make a new leaf
		right = new Leaf(start+nleft);
		// rl->start = start + nleft;

		right->parent = this;
		right->tree = tree;
		right->leftRight = true;
		right->depth = depth+1;
	
		right->train(&subData, sumWeights/wright);

		// right = rl;
		// right->print(depth+1, 0, false);
	}

	//unwind the rescalintg
	wright = 1./wright;
	for (int i = nleft;i<N;i++) {
		sampleWeights[i] *= wright; //renormalize the weights
	}
	
	//let the accuracy and misclassificaiton count trickle up
	accuracy = (left->accuracy*wleft + right->accuracy*wright);
	misclassified = left->misclassified + right->misclassified;

	// std::cout << "total weight at node: " << sumWeights << ", acc (weighted) " << accuracy<<", acc (raw): " <<  1 - (double)misclassified/N <<std::endl;
	//and we're done!
	complete = true;
	return 0;
}

void decisionNode::save(std::ofstream& datafile) {
	//first save a char 1 indicating that this is a decisionNode
	datafile.put(1);

	//now save the node parameters
	datafile.write(reinterpret_cast<char*>(&N), sizeof(unsigned int));
	datafile.write(reinterpret_cast<char*>(&nnz), sizeof(unsigned int));
	
	datafile.write(reinterpret_cast<char*>(&(t[0])), nnz*sizeof(double));
	datafile.write(reinterpret_cast<char*>(&(s[0])), nnz*sizeof(double));
	datafile.write(reinterpret_cast<char*>(&(beta1[0])), nnz*sizeof(double));
	datafile.write(reinterpret_cast<char*>(&(mean[0])), nnz*sizeof(double));
	datafile.write(reinterpret_cast<char*>(&(istd[0])), nnz*sizeof(double));
	
	datafile.write(reinterpret_cast<char*>(&(qualityT[0])), nnz*sizeof(double));
	datafile.write(reinterpret_cast<char*>(&(qualityS[0])), nnz*sizeof(double));
	datafile.write(reinterpret_cast<char*>(&(quality[0])), nnz*sizeof(double));
	
	datafile.write(reinterpret_cast<char*>(&beta0), sizeof(double));
	datafile.write(reinterpret_cast<char*>(&lambda), sizeof(double));
	datafile.write(reinterpret_cast<char*>(&deviance), sizeof(double));
	datafile.write(reinterpret_cast<char*>(&classU), sizeof(double));
	datafile.write(reinterpret_cast<char*>(&support), sizeof(double));
	
	//save the daughter nodes 
	left->save(datafile);
	right->save(datafile);
}

void decisionNode::load(std::ifstream& datafile) {
	complete = true;

	datafile.read(reinterpret_cast<char*>(&N), sizeof(unsigned int));
	datafile.read(reinterpret_cast<char*>(&nnz), sizeof(unsigned int));

	t = new double[nnz];
	s = new double[nnz];
	qualityT = new double[nnz];
	qualityS = new double[nnz];
	quality = new double[nnz];
	beta1 = new double[nnz];
	mean = new double[nnz];
	istd = new double[nnz];

	datafile.read(reinterpret_cast<char*>(&(t[0])), nnz*sizeof(double));
	datafile.read(reinterpret_cast<char*>(&(s[0])), nnz*sizeof(double));
	datafile.read(reinterpret_cast<char*>(&(beta1[0])), nnz*sizeof(double));
	datafile.read(reinterpret_cast<char*>(&(mean[0])), nnz*sizeof(double));
	datafile.read(reinterpret_cast<char*>(&(istd[0])), nnz*sizeof(double));
	
	datafile.read(reinterpret_cast<char*>(&(qualityT[0])), nnz*sizeof(double));
	datafile.read(reinterpret_cast<char*>(&(qualityS[0])), nnz*sizeof(double));
	datafile.read(reinterpret_cast<char*>(&(quality[0])), nnz*sizeof(double));

	datafile.read(reinterpret_cast<char*>(&beta0), sizeof(double));
	datafile.read(reinterpret_cast<char*>(&lambda), sizeof(double));
	datafile.read(reinterpret_cast<char*>(&deviance), sizeof(double));
	datafile.read(reinterpret_cast<char*>(&classU), sizeof(double));
	datafile.read(reinterpret_cast<char*>(&support), sizeof(double));

	//load the daughter nodes 
	int nodeType = datafile.get();

	decisionNode *ln;
	Leaf *ll;
	if (nodeType==1) { //we have a decision node
		ln = new decisionNode;
		left = ln;
	} else if (nodeType==2) { // we have a leaf
		ll = new Leaf;
		left = ll;
	}

	left->parent = this;
	left->tree = tree;
	left->leftRight = false;
	left->depth = depth+1;

	left->load(datafile);

	//now move to the right
	nodeType = datafile.get();
	decisionNode *rn;
	Leaf *rl;
	if (nodeType==1) { //we have a decision node
		rn = new decisionNode;
		right = rn;
	} else if (nodeType==2) { // we have a leaf
		rl = new Leaf;
		right = rl;
	}

	right->parent = this;
	right->tree = tree;
	right->leftRight = true;
	right->depth = depth+1;	

	right->load(datafile);

	accuracy = (left->accuracy*left->N + right->accuracy*right->N)/N;
	misclassified = left->misclassified + right->misclassified;
}

void decisionNode::print(int d, int o, bool v) {
	if (d==depth) {
		std::cout << std::showpos << std::fixed << std::setprecision(2);
		std::cout << "Node (depth "<< depth <<", parent " << o/2<< ", this "<< o << "): " << N << "cells, uncertainty: " << classU << ", accuracy: " << accuracy << ", misclassified: " << misclassified <<std::endl;
		
		if (v) {
			std::cout << std::setw(18)<< "weights: [" << std::setprecision(4); 
			for (int i =0; i<nnz; i++) {
				std::cout << std::setw(15)<< beta1[i] << ", ";
			}
			std::cout << "]" <<std::endl;

			std::cout << std::setw(18)<< "mean: [";
			for (int i =0; i<nnz; i++) {
				std::cout << std::setw(15)<< mean[i] << ", ";
			}
			std::cout << "]" <<std::endl;

			std::cout << std::setw(13)<<"std" << "\u207b"<< "\u00b9"<<": [";
			for (int i =0; i<nnz; i++) {
				std::cout << std::setw(15) << istd[i] << ", ";
			}
			std::cout << "]" <<std::endl;

			std::cout << std::setw(18) <<"index: [" << std::setprecision(0);
			for (int i =0; i<nnz; i++) {
				std::cout << "(" << std::setw(6) << s[i] << "," << std::setw(6) << t[i] << "), ";
			}
			std::cout << "]" <<std::endl;

			std::cout << std::setw(18)<<"quality: [" <<std::setprecision(2);
			for (int i =0; i<nnz; i++) {
				std::cout << "(" <<std::setw(6)<< qualityS[i] << "," << std::setw(6)<<qualityT[i] << "), ";
			}
			std::cout << "]" <<std::endl <<std::setprecision(4);

			std::cout << std::setw(18)<< "bias:  " << beta0 << std::endl;
			std::cout << std::setw(16)<< "\u03bb" << ":  " << lambda << std::endl;
			std::cout << std::setw(18)<< "deviance:  " << deviance << std::endl;
			std::cout << std::setw(18)<< "uncertainty:  " << classU << std::endl;
			std::cout << std::setw(18)<< "support:  " << support << std::endl;
		}
	} else {
		left->print(d,2*o,v);
		right->print(d,2*o+1,v);
	}
}

decisionNode::decisionNode(int N_, Tree* tree_, int start_) {

	N=N_;
	tree=tree_;
	start = start_;

	//allocate member variables
	t = new double[tree->params->nFeat];
	s = new double[tree->params->nFeat];
	qualityT = new double[tree->params->nFeat];
	qualityS = new double[tree->params->nFeat];
	quality = new double[tree->params->nFeat];
	beta1 = new double[tree->params->nFeat];
	mean = new double[tree->params->nFeat];
	istd = new double[tree->params->nFeat];
}

decisionNode::~decisionNode() {
	// std::cout <<"attempting to delete node" <<std::endl;
	delete[] t;
	delete[] s;
	delete[] qualityT;
	delete[] qualityS;
	delete[] quality;
	delete[] beta1;
	delete[] mean;
	delete[] istd;

	if (complete) {
		delete left;
		delete right;
	}
}

Leaf::Leaf(int start_) {
	start = start_;
}


int Leaf::train(psthSet* data, double sumWeights) {
	N = data->N;

	//can't we just inherit these from the node now?
	pos = 0;
	neg = 0;
	fracPos = 0.;
	fracNeg = 0.;

	char* labels = &(tree->labels[start]);
	bool* error = &(tree->error[start]);
	double* scoreP = &(tree->scoreP[start]);
	double* scoreN = &(tree->scoreN[start]);
	double* sampleWeights = &(tree->sampleWeights[start]);

	accuracy = 0.;

	for (int i=0; i<N; i++) {
		if (labels[i]==1) {
			pos += 1;
			fracPos += sampleWeights[i];
		}
	}
	neg = N - pos;
	fracNeg = (1 - fracPos) * sumWeights;
	fracPos *= sumWeights;

	for (int i=0; i<N; i++) {
		scoreP[i] = fracPos;
		scoreN[i] = fracNeg;
	}

	// fracPos = ((double) pos) / N;
	// fracNeg = ((double) neg) / N;

	if (fracPos > fracNeg) {
		MLE = 1;
		accuracy = fracPos / sumWeights;
		misclassified = neg;

		for (int i=0; i<N; i++) {
			if (labels[i]==1) {
				error[i] = false;
			} else {
				error[i] = true;
			}
		}

	} else if (fracNeg > fracPos) {
		MLE = -1;
		accuracy = fracNeg / sumWeights;
		misclassified = pos;

		for (int i=0; i<N; i++) {
			if (labels[i]==1) {
				error[i] = true;
			} else {
				error[i] = false;
			}
		}

	} else { //pos and neg are equal, so our best bet is to just randomly pick one
		MLE = 0;
		accuracy = .5;
		misclassified = N/2;

		for (int i=0; i<N; i++) {
			error[i] = false; //just repot that everything is wrong
		}
	}

	// std::cout << "total weight at leaf: " << sumWeights << ", acc (weighted): " << accuracy << ", score(+): " << fracPos << std::endl;

	return 0;
}

void Leaf::save(std::ofstream& datafile) {
	//save a char 2 indicating this is a leaf
	datafile.put(2);
	
	//now save the leaf parameters
	datafile.write(reinterpret_cast<char*>(&pos), sizeof(unsigned int));
	datafile.write(reinterpret_cast<char*>(&neg), sizeof(unsigned int));
}

void Leaf::load(std::ifstream& datafile) {

	//now save the leaf parameters
	datafile.read(reinterpret_cast<char*>(&pos), sizeof(unsigned int));
	datafile.read(reinterpret_cast<char*>(&neg), sizeof(unsigned int));

	N = pos + neg;
	fracPos = (double) pos / N;
	fracNeg = 1-fracPos;

	if (pos>neg) {
		MLE = 1;
		misclassified = neg;
		accuracy = fracPos;
	} else if (pos<neg) {
		MLE = -1;
		misclassified = pos;
		accuracy = fracNeg;
	} else {
		MLE = 0;
		misclassified = N/2.;
		accuracy = .5;
	}
}

void Leaf::print(int d, int o, bool v) {
	if(d==depth) {
		std::cout << std::showpos << std::fixed << std::setprecision(2);
		std::cout << "Leaf (depth "<< depth <<", parent " << o/2 << ", this"<< o << "): " << pos << "pos, " << neg << "neg, accuracy: " <<accuracy << ", misclassified: " << misclassified << std::endl;
	}	
}

regression::regression(paramSet* params, int maxN) {
	// note X is nFeat-by-N
	//assumes features are normalized

	F = params->nFeat;
	nFolds = params->nFolds;
	int ne = F*maxN;

	B = new double[F]; //start with 0s
	beta1 = new double[F];

	Bo = new double[F];
	Bi = new double[F]; 
	muX = new double[F];
	deno = new double[F];
	active = new bool[F]; //start false
	constant = new bool[F];
	wX2 = new double[F];
	// mean = new double[F];
	// istd = new double[F];

	//allocate memory -- we're saving more than we need here but shouldn't be an issue
	P = new double[maxN];
	W = new double[maxN];
	//nW = new double[maxN];
	ones = new double[maxN];
	negones = new double[maxN];
	Z = new double[maxN];
	XB = new double[maxN];

	Xc = new double[ne];
	wX = new double[ne];
	X2 = new double[ne];

	for (int i=0;i<maxN;i++) {
		ones[i] = 1.;
		negones[i] = -1.;
	}

	nLambda = params-> nLambda;
	alpha = params -> alpha;

	lastDev = new double[nLambda];

	if (nFolds>0) {

		seDev = new double[nLambda];
		inds = new int[maxN];
		stdX = new double[F];
		lmuX = new double[F];
	

		localX = new double[(maxN - (maxN / nFolds) + 1)*F];
		localY = new char[maxN - (maxN / nFolds) + 1];
		localSW = new double[maxN - (maxN / nFolds) + 1];

	}
	col_major = CblasColMajor;
	y = CblasTrans;
	n = CblasNoTrans;
}

regression::~regression() {
	delete[] B;//start with 0s
	delete[] beta1;

	delete[] Bo;
	delete[] Bi;
	delete[] muX;
	delete[] deno;
	delete[] active; //start false
	delete[] constant;
	delete[] wX2;
	// delete[] mean;
	// delete[] istd;

	delete[] P;
	delete[] W;
	delete[] ones;
	delete[] negones;
	delete[] Z;
	delete[] XB;

	delete[] Xc;
	delete[] wX;
	delete[] X2;

	delete[] lastDev;

	if (nFolds>0) {

		delete[] seDev;
		delete[] inds;
		delete[] stdX;
		delete[] lmuX;
	

		delete[] localX;
		delete[] localY;
		delete[] localSW;

	}
}

void regression::do_fit(double* X, char* Y, double* sw, double* fw, int N, int npos, randomDraws* r) {
	
	//Print the labels and predictors, for copying into matlab...
	// std::cout << std::fixed << "X = [";
	// for (int i=0;i<N;i++) {
	// 	for (int j=0; j<F;j++) {
	// 		std::cout << X[i*F + j] << ",";
	// 	}
	// 	std::cout << ";" << std::endl;
	// }
	// std::cout << "];" <<std::endl;

	// std::cout << std::fixed << "Y = [";
	// for (int i=0;i<N;i++) {
	// 	std::cout << (int) Y[i] << ",";
	// }
	// std::cout << "];" <<std::endl;


	//compute the lambda sequence

	//assume sw is normalized to sum 1

	// cblas_dgemv(col_major, n, F, N, 1., X, F, sw, 1, 0., mean, 1); //update muX, the weighted average
	// for (int i=0;i<N;i++) {
	// 	stride = i*F;
	// 	for (int j=0;j<F;j++) {
	// 		ind = stride+j;
	// 		X[ind] = X[ind] - mean[j];
	// 		X2[ind] = X[ind] * X[ind];
	// 	}
	// }

	// cblas_dgemv(col_major, n, F, N, 1., X2, F, sw, 1, 0., istd, 1);
	// for (int j=0;j<F;j++) {
	// 	istd[j] = fw[j] / std::sqrt(istd[j]);
	// }

	muY = 0;
	for (int i=0;i<N;i++) {
		stride = i*F;
		for (int j=0;j<F;j++) {
			ind = stride+j;
			// X[ind] = X[ind] * istd[j]; //standardize X
			wX[ind] = X[ind] * sw[i];
		}
		muY += Y[i] * sw[i];
	}

	for (int i=0;i<N;i++) {
		Z[i] = Y[i] - muY; //center Y, placeholder
	}

	cblas_dgemv(col_major, n, F, N, 1., wX, F, Z, 1, 0., B, 1); // Yb = X*(Y-muY), placeholder

	lambdaMax = 0;

	for (int i=0;i<F;i++) {
		lambdaMax = std::max( abs(B[i]), lambdaMax );
	}

	lambdaMax = lambdaMax / alpha;
	lambdaMin = log(lambdaMax* .0001);
	lambdaMax = log(lambdaMax); 

	deviance = INFINITY;
	if (nFolds > 0) {

		//randomly permute the samples
		stride=0;
		stride2=0;
		for (int i=0; i<N; i++) {
			if (Y[i]==1) { 
				inds[stride] = i;
				stride++;
			} else {
				inds[stride2+npos] = i;
				stride2++; 
			}
		} //now the first npos elements in inds are the positive indices, followed by the negative indices
		
		//shuffle the positive and negative indices independently
		r->getP(inds,npos);
		r->getP(&(inds[npos]),N-npos);

		//prepare to store feature statistics
		for (int i=0;i<nLambda; i++) {
			lastDev[i] = 0;	
			seDev[i] = 0;
		}

		for (int j=0;j<nFolds; j++) {

			//get this fold
			//localN = N - (N / nFolds + (N % nFolds > j ? 1 : 0));

			//test set:
			//num pos: npos/nFolds + (npos%nFolds > j ? 1 : 0)
			//num neg: (N-npos)/nFolds + ((N-npos) %nFolds > j ? 1 : 0);
			localNp = npos - (npos/nFolds + (npos%nFolds > j ? 1 : 0));
			localNn = (N-npos) - ((N-npos)/nFolds + ((N-npos) %nFolds > j ? 1 : 0));
			localN = localNp + localNn;

			for (int k=0;k<F; k++) {
					// lmuX[k] = 0.;
					// stdX[k] = 0.;
					B[k] = 0.;
					active[k] = false;
			}
			
			//grab the training data
			//muY = 0;
			ind=0;
			sumWi = 0.;
			for (int i=0; i<N; i++) {
				if ((i<npos && i%nFolds == j) || i>=npos && ((i-npos)%nFolds == j)) { //the cell is in the test fold
					continue;
				}
				stride = inds[i]*F;
				stride2 = ind*F;
				localY[ind] = Y[inds[i]]; 

				localSW[ind] = sw[inds[i]]; //need to renormalize!!!
				sumWi += localSW[ind];

				for (int k=0; k<F; k++) {
					localX[stride2+k] = X[stride+k];
				}
				ind++; //increment the index
			}
			sumWi = 1./sumWi;

			//calculate the predictor statistics and normalize the data
			cblas_dgemv(col_major, n, F, localN, sumWi, localX, F, localSW, 1, 0., lmuX, 1); //update muX, the weighted average
			for (int i=0;i<localN;i++) {
				stride = i*F;
				for (int k=0;k<F;k++) {
					ind = stride+k;
					localX[ind] = localX[ind] - lmuX[k];
					X2[ind] = X[ind] * X[ind];
				}
			}
			cblas_dgemv(col_major, n, F, localN, sumWi, X2, F, localSW, 1, 0., stdX, 1);
			for (int k=0;k<F;k++) {
				stdX[k] = std::sqrt(stdX[k]);

				if (abs(stdX[k]) > 1e-10) {
					stdX[k] = fw[k] / stdX[k];
					constant[k] = false;
				} else {//this feature is constant; we don't want to divide by zero or perform updates on it
					stdX[k] = 1;
					constant[k] = true;
				}
			}

			//standardize the predictors
			for (int i=0; i<localN; i++) {
				stride = i*F;
				for (int k=0; k<F; k++) {
					ind = stride + k;
					localX[ind] = localX[ind] * stdX[k];
				}
			
				//populate the internal state with the initial conditions
				W[i] = .1875;

				if (localY[i]==1) {
					P[i] = .75;
					XB[i] = 1.09861228867;
				} else {
					P[i] = .25;
					XB[i] = -1.09861228867;
				}
			}

			// perform irls
			for (int i = nLambda - 1; i >= 0; i--) {
				//log space the lambdas
				thisLambda = exp( lambdaMin + ((double)i)/(nLambda-1) * (lambdaMax - lambdaMin));

				gamma = thisLambda*alpha;
				denoL = thisLambda*(1-alpha);

				irls(localX,localY,localSW,localN); //calculate B

				// std::cout << std::fixed << "B = ["; //print the results
				// for (int j=0; j<F;j++) {
				// 	std::cout << B[j] << ",";
				// }
				// std::cout << "];" <<std::endl;
				
				// unwind the scaling
				for (int k=0;k<F;k++) {
					Bo[k] = B[k] * stdX[k]; //we don't want to change B because we will use it for warm starting the next solution
				}

				intercept -= cblas_ddot(F, lmuX, 1, Bo, 1);
				//okay to change intercept since XB was already calculated internally

				//now we want to apply B to the test set to calculate the deviance
				delta = 0.;
				for (int k=0; k<N; k++) {
					if ((k<npos && k%nFolds != j) || k>=npos && ((k-npos)%nFolds != j)) { //the cell is in the test fold
						continue;
					}
					// if (inds[k]%nFolds != j) { //this sample is in the training fold
					// 	continue;
					// }
					//predict the label for this sample
					lmuY = 1. / (1 + exp( -intercept - cblas_ddot(F, &(X[inds[k]*F]), 1, Bo, 1) ));
					//rectify the prediction to avoid divide by 0
					lmuY = lmuY < pe ? pe : (lmuY > mpe ? mpe : lmuY);

					//accumulate the deviance of the model from perfect solution
					delta += (Y[inds[k]]==1 ? log(1./lmuY) : log(1./(1-lmuY))) * sw[inds[k]];
				}
				lastDev[i] += delta;
				seDev[i] += delta*delta;

			}
		} //endfor (folds)

		sqrtN = 1. /std::sqrt(nFolds);
		for (int i=0; i<nLambda-1; i++) {//note we're avoiding maxLambda since we don't want to return nonsense, even if the model has no predictive power
			
			// std::cout << 2*lastDev[i]/nFolds << "(+/- " << std::sqrt((seDev[i] - (lastDev[i] * lastDev[i]) / nFolds) / (nFolds-1)) * sqrtN << ")," <<std::endl;
			
			if (2*lastDev[i] < deviance) {
				seDev[i] = std::sqrt((seDev[i] - (lastDev[i] * lastDev[i]) / nFolds) / (nFolds-1)) * sqrtN;
				deviance  =	2*lastDev[i];
				devThresh = 2*lastDev[i] / nFolds + seDev[i];
				ind = i;
			} else if (2*lastDev[i] /nFolds <= devThresh) {
				ind = i; //this is the index of min + 1se
			}
			
		}
		//just print the dev corresponding to max Lambda but we don't want to use it
		// std::cout << 2*lastDev[nLambda-1]/nFolds << "(+/- " << std::sqrt((seDev[nLambda-1] - (lastDev[nLambda-1] * lastDev[nLambda-1]) / nFolds) / (nFolds-1)) * sqrtN << ")," <<std::endl;
			
		// std::cout << std::endl << ind << ": " << 2*lastDev[ind]/nFolds << " ( < " << devThresh << ")";
		// std::cout << std::endl <<std::endl;
	} else {
		ind = nLambda - 1;
	}

	//calculate B using the full dataset

	//populate the internal state with the initial conditions
	for (int i=0;i<N;i++) {
		W[i] = .1875;
		if (Y[i]==1) {
			P[i] = .75;
			XB[i] = 1.09861228867;
		} else {
			P[i] = .25;
			XB[i] = -1.09861228867;
		}
	}
	// sumWi = 1./(.1875*N);
	
	//muY = muYin;

	for (int k=0;k<F; k++) {
		//muX[k] = 0.;
		constant[k] = false; //features were picked to guarantee this	
		B[k] = 0.;
		active[k] = false;				
	}

	for (int i = ind; i >= 0; i--) {
		//log-space the lambda values
		thisLambda = exp( lambdaMin + ((double)i)/(nLambda-1) * (lambdaMax - lambdaMin));
		gamma = thisLambda*alpha;
		denoL = thisLambda*(1-alpha);

		//std::cout << "L(" << i << ") = " << thisLambda / N << ";" << std::endl;
		
		irls(X,Y,sw,N); //update the internal state for the given lambda

		// std::cout << std::fixed << "B = ["; //print the results
		// for (int j=0; j<F;j++) {
		// 	std::cout << B[j] << ",";
		// }
		// std::cout << "];" <<std::endl;

		if (nFolds == 0) {
			//calculate the deviance
			lastDev[i] = 0;
			for (int j=0; j<N; j++) {
				lastDev[i] += (Y[j]==1 ? log(1./P[j]) : log(1./(1-P[j])))*sw[j];
			}
			lastDev[i] = 2*lastDev[i];

			if (lastDev[i] < deviance) {
				beta0 = intercept;
				for (int j=0; j<F; j++) {
					beta1[j] = B[j];
				}
				lambda = thisLambda / N;
				deviance = lastDev[i];
			}
		} else {
			//make sure B is not zeros
			activeFlag = false;
			for (int j=0;j<F;j++) {
				if (abs(B[j]) > 1e-7) {
					break;
				}
				if (j == F-1) {
					activeFlag = true;
				}
			}
			if (activeFlag) {
				continue; //lower lambda and try again
			}

			//store B
			for (int j=0;j<F;j++) {
				beta1[j] = B[j];
			}
			beta0 = intercept;
			lambda = thisLambda / N;
			return; //no need to loop through other lambda
		}
		//repeat the loop, using the internal state to warm-start the solution
	}
}

//the brains behind the operation
void regression::irls(double* X, char* Y, double* Ws, int Nlocal) {

	for (int m=0; m<maxIterO; m++) {
		//outer update

		//store B
		cblas_dcopy(F, B, 1, Bo, 1);

		sumWi = 0;
		for (int i=0;i<Nlocal;i++) {
			Z[i] = XB[i] + ((Y[i]-P[i]) / W[i]); //update Z, locally linear

			W[i] = W[i] * Ws[i];
			sumWi += W[i];
		}
		sumWi = 1. /sumWi;

		//W <- W .* Ws
		// sumWi <- sum(W)

		// muZ <- <nW, Z> the weighted average
		muZ = cblas_ddot(Nlocal, W, 1, Z, 1) * sumWi;

		//calculate muX here...
		cblas_dgemv(col_major, n, F, Nlocal, sumWi, X, F, W, 1, 0., muX, 1); //update muX, the weighted average


		//std::cout << sumWi << std::endl;

		// Z <- Z - muZ* (ones) <- the centered response
		cblas_daxpy(Nlocal, muZ, negones, 1, Z, 1);

		// update X-dependencies...
		for (int i=0; i<Nlocal; i++) {
			stride = i*F;
			for (int j=0; j<F; j++) {
				ind = j + stride;
				Xc[ind] = X[ind] - muX[j]; //the centered predictors
				wX[ind] = Xc[ind] * W[i]; //the weighted centered predictors
				X2[ind] = Xc[ind] * Xc[ind]; //squared centered predictors
			}
		}

		cblas_dgemv(col_major, y, F, Nlocal, 1., Xc, F, B, 1, 0., XB, 1); // XB <- X0*B (the centered model)
		cblas_daxpy(Nlocal, -1, XB, 1, Z, 1);// Z <- Z0 - X0*B (the model residual)
		cblas_dgemv(col_major, n, F, Nlocal, 1., X2, F, W, 1, 0., wX2, 1); //wX2 = X2 * W 

		//calculate the denominator
		for (int i=0; i<F; i++) {
			deno[i] = 1./(wX2[i] + denoL);
		}

		//inner update -- Z is linear on X
		for (int k=0; k<maxIterI; k++) {
			maxD = 0;
			cblas_dcopy(F, B, 1, Bi, 1);

			//active set loop
			for (int i=0; i<F; i++) {
				if (constant[i] || !active[i]) {
					continue;
				}
				bp = cblas_ddot(Nlocal, Z, 1, &(wX[i]), F) + (B[i]*wX2[i]); //regress the residual onto this predictor

				margin = abs(bp) - gamma;

				if (margin <= 0) { //lasso penalty -> B goes to 0 and the predictor leaves the active set
					active[i] = false;
					delta = B[i];
					B[i] = 0;
				} else {
					B[i] = margin * deno[i] * (bp > 0 ? 1 : -1);
					delta = Bi[i] - B[i];
				}
				maxD = std::max(maxD, abs(delta));

				//Z = Z - X0*delta
				cblas_daxpy(Nlocal, delta, &(Xc[i]), F, Z, 1);
			}

			if (maxD < tol) { //the active set has converged
				activeFlag = true;

				//inactive set loop
				for (int i=0; i<F; i++) {
					if (constant[i] || active[i]) {
						continue;
					}

					bp = cblas_ddot(Nlocal, Z, 1, &(wX[i]), F); //regress the residual onto this predictor
					margin = abs(bp) - gamma;

					if (margin > 0) { // due to updates in the active set, the predictor should re-enter the model 
						active[i] = true;
						activeFlag = false;
						B[i] = margin * deno[i] * (bp > 0 ? 1 : -1);
						maxD = std::max(maxD, abs(B[i]));

						//recall that Bi[i] === 0, since it was inactive
						// Z = Z - X0*B
						cblas_daxpy(Nlocal, -B[i], &(Xc[i]), F, Z, 1);
					}
					
				}
				if (activeFlag || maxD < tol) { //the active set hasn't changed enough
					//the solution has converged for this Z; we will relinearize and start again with the new B 
					break;
				}
			}
		}

		intercept = muZ - cblas_ddot(F, B, 1, muX, 1);

		//XB = X*B + intercept
		cblas_dgemv(col_major, y, F, Nlocal, 1., X, F, B, 1, 0., XB, 1);
		cblas_daxpy(Nlocal, intercept,  ones, 1, XB, 1);

		//update P and W using the link function
		//sumWi = 0;
		for(int i=0; i<Nlocal; i++) {
			P[i] = 1./(1+exp(-XB[i]));
			if (P[i] < pe) {
				P[i] = 0;
				W[i] = pe;
				sumWi += pe;
			} else if (P[i] > mpe) {
				P[i] = 1;
				W[i] = pe;
				sumWi += pe;
			} else {
				W[i] = P[i] * (1 - P[i]);
				sumWi += W[i];
			}
		}
		//sumWi = 1. / sumWi;				

		for (int i=0; i<F; i++) {
			if (abs(Bo[i] - B[i]) > tol) {
				break;
			}
			if (i == F-1) { //none of the B have changed enough since the last time we updated Z
				//the model has converged
				return;
			}
		}
	}
}

uncertainty::uncertainty() {
	minU = INFINITY;
	r=-1;
}

void uncertainty::set_params(char* Y, double* W, int thisN, int minSize) {
	N = thisN;
	m  = minSize;
	WP = 0.; //ought to just inherit this from the parent node
	TP = 0;
	for (int i=0; i<N; i++) {
		if (Y[i] == 1) {
			WP+= W[i];
			TP++;
		}
	}
	current = 0;
}

void uncertainty::get_uncertainty(double* X, char* Y, double* W, int* inds) { //inputs: X, Y, W, N
	// std::cout << "in getu: " << N << "," << m << "," << WP << std::endl;
	for (int i=0;i<N;i++) {
		inds[i] = i; //initialize sorting...
	}
	std::stable_sort(inds, &(inds[N]), [&X](size_t a, size_t b) {return X[a] < X[b];});
	//requires cpp11 compatible compiler?
	
	// nx = 0;
	cl = 0;
	// tn = 0;

	// sw = 0.;
	tw = 0.;
	twp = 0.;
	// twn = 0.;
	for (int i=0; i<m; i++) {
		ii = inds[i];
		wi = W[ii];
		tw+=wi;
		// std::cout <<"(" <<i << "," << inds[i] << ","<< X[inds[i]] << "),";
		
		if (Y[ii] == 1) {
			cl++; //number of positive labels so far
			twp+=wi; //total weight of positive labels so far
		}
	}

	for (int i=m; i<(N-m); i++) {
		ii = inds[i];
		wi = W[ii];
		tw+=wi;

		if (Y[ii] == 1) {
			cl++; //number of positive labels so far
			twp+=wi; //total weight of positive labels so far
		}

		thisX = X[ii];
		// std::cout <<"(" <<i << "," << inds[i] << ","<< thisX << "),";
		nextX = X[inds[i+1]]; 
		if (thisX == nextX) {
			continue;
		}
		
		fw = 1.-tw;

		pl = twp / tw; //weighted probability of positive label given X <= this x
		pg = (WP - twp) / fw; //weighted prob. of pos. label given X > this x
		
		nl = 1.-pl;
		ng = 1.-pg;
		// nl = twn / tw;
		// ng = (1 - WP - twn) / fw;
		ul = - tw*((pl<tol ? 0 : pl*log(pl)) + (nl<tol ? 0 : (nl)*log(nl)));
		ur =  - fw*((pg<tol ? 0 : pg*log(pg)) + (ng<tol ? 0 : ng*log(ng)));
		u = ul + ur;
		
		if (u<minU || (abs(u-minU)<1e-4 && abs(tw-fw)<abs((1./wleft) - (1./wright)))) {
			//update the best fit if the uncertainty improves
			//or if it stays the same with more even partitioning of the weights

			minU = u;
			support = nextX - thisX;
			bias = thisX + support/2;
			r = current;
			nleft = i+1;
			nright = N-nleft;
			wleft = 1./tw;
			wright = 1./fw;
			uleft = ul;
			uright = ur;

			pleft = cl;
			pright = TP-cl;
		}

	}

	current++;
}

void paramSet::print() {
	std::cout << "Training parameters: " << std::endl;
	std::cout << "\tLogistic Regression:" << std::endl;
	std::cout << "\t\tNumber of features: " << nFeat <<std::endl;
	std::cout << "\t\tNumber of cross-validation folds: " << nFolds << std::endl;
	std::cout << "\t\tNumber of lambda values to try: " <<nLambda << std::endl;
	std::cout << "\t\tAlpha value for elastic net: " <<alpha << std::endl;

	std::cout << "\tTree parameters:" << std::endl;
	std::cout << "\t\tNumber of node repeats: " << nRepeats << std::endl;
	std::cout << "\t\tMax depth: "	<< maxDepth << std::endl;
	std::cout << "\t\tMinimum size: " <<minSize << std::endl;

	std::cout << "\tForest parameters:" << std::endl;
	std::cout << "\t\tMax tree count: " << maxTrees << std::endl;
	std::cout << "\t\tMin tree count: "	<< minTrees << std::endl;
	std::cout << "\t\tMinimum improvement: " << treeStopThresh*100 << "% over last "<< treeStopCount << " trees" << std::endl;

	std::cout << "\tEnsemble parameters:" << std::endl;
	std::cout << "\t\tNumber of forests: " << ensembleSize << std::endl;
	std::cout << "\t\tNumber of cell types: " << nLabels << std::endl;
	std::cout << "\t\tPercent of cell types in positive group: " << probNeg*100 << std::endl;
	std::cout << "\t\tPercent of cell types in negative group: " << probPos*100 << std::endl;
}

paramSet::paramSet(std::fstream& paramfile) {
	paramfile >> nFeat >> nFolds >> nLambda >> alpha; 
	paramfile >> nRepeats >> maxDepth >> minSize;
	paramfile >> maxTrees >> minTrees >> treeStopCount >> treeStopThresh;
	paramfile >> ensembleSize >> probNeg >> probPos >> nLabels;
}

Ensemble::Ensemble(paramSet* params_, short unsigned int N_, int thread_, std::filesystem::path rootDir_, int* ecoc_):params(params_),N(N_),thread(thread_),rootDir(rootDir_),ecoc(ecoc_) {
	// CreateDirectoryA(&rootDir); //create folder if it doesn't exist
	
	forests = new Forest[params->ensembleSize];//(params, this);

	inds = new int[N];
	feat = new double[params->nFeat*N];
	pred = new double[N];
	bestPred = new double[N];
	error = new bool[N];
	scoreP = new double[N];
	scoreN = new double[N];

	allPos = new double[N];
	allNeg = new double[N];


	labels = new char[N];
	sampleWeights = new double[N];
	featureWeights = new double[params->nFeat];

	// treeWeights = new double[params->nTrees];
	accuracy = new double[params->maxTrees];

	rd = new randomDraws(thread);
	features = new featureList(params->nFeat*params->nRepeats); //randomly draws features from the dataset

	reg = new regression(params, N);
}

int Ensemble::train(psthSet* data, int i) {
	// for (int i=0; i<params->ensembleSize; i++) {
	current = i;

	//assign the coding for this tree...
	data->npos = 0;	
	int stride = i*params->nLabels;
	int badCount = 0;
	int j=0;
	int err;

	while (j<N-badCount) {
		// std::cout <<(int) data->data[j]->label << ",";
		while (ecoc[stride+data->data[j]->label - 1] == 0) { //we ignore these cells

			//swap the datum at this index with another datum
			badCount++;
			if (j>=N-badCount) {
				break;
			}
			std::swap(data->data[j], data->data[N-badCount]); //swap the pointers
			// std::cout << (int) data->data[j]->label << ",";
		} 
		j++;
	}
	std::cout << std::endl;

	data->N = N-badCount; // we will ignore the cells at the tail end 
	double a = 1./(N-badCount);
	for (j=0;j<N-badCount;j++) {
		if (ecoc[stride+data->data[j]->label-1] == 1) {
			data->npos++;
			labels[j] = 1;
		} else {
			labels[j] = 0;
		}

		allPos[j] = 0.;
		allNeg[j] = 0.;
		sampleWeights[j] = a; //initialize the weights
	}

	//name/mkdir forests[i].treeDir
	forests[i].params = params;
	forests[i].ensemble = this;
	forests[i].nTrees = 0;
	forests[i].treeDir = rootDir / ("forest" + std::to_string(i+1));
	std::filesystem::create_directories(forests[i].treeDir);
	err = forests[i].train(data, i);
	if (err) {
		std::cout << "Thread[" << thread+1 << "], " << "Forest["<< i+1 << "/" << params->ensembleSize << "]  -> error!!" << std::endl;
	
		return 1;
	}

	std::cout << "Thread[" << thread+1 << "], " << "Forest["<< i+1 << "/" << params->ensembleSize << "]  -> completed" << std::endl;
	return 0;
	// }
}

Ensemble::~Ensemble() {
	delete[] inds;
	delete[] feat;
	delete[] pred;
	delete[] bestPred;
	delete[] error;
	delete[] labels;
	delete[] sampleWeights;
	delete[] featureWeights;
	// delete[] treeWeights;
	delete[] accuracy;
	delete[] scoreP;
	delete[] scoreN;
	delete[] allPos;
	delete[] allNeg;
	delete rd;
	delete features;
	delete reg;

	delete[] forests;
}

Forest::Forest() {
}

int Forest::train(psthSet* data, int f) {
	// trees = new Tree*[params->maxTrees];//(params, ensemble);
	// trees = Tree*;

	sampleWeights = ensemble->sampleWeights;

	double sumWeights;
	int comp;
	double compAcc;
	bool improved;

	double* accuracy = ensemble->accuracy;
	double* allPos = ensemble->allPos;
	double* allNeg = ensemble->allNeg;


	// double allTreeWeights = 0.;
	std::ofstream treefile;

	double* scoreP = ensemble->scoreP;
	double* scoreN = ensemble->scoreN;
	char* labels = ensemble->labels;

	double a = 1./data->N;

	int err = 0;

	//do training
	while (true) {
		//check stopping criteria -- want it at top in case we interrupted training...
		if (nTrees == params->maxTrees) {
			break;
		} else if (nTrees >= params->minTrees) { //we're allowed to stop if the tree stopped improving
			improved = false;
			comp = nTrees-params->treeStopCount;
			compAcc = accuracy[comp] + params->treeStopThresh;
			for (int j=nTrees-1; j>comp; j--) {
				if (accuracy[j] >= compAcc) {
					improved = true;
					break;
				}
			}
			if (!improved) {
				break;
			}
		}
		trees = new Tree;

		//train the next tree
		trees->params = params;
		trees->ensemble = ensemble;
		trees->init();
		err = trees->train(data);
		if (err) {
			return 1;
		}

		//get the cells that were incorrectly classsified
		eta = 0.;
		for(int j=0; j<data->N; j++) {
			if (ensemble->error[j]) {
				eta += sampleWeights[j];	
			}
		}

		// ensemble->treeWeights[nTrees] = log((1-eta+1e-9)/(eta+1e-9)); 
		trees->treeWeight = log((1-eta+1e-9)/(eta+1e-9));//add an epsilon to avoid divide by zero

		//save tree with treeweight (before updating sampleweights)
		treefile.open(treeDir / ("tree" + std::to_string(nTrees+1) + ".out"), std::ios::binary);
		trees->save(treefile);
		treefile.close();

		accuracy[nTrees] = 0.;
		for (int i=0; i<data->N; i++) {
			allPos[i] += scoreP[i]*trees->treeWeight; //if tree weight is negative it subtracts from both pos and neg...
			allNeg[i] += scoreN[i]*trees->treeWeight;

			// std::cout <<"("<< allPos[i] << "," << allNeg[i] << "," << (int) labels[i] << "),"; 
			accuracy[nTrees] += ((labels[i]==1 && allPos[i]>allNeg[i]) || (labels[i]==0 && allPos[i]<allNeg[i]))? 1 : 0;
		}
		// std::cout << std::endl;
		accuracy[nTrees] *=a;
		
		std::cout << "Thread[" << ensemble->thread+1 << "], Forest[" << f+1 <<"/" << params->ensembleSize << "], Tree[" << nTrees+1 <<"/"<< params->maxTrees << "] -> " << accuracy[nTrees]*100 <<"% accurate " << std::endl;
		// std::cout << "this tree acc: " << trees->accuracy << ", overall: " << accuracy[nTrees] << std::endl;
		
		sumWeights = 0.;
		for(int j=0; j<data->N; j++) {
			if (ensemble->error[j]) {
				sampleWeights[j] *=exp(trees->treeWeight);
			}
			sumWeights+= sampleWeights[j];
		}
		
		sumWeights = 1./sumWeights;
		for (int j=0; j<data->N; j++) {
			sampleWeights[j] *= sumWeights; //re-normalize weights
		}
		delete trees;
		nTrees++;
	}
	return 0;
}

Forest::~Forest() {
	// delete[] trees;
}