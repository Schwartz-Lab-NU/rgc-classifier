#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <cblas.h>
#include <filesystem>
#include <sstream>

class randomDraws {
public:
	std::mt19937 mt; //creates mersenne twister
	std::uniform_real_distribution<double> dist{std::uniform_real_distribution<double>(0.0, 1.0)};

	randomDraws(int thread);
	randomDraws(const randomDraws& rd);
	randomDraws(){};
	std::mt19937 make_twister(int thread);

	double getR();
	void getP(int* ind, int N);
};

class psth { // stores the psth data and metadata for a particular cell
	//all of the data should be binned the same, meaning we want to impute the data at the beginning
	//as we move down the tree we only want to copy the pointer to each cell's data, since it is read-only

public: 
	std::string cellName; //helps us keep track of which cells used to train
	char label; // stores labels as signed integer

	double* data;

	unsigned int* qualityS;
	char* qualityT;

	~psth();
	psth(){};
	// psth(const psth& indata, int nt, int ns);
};

class psthSet {//stores a collection of cells
public:
	unsigned short int N;
	unsigned short int npos; //number of positive labels
	psth** data; // a vector of pointers to the individual cells 

	double dt; //number of time points
	double t0;
	double tf;
	unsigned short int nt;

	unsigned short int ns; //number of spot sizes
	unsigned short int s0;
	unsigned short int sf;

	unsigned long int ne;

	double* pdfS; //reflects the counts of spots of different sizes 
	double* pdfT;

	double* cdfS; //integral of pdf with nearest-neighbors blurring
	double* cdfT;

	unsigned int* qualityT;
	unsigned long long int* qualityS;

	// double* weights;

	bool light; // true if we don't need to alloc/delete the individual psth's
	bool empty; //if true then data is a shallow copy of another psthset's data

	void stats();
	void append(psthSet* indata);
	void copy(psthSet* data);
	psthSet(int N);
	psthSet(const psthSet* indata);
	psthSet(std::ifstream& datafile);
	~psthSet();
};

class featureItem {
public:
	int t; //the time points selected
	int s; //the spot sizes selected

	double mean = 0;
	double istd = 0;

	double qualityT;
	double qualityS;
	double quality;
};

class featureList {
private:
	double j;
	int indT;
	int indS;
	int ind;
	double delta;

public:
	featureItem* feats;
	unsigned int F;

	featureList(int F_);
	~featureList();

	void get(int N, psthSet* data, randomDraws* r, double* sampleWeights);
};

class paramSet {//stores the parameters for the tree
public:
	//int nCodes = 0; //number of coding schemes to try at each node
	
	int nFeat = 20; //number of features to include in logistic regression
	int nFolds = 5; //number of cross-validation folds for regression
	int nLambda = 25; //number of lambda values to try
	double alpha = 0.5; //the balance between lasso and ridge regression
	// double nullDevFrac = 0.05; //if not cross-validating, use the largest lambda that is less than this fraction of the null

	int nRepeats = 3; //number of repetitions at each node
	int maxDepth = 11; //how deep to grow the tree
	int minSize = 3; //the smallest possible node size

	int maxTrees = 100; //always stop when this many trees are trained
	int minTrees = 80; //can't stop before this many trees are trained
	int treeStopCount = 20; //look at the delta over the last x trees
	double treeStopThresh = .02; //if the accuracy doesn't improve by this much we assume we're done

	int ensembleSize = 100;	
	double probNeg = 0.5; //probability of assigning label to negative class
	double probPos = 0.5;

	int nLabels = 31;

	paramSet() {};
	paramSet(std::fstream& paramfile);
	void print();
	void print(std::fstream&);
};

class regression {
private:
	//scalars
	double intercept;
	double muY;
	double lmuY;
	double lnmuY;
	double muZ;
	double alpha;
	double denoL;
	double gamma;
	double thisLambda;
	double lambdaMax;
	double lambdaMin;

	double sumWi;
	int stride, stride2;
	int ind;
	bool activeFlag;
	double maxD;
	double delta;
	double bp;
	double margin;

	int maxIterI = 1e4;
	int maxIterO = 1e2;
	double pe = 1e-5;
	double mpe = .99999;
	double tol = 1e-4;

	//1-by-F
	double* B;
	double* Bo;
	double* Bi;
	double* muX;
	double* lmuX;
	double* deno;
	bool* active;
	bool* constant;

	double* wX2;

	//1-by-N
	double* P;
	double* W;
	double* ones; //double for compatibility with BLAS?
	double* negones;
	double* Z;
	double* XB;
	double* stdX;

	//F-by-N
	double* Xc;
	double* wX;
	double* X2;

	//other
	int localN;
	int localNp, localNn;
	double* lastDev;
	double* seDev;
	int* inds;
	double devThresh;
	double sqrtN;
	double* localX;
	char* localY;
	double* localSW;

	#ifdef CBLAS_LAYOUT
	CBLAS_LAYOUT col_major;
	#else
	CBLAS_ORDER col_major;
	#endif

	CBLAS_TRANSPOSE y, n;
public:
	int F;
	double lambda;
	int nLambda;
	double deviance;
	int nFolds;

	double beta0;
	double* beta1;

	regression(paramSet* params, int maxN);
	~regression();

	void do_fit(double* X, char* Y, double* sw, double* fw, int N, int npos,randomDraws* r);
	void irls(double* X, char* Y, double* Ws, int Nlocal);
};

class uncertainty {
private:
	int m;
	int* inds;
	int current;

	double WP;
	double wi;
	double tw;
	double twp;
	
	int ii;
	int TP;
	double pl;
	double pg;
	double nl;
	double ng;

	double fw;
	double cl;

	double u;
	double ul;
	double ur;

	double thisX;
	double nextX;

	double tol = 1e-4;
public:
	unsigned int N;
	double bias; //the primary output
	
	double minU; // bias <- argmin ( U ), reflects the purity of the daughter nodes
	double support; // the distance between the samples that immediately neighbor the threshold
	//support reflects the stability of the threshold

	unsigned int r; //the index of the chosen repetition
	unsigned int nleft; //number of samples in resultant left daughter node
	unsigned int nright;

	unsigned int pleft;
	unsigned int pright;

	double wleft; //inverse of the total weight in the left daughter node, for renormalization
	double wright;
	double uleft;
	double uright;

	uncertainty();
	void set_params(char* Y, double* W, int thisN, int minSize);
	void get_uncertainty(double* X, char* Y, double* W, int* inds);
};

class Tree;
class Ensemble;
class Node {
public:
	unsigned int N; //number of cells at this node
	int misclassified; //number of misclassified cells at this node
	double accuracy; //accuracy at this node
	char MLE; //-1 for neg, 0 for equal, 1 for pos
	double sumWeights;

	Node* parent; //pointer to parent node
	Tree* tree;
	bool leftRight; //0 if node is left of parent, 1 if right
	int depth;


	Node() {};
	// virtual void train(psthSet* data, paramSet* params, randomDraws* rd, regression* reg);
	virtual int train(psthSet* data, double sumWeights_);
	virtual void test(psth** data, double* results, double* vals, int* inds, int end, int nt){};
	virtual void save(std::fstream& datafile);
	virtual void load(std::fstream& datafile);
	virtual void print(int d, int o, bool v);

	virtual ~Node(){};
};

class Leaf : public Node{
public:
	unsigned int pos;
	unsigned int neg;
	double score;
	// double fracPos;
	// double fracNeg;
	int start; //for training...
	// double sumWeights;

	Leaf(){};
	Leaf(int start_);
	int train(psthSet* data, double sumWeights_);
	void test(psth** data, double* results, double* vals, int* inds, int end, int nt);
	void save(std::fstream& datafile);
	void load(std::fstream& datafile);
	void print(int d, int o, bool v);
	~Leaf(){};
};
// 
class compactDecisionNode : public Node{
public:
	double* t; //shoud be int!!!
	double* s;

	double* beta1;
	double* mean;
	double* istd;
	double beta0;

	unsigned int nnz;

	Node* left; //pointer to left child
	Node* right; //pointer to right child

	bool complete = false;
};

class decisionNode : public compactDecisionNode {
public:
	double* qualityT;
	double* qualityS;
	double* quality;

	double lambda; // value of chosen lambda, the regression penalizer
	double deviance; // deviance of regression fit from binomial model
	double classU; // reflects the entropy gain of the chosed split wrt the labels
	double support; // the distance from the threshold to the nearest measurements

	int start;

	decisionNode(int N_, Tree* tree_, int start_);
	decisionNode(){};
	int train(psthSet* data, double sumWeights_);
	void test(psth** data, double* results, double* vals, int* inds, int end, int nt);
	void save(std::fstream& datafile);
	void load(std::fstream& datafile);
	void print(int d, int o, bool v);
	~decisionNode();
};

class Tree {
public:
	Ensemble* ensemble;
	int N;
	int* inds;
	double* feat;
	double* pred;
	double* bestPred;
	bool* error;
	double* score;
	// double* scoreP;
	// double* scoreN;
	// double* allPos;
	// double* allNeg;
	double* allScores;

	char* labels;
	double* sampleWeights;
	double* featureWeights;
	bool light = true;

	double treeWeight;
	double accuracy;

	bool alloced = false;

	featureList* features;

	Node* root;
	paramSet* params;
	regression* reg;
	randomDraws* rd;
	randomDraws* rd_original;

	int minSize;

	Tree(paramSet* params_);
	Tree(std::fstream& datafile);
	Tree();
	~Tree();

	void init();
	int train(psthSet* data);
	void test(psthSet* data, double* results, double* vals, int* inds);
	void save(std::fstream& datafile);
	void print(bool o, bool v);
};

class Forest {
private:
	double eta; //the total weight in error

	// featureList* features;
	// regression* reg;
	// randomDraws* rd;

public:
	int nTrees;

	Ensemble* ensemble;
	Tree* trees;
	paramSet* params;

	double* sampleWeights;
	// double* treeWeights;

	// double* accuracy; // for tracking improvements over time...

	std::filesystem::path treeDir; //for saving 

	Forest();
	Forest(paramSet* params_, Ensemble* ens_);
	int train(psthSet* data, int f, double* vals); //and some more stuff...
	void test(psthSet* data, double* results, double* vals, int* inds);
	void test(psthSet* data, double* results, double* vals, int* inds, int* saveN);
	~Forest();
};

class Ensemble {
public:
	// char* coding;
	// int nLabels;

	paramSet* params;
	int N;
	int thread;
	std::filesystem::path rootDir;
	int* ecoc;

	bool light = true;

	int current;

	Forest* forests;

	int* inds;
	double* feat;
	double* pred;
	double* bestPred;
	bool* error;
	double* score;
	// double* scoreP;
	// double* scoreN;

	// double* allPos;
	// double* allNeg;
	double* allScores;

	char* labels; //really can be bool at this point, but same memory space
	double* sampleWeights;
	// double* treeWeights;
	double* accuracy;

	double* featureWeights;

	randomDraws* rd;
	featureList* features;
	regression* reg;

	std::stringstream buffer;


	//Ensemble(paramSet* params_, int nLabels_);
	Ensemble(paramSet* params_, short unsigned int N_, int thread_, std::filesystem::path rootDir_, int* ecoc_);
	Ensemble(){};
	int train(psthSet* data, int* sampleCounts, int i);
	double test(psthSet* data, char* pred, int* sampleCounts, char transform);
	void normalize(double* R, double* X, double* Y, size_t N_);
	void normalize(psthSet* data, int* sampleCounts, int f);
	double normalize(double*, int*, int*, int*, char*);
	void pav(char* Y_, double* Y, double* W, size_t N_);
	void get_r_hat(double*, double*, double*, int*, int*);

	~Ensemble();

};