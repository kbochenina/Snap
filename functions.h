


// number of additional functions
const int NFUNC = 9;
const int NPARCOPY = 4;
const string FUNCNAMES[] = {"KRONTEST", "GRAPHGEN_M", "MTXGEN_M", "KRONFIT_M", "KRONGEN_M", "GRAPHGEN_MS", "MTXGEN_MS", "KRONFIT_MS", "KRONGEN_MS"};
const enum CMDARGS { KRONTEST = 0, GRAPHGEN_M = 1, MTXGEN_M = 2, KRONFIT_M = 3, KRONGEN_M = 4, GRAPHGEN_MS = 5, MTXGEN_MS = 6, KRONFIT_MS = 7, KRONGEN_MS = 8};
const enum ARGS { GRAPHGEN = 0, MTXGEN = 1, KRONFIT = 2, KRONGEN = 3, PTYPE = 4, PLT = 5, NAME = 6, NEIGEN = 7, BINRADIX = 8, HOPS = 9, CLUST = 10};


// generates Kronecker model on the base of existing graph
void KroneckerTest(vector<TStr> commandLineArgs);
// generates Kronecker model using configuration model of small model network
// and compare it to big network
void KroneckerByConf(vector<TStr> commandLineArgs);
// generates set of graphs + calculates the metrics + create plots
void GetGraphs(const vector <TStr>& parameters, const TStr& ModelGen, const TStr&ModelPlt);
// get FitMtx and scaling coefficient from small model
void GetFitMtxFromMS(TKronMtx& FitMtx, TFlt& ScalingCoeff, vector<BaseDiap>& BaseDiaps, const vector<TStr>& Parameters, vector<int>& Prev);
// graph generator by Lescovec
int BasicGraphGen(const TStr args, PNGraph &GD);
// create string with parameters of model graph as an input to GenKron() function
TStr GetModelParamsStr(const PNGraph& G);
// generates Kronecker graphs
void GenKron(const TStr& Args, TKronMtx& FitMtx, TFltPrV& KronDegAvgIn, TFltPrV& KronDegAvgOut, vector<BaseDiap>& BaseDiaps, PNGraph& out);
// estimation of initiator matrix
int InitKronecker(const TStr args, PNGraph &G, TKronMtx& FitMtx);
// generates one instance of Kronecker graphs 
void KroneckerGen(PNGraph& out, const TKronMtx& FitMtx, const TInt NIter, const TStr& IsDir, const TIntPr& InDegR, const TIntPr& OutDegR, double NoisePart);
// get average estimates of out-degree of NKron Kronecker graphs
void GetAvgKronDeg(const TKronMtx& NewMtx, const TInt& NIter, const TStr& IsDir, const TInt& NKron, const TIntPr& ModelDegR, TFltPrV& KronDeg);


