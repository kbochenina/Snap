#include <fstream>
#include <string>

// number of additional functions
const int NFUNC = 9;
const int NPARCOPY = 4;
const string FUNCNAMES[] = {"KRONTEST", "GRAPHGEN_M", "MTXGEN_M", "KRONFIT_M", "KRONGEN_M", "GRAPHGEN_MS", "MTXGEN_MS", "KRONFIT_MS", "KRONGEN_MS"};
const enum CMDARGS { KRONTEST = 0, GRAPHGEN_M = 1, MTXGEN_M = 2, KRONFIT_M = 3, KRONGEN_M = 4, GRAPHGEN_MS = 5, MTXGEN_MS = 6, KRONFIT_MS = 7, KRONGEN_MS = 8};
const enum ARGS { GRAPHGEN = 0, MTXGEN = 1, KRONFIT = 2, KRONGEN = 3, PTYPE = 4, PLT = 5, NAME = 6, NEIGEN = 7, BINRADIX = 8, HOPS = 9, CLUST = 10};

void DefaultSettings(vector<TStr>& out);
void ReadParameters(TStr settingsFN, vector<TStr>& out);
int GraphGen(const TStr args, PNGraph &GD);
int InitKronecker(const TStr args, PNGraph &G, TKronMtx& FitMtx);
double KroneckerGen(const TStr args, const TKronMtx& FitMtx, PNGraph& out, double ModelClustCf = 0.0);
void ReadPNGraphFromFile(const TStr args, PNGraph& out);
int FindVal1Elem(const TIntPrV& vec, const TInt& elem, bool& isExact);
int FindVal1Elem(const TFltPrV& vec, const TFlt& elem, bool& isExact);
// generates Kronecker model on the base of existing graph
void KroneckerTest(vector<TStr> commandLineArgs);
// generates Kronecker model using configuration model of small model network
// and compare it to big network
void KroneckerByConf(vector<TStr> commandLineArgs);
// generates set of graphs + calculates the metrics + create plots
void GetGraphs(const vector <TStr>& parameters, const TStr& ModelGen, const TStr&ModelPlt);
// create string with parameters of model graph as an input to GenKron() function
TStr GetModelParamsStr(const PNGraph& G);
// generates Kronecker graphs
void GenKron(const TStr& Args, TKronMtx& FitMtx, TFltPrV& KronDegAvgIn, TFltPrV& KronDegAvgOut);
// generates one instance of Kronecker graphs 
void KroneckerGen(PNGraph& out, const TKronMtx& FitMtx, const TInt NIter, const TStr& IsDir, const TIntPr& InDegR, const TIntPr& OutDegR, double NoisePart);

void RemoveUnusedNodes(PNGraph& out, const int& MinDeg);
// scaling initiator matrix
void ScaleFitMtxForUnDir(TKronMtx& FitMtx);
void ScaleFitMtx(TKronMtx& FitMtx, const TInt& NIter, const int& InitModelNodes, const int& MaxModelDeg, const TStr& IsDir);