#include <fstream>
#include <string>
#include <algorithm>
#include <map>

// number of additional functions
const int NFUNC = 9;
const int NPARCOPY = 4;
const string FUNCNAMES[] = {"KRONTEST", "GRAPHGEN_M", "MTXGEN_M", "KRONFIT_M", "KRONGEN_M", "GRAPHGEN_MS", "MTXGEN_MS", "KRONFIT_MS", "KRONGEN_MS"};
const enum CMDARGS { KRONTEST = 0, GRAPHGEN_M = 1, MTXGEN_M = 2, KRONFIT_M = 3, KRONGEN_M = 4, GRAPHGEN_MS = 5, MTXGEN_MS = 6, KRONFIT_MS = 7, KRONGEN_MS = 8};
const enum ARGS { GRAPHGEN = 0, MTXGEN = 1, KRONFIT = 2, KRONGEN = 3, PTYPE = 4, PLT = 5, NAME = 6, NEIGEN = 7, BINRADIX = 8, HOPS = 9, CLUST = 10};

typedef pair<TFltPr, TFlt> Diap;
typedef map<int, vector<pair<int, int>>> RewireDiap;
typedef map<int, map<pair<int,int>, vector<int>>> ClusterMap;

// generates Kronecker model on the base of existing graph
void KroneckerTest(vector<TStr> commandLineArgs);
// generates Kronecker model using configuration model of small model network
// and compare it to big network
void KroneckerByConf(vector<TStr> commandLineArgs);
// generates set of graphs + calculates the metrics + create plots
void GetGraphs(const vector <TStr>& parameters, const TStr& ModelGen, const TStr&ModelPlt);
// get FitMtx and scaling coefficient from small model
void GetFitMtxFromMS(TKronMtx& FitMtx, TFlt& ScalingCoeff, vector<Diap>& SmoothedDiaps, const vector<TStr>& Parameters, vector<int>& Prev);
// graph generator by Lescovec
int BasicGraphGen(const TStr args, PNGraph &GD);
// create string with parameters of model graph as an input to GenKron() function
TStr GetModelParamsStr(const PNGraph& G);
// generates Kronecker graphs
void GenKron(const TStr& Args, TKronMtx& FitMtx, TFltPrV& KronDegAvgIn, TFltPrV& KronDegAvgOut, vector<Diap>& SmoothedDiaps, vector<int>& Prev);
// estimation of initiator matrix
int InitKronecker(const TStr args, PNGraph &G, TKronMtx& FitMtx);
// generates one instance of Kronecker graphs 
void KroneckerGen(PNGraph& out, const TKronMtx& FitMtx, const TInt NIter, const TStr& IsDir, const TIntPr& InDegR, const TIntPr& OutDegR, double NoisePart);
// scaling initiator matrix
void ScaleFitMtx(int ModelNodes, int ModelEdges, TFlt MaxDegInModel, TFlt MaxDegOutModel, TKronMtx& FitMtx, const TInt& NIter, TStr IsDir, const TStr& ScaleMtx);
void ScaleFitMtxForUnDir(TKronMtx& FitMtx);
void ScaleFitMtx(TKronMtx& FitMtx, const TInt& NIter, const int& InitModelNodes, const int& MaxModelDeg, const TStr& IsDir);
void ScaleFitMtxForEdges(TKronMtx& FitMtx, const TInt& NIter, const int& ExpectedModelEdges);
// estimate scaling coefficient
double GetScalingCoefficient(const TFltPrV& InDegCnt, const TFltPrV& OutDegCnt, const TKronMtx& FitMtx, const TInt& NIter, const TStr& IsDir);
double GetAvgDeviation(const TFltPrV& ModelDegCnt, const TFltPrV& KronDegCnt);
// get average estimates of out-degree of NKron Kronecker graphs
void GetAvgKronDeg(const TKronMtx& NewMtx, const TInt& NIter, const TStr& IsDir, const TInt& NKron, const TIntPr& ModelDegR, TFltPrV& KronDeg);
// get relative differences of degrees
void GetRelativeDiff(const TFltPrV& MDeg, const TFltPrV& KronDeg, TFltPrV&  RelDiffNonCum, bool NonCum = true);
// get smoothed diapasons for scaling
void GetSmoothedDiaps(const TFltPrV& RelDiffNonCum, vector<Diap>& SmoothedDiaps, vector<int>& Prev);
// rewire edges according to smoothed diaps
void Rewire(PNGraph& Kron, const vector<Diap>& SmoothedDiaps, const TIntPr& OutDegR, vector<int>& Prev);
void Rewire(PNGraph& Kron, RewireDiap& DiapsToCluster, RewireDiap& DiapsToDel, const TIntPrV& DiapBorders, const int DegMin, const int DegMax);
// add missing or delete excess edges
void AddEdges(PNGraph&Kron, int Diff, int DegMin, int DegMax, int ModelEdges);
// get appropriate number of nodes for each diapason and its average degree
void GetDiapNodes(TIntPrV& DiapNodes, TIntPrV& DiapBorders, const vector<Diap>& SmoothedDiaps, const TFltPrV& KronDeg, const TInt& DegMin, const TInt& DegMax, const vector<int>& Prev);
// get random number from diap
int GetRndDeg(TRnd& Rnd, const TIntPr& Borders);
// get rewire strategies
int GetRewireStrategies(RewireDiap& DiapsToCluster, RewireDiap& DiapsToDel, TIntPrV& DiapNodes, const TIntPrV& DiapBorders);
// get diap index
bool GetDiap(TInt& Deg, const TIntPrV& DiapBorders, TInt& DegIndex);
// add random edge
bool AddRndEdge(TRnd& Rnd, PNGraph&Kron, int Node, int DegMax);

