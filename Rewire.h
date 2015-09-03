typedef pair<TFltPr, TFltV> Diap;
typedef map<int, vector<pair<int, int>>> RewireDiap;
typedef map<int, map<pair<int,int>, vector<int>>> ClusterMap;

// get relative differences of degrees
void GetRelativeDiff(const TFltPrV& MDeg, const TFltPrV& KronDeg, TFltPrV&  RelDiff);
// get base diapasons for scaling
void GetBaseDiaps(const TFltPrV& MDeg, const TFltPrV& KronDeg, vector<BaseDiap>& BaseDiaps);
// rewire edges according to smoothed diaps
void Rewire(PNGraph& Kron, const vector<BaseDiap>& BaseDiaps, const TIntPr& OutDegR, ofstream& TFile);
void Rewire(PNGraph& Kron, vector<Diaps>& DPlus, vector<Diaps>& DMinus, const int DegMin, const int DegMax);
// add missing or delete excess edges
void AddEdges(PNGraph&Kron, int Diff, int DegMin, int DegMax, int ModelEdges, vector<Diaps>& DMinus, vector<Diaps>& DPlus);
// get appropriate number of nodes for each diapason and its average degree
void GetDiaps(vector<Diaps>& DPlus, vector<Diaps>& DMinus, const vector<BaseDiap>& BaseDiaps, const TFltPrV& KronDeg, const TInt& DegMin, const TInt& DegMax);
// get random number from diap
int GetRndDeg(TRnd& Rnd, const TIntPr& Borders);
// get random number from diap with account of possibilities
int GetRndDeg(TRnd&Rnd, const TIntPr& Borders, const TFltV& DiapPossib);
// get rewire strategies
int GetRewireStrategies(vector<Diaps>& DPlus, vector<Diaps>& DMinus);
// get diap index
bool GetDiap(int Deg, vector<Diaps>& DPlus, vector<Diaps>& DMinus, int& DiapIndex, bool& IsDPlus);
// add random edge
bool AddRndEdge(TRnd& Rnd, PNGraph&Kron, int Node, int DegMax);
// get average deviation
double GetAvgDeviation(const TFltPrV& ModelDegCnt, const TFltPrV& KronDegCnt);