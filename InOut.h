void DefaultSettings(vector<TStr>& out);
void ReadParameters(TStr settingsFN, vector<TStr>& out);
void CheckParams(const TStr& model_gen, const TStr& model_plt);
void GetMtxFromSepLine(const TStr& line, const TStr& separator, TFltV& matrix);
void ReadPNGraphFromFile(const TStr args, PNGraph& G);
ofstream OpenFile(const TStr& fileName);
void PrintDegDistr(const TFltPrV& distr, const TStr& OutFNm);
void PrintDegDistr(const TIntPrV& distr, const TStr& OutFNm);
void PrintDegDistr(const PNGraph& G, const TStr& OutFNm);
void ReadMtx(const TStr& Mtx, const TInt& MtxSize, TKronMtx& FitMtx);
double PrintLargestEigenVal(const PNGraph& G, ofstream& F, const TStr& GName);
void MakeDatFile(const TStr& Name, const TStr& AddStr, const TStrV& ColumnNames, const vector<vector<double>>& Data, const int& Nodes = 0, const int& Edges = 0);
void PrintNodeDegrees(const PNGraph& G, const TKronMtx& FitMtx, const int& NIter);
void PrintDegDistr(const TFltPrV& distr, const TStr& OutFNm);
void PrintDegDistr(const TIntPrV& distr, const TStr& OutFNm);
void GetParameters(const vector<TStr>& CommandLineArgs, const TStr& Type, vector<TStr>& Parameters);
void PrintRelDiff(const TFltPrV& RelDiffV, const TStr& OutFNm);
void PrintSmoothedDiaps(const vector<pair<TFltPr, TFlt>>& Diaps, const TStr& OutFNm);

