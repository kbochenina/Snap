void DefaultSettings(vector<TStr>& out);
void ReadParameters(TStr settingsFN, vector<TStr>& out);
void CheckParams(const TStr& model_gen, const TStr& model_plt);
void GetMtxFromSepLine(const TStr& line, const TStr& separator, TFltV& matrix);
void ReadPNGraphFromFile(const TStr args, PNGraph& G);
ofstream OpenFile(const TStr& fileName);
void PrintDegDistr(const TFltPrV& distr, const TStr& OutFNm);
void PrintDegDistr(const TIntPrV& distr, const TStr& OutFNm);
void ReadMtx(const TStr& Mtx, const TInt& MtxSize, TKronMtx& FitMtx);