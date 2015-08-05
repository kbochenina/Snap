
void AddDegreeStat(const TFltPrV& current, TFltPrV& degrees, TIntPrV& samples);
void AddDegreesStat(TFltPrV& deg, TIntPrV& samples, const PNGraph& G, bool isIn);
void GetAvgDegreeStat (TFltPrV& deg, const TIntPrV& samples);
void GetAvgDegreeStat (TFltPrV& deg, const TInt& NKron);
void GetPoints(const TFlt& maxDegLog, const TFlt& minDegLog, const int& NInt, const TFltPrV& base, TFltPrV& points);
int GetMaxMinDeg(const PNGraph& G, const TStr& IsDir, const TStr& IsIn, const TStr& IsMax);
void CompareDeg(const int i, const int MaxDeg, int& MinMaxDeg, int& MaxMaxDeg, int& AvgMaxDeg);
int GetExpectedModelEdges(const PNGraph& G, const int k, const TStr& order);
bool CheckReciprocity(const PNGraph& G);