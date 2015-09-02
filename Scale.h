// scaling initiator matrix
void ScaleFitMtx(int ModelNodes, int ModelEdges, TFlt MaxDegInModel, TFlt MaxDegOutModel, TKronMtx& FitMtx, const TInt& NIter, TStr IsDir, const TStr& ScaleMtx, ofstream& TFile);
void ScaleFitMtxForUnDir(TKronMtx& FitMtx, ofstream& TFile);
void ScaleFitMtx(TKronMtx& FitMtx, const TInt& NIter, const int& InitModelNodes, const int& MaxModelDeg, const TStr& IsDir);
void ScaleFitMtxForEdges(TKronMtx& FitMtx, const TInt& NIter, const int& ExpectedModelEdges, ofstream& TFile);
// estimate scaling coefficient
double GetScalingCoefficient(const TFltPrV& InDegCnt, const TFltPrV& OutDegCnt, const TKronMtx& FitMtx, const TInt& NIter, const TStr& IsDir);