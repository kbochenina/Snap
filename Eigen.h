void GetEigVProbMtx(const double& EigMax, const double& EigMin, const int& NIter, vector<double>& EigV, vector<double>& Mult);
int GetBinomCoeff(const int& N, const int& K);
void PrintEigen(const TKronMtx& FitMtx, const int &NIter, const int& NEigen);
// ModelEigValV - output parameter
void PlotEigen(const PNGraph& G, const TStr& NEigenStr, const TStr& Name, TFltV& ModelEigValV);
