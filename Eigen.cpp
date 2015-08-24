#include "stdafx.h"
#include "Error.h"
#include <algorithm>

// calculate binomial coefficient
int GetBinomCoeff(const int& N, const int& K){
	if (N < 0 || K < 0)
		Error("GetBinomCoeff", "N or K is negative");
	if (N < K)
		Error("GetBinomCoeff", "N < K");
	long long int Term = 1;
	for (int i = 0; i < K; i++){
		Term *= (N - i);
	}
	long long int Denom = 1;
	for (int i = 2; i <= K; i++)
		Denom *= i;
	int Res = Term / Denom;
	if (Res < 0)
		Error("GetBinomCoeff", "Res is negative");
	return Term / Denom;
}

// return eigenvalues of 2x2 probability matrix for NIter Kronecker multiplications
// eigenvalues are placed in EigV vector in descending order
void GetEigVProbMtx(const double& EigMax, const double& EigMin, const int& NIter, vector<double>& EigV, vector<double>& Mult){
	for (int i = 0; i <= NIter; i++){
		int A1 = NIter - i, A2 = i;
		int K = A1 > A2 ? A1 : A2;
		double EigVal = pow(EigMax, A1) * pow(EigMin, A2);
		EigV.push_back(EigVal);
		Mult.push_back(GetBinomCoeff(NIter, K));
	}
	for (int i = 0; i < EigV.size() - 1; i++){
		for (int j = i+1; j < EigV.size(); j++){
			if (EigV[i] < EigV[j]){
				double Tmp = EigV[i];
				EigV[i] = EigV[j];
				EigV[j] = Tmp;
				Tmp = Mult[i];
				Mult[i] = Mult[j];
				Mult[j] = Tmp;
			}
		}
	}
	//sort(EigV.begin(), EigV.end(), std::greater<double>());
}

void PrintEigen(const TKronMtx& FitMtx, const int &NIter, const int& NEigen){
	vector<double> EigProbMtx, Mult;
	TFlt EigMax = FitMtx.GetEigMax(), EigMin = FitMtx.GetEigMin();
	printf("EigMax = %f, EigMin = %f\n", EigMax, EigMin);
	GetEigVProbMtx(EigMax, EigMin, NIter, EigProbMtx, Mult);
	if (EigProbMtx.size() == 0)
		Error("GenKron", "EigProbMtx has zero size");
	TStrV ColumnNames; ColumnNames.Add("Rank");ColumnNames.Add("EigVal");
	vector<vector<double>> Data;

	vector<double> EigMult;
	
	int ValuesAdded = 0;

	for (int i = 0; i < EigProbMtx.size(); i++){
		for (int j = 0; j < Mult[i]; j++){
			ValuesAdded++;
			if (ValuesAdded == NEigen / 2)
				break;
			EigMult.push_back(EigProbMtx[i]);
		}
		if (ValuesAdded == NEigen / 2)
			break;
	}

	vector<double> Rank;
	for (int i = 0; i < EigMult.size(); i++)
		Rank.push_back(i + 1);
	Data.push_back(Rank);
	Data.push_back(EigMult);

	TFlt EigMaxProb = pow(EigMax, NIter);
	TStr AddStr("Largest eigen val = "), AddVal(EigMaxProb.GetStr());
	AddStr.InsStr(AddStr.Len(), AddVal);

	MakeDatFile("eigVal.ProbMtx", AddStr, ColumnNames, Data);
}

void PlotEigen(const PNGraph&G, const TStr& NEigenStr, const TStr& Name){
	TFltV EigValV;
	TSnap::PlotEigValRank(TSnap::ConvertGraph<PUNGraph>(G), NEigenStr.GetInt(), Name + "Eigen");
}