#include "stdafx.h"
#include "Error.h"
#include <algorithm>

// calculate binomial coefficient
int GetBinomCoeff(const int& N, const int& K){
	if (N < 0 || K < 0)
		Error("GetBinomCoeff", "N or K is negative");
	if (N < K)
		Error("GetBinomCoeff", "N < K");
	int Term = 1;
	for (int i = 0; i < K; i++){
		Term *= (N - i);
	}
	int Denom = 1;
	for (int i = 2; i <= K; i++)
		Denom *= i;
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

