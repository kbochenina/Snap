#include "StdAfx.h"
#include "Diaps.h"
#include "Error.h"

Diaps::Diaps(int I, pair<int, int> B, int BL) : L(Borders.first), R(Borders.second){
	if (I < 0)
		Error("Diaps::Diaps", "I < 0");
	if (B.first < 0 || B.second < 0)
		Error("Diaps::Diaps", "Negative border");
	if (B.first > B.second)
		Error("Diaps::Diaps", "B.first > B.second");
	if (BL > B.second - B.first + 1)
		Error("Diaps::Diaps", "Base length is more that real length");
	Index = I;
	Borders.first = B.first;
	Borders.second = B.second;
	L = Borders.first;
	R = Borders.second;
	BaseLen = BL;
	Len = R - L + 1;
	SetSubB();
	StratNodes = 0;
}

void Diaps::SetNodes(int N){
	if (N == 0)
		Error("Diaps::SetNodes", "N == 0");
	Nodes = N;
}

// set subborders
void Diaps::SetSubB(){
	double NodesPerDiap = Len / BaseLen;
	double AccNodes = NodesPerDiap;
	int SubBBegin = L, SubBEnd = L;
	while (SubBEnd <= R){
		while (SubBEnd - L + 1 < static_cast<int>(AccNodes + 0.5))
			SubBEnd++;
		SubB.push_back(make_pair(SubBBegin, SubBEnd));
		AccNodes += NodesPerDiap;
		SubBBegin = SubBEnd + 1;
		SubBEnd = SubBBegin;
	}
	TestSubB();
}

// test
void Diaps::TestSubB(){
	int Count = 0;
	if (SubB.size() != BaseLen)
		Error("Diaps::TestSubB", "Subborders count != BaseLen");
	for (auto it = SubB.begin(); it != SubB.end(); it++){
		if (it - SubB.begin() == 0 && it->first != L)
			Error("Diaps::TestSubB", "Left border of 1st subdiapason != L");
		if (it - SubB.begin() == SubB.size() && it->second != R)
			Error("Diaps::TestSubB", "Left border of 1st subdiapason != L");
		if (it - SubB.begin() != 0 && it->first != (it-1)->second + 1)
			Error("Diaps::TestSubB", "Left border of current subdiapason != right border of previous diapason + 1");
		if (it->first > it->second)
			Error("Diaps::TestSubB", "B.first > B.second");
		Count += it->second - it->first + 1;
	}
	if (Count != Len)
		Error("Diaps::TestSubB", "Count != Len");
}

int Diaps::GetSubBIndex(int Deg){
	if (SubB.size() == 0)
		Error("Diaps::GetSubBIndex", "SubB size == 0");
	if (Deg < SubB[0].first || Deg > SubB[SubB.size()-1].second)
		Error("Diaps::GetSubBIndex", "Deg is out of range");
	for (int i = 0; i < SubB.size(); i++){
		if (Deg >= SubB[i].first && Deg <= SubB[i].second)
			return i;
	}
	Error("Diaps::GetSubBIndex", "Deg not found");
}

void Diaps::SetProb(vector<double> P){
	if (P.size() != SubB.size())
		Error("Diaps::SetProb", "P.size() != SubB.size()");
	for (size_t i = 0; i < P.size(); i++)
		Prob.push_back(P[i]);
}

int Diaps::AddStrat(int I, int N){
	/*if (I == Index)
		Error("Diaps::AddStrat", "Strategy (x, (x,y) is not allowed");*/
	/*if (Nodes > 0 && I > Index)
		Error("Diaps::AddStrat", "Strategy is not allowed");
	if (Nodes < 0 && I < Index)
		Error("Diaps::AddStrat", "Strategy is not allowed");*/
	int NodesToStrat = StratNodes + N <= abs(Nodes) ? N : abs(Nodes) - StratNodes;
	Strat.push_back(make_pair(I, NodesToStrat));
	StratNodes += NodesToStrat;
	return NodesToStrat;
}

// get random degree with account of probabilities
int Diaps::GetRndDeg(TRnd& Rnd){
	double RndVal = Rnd.GetUniDev();
	int SubBIndex = 0;
	while (Prob[SubBIndex] <= RndVal)
		SubBIndex++;
	if (SubBIndex > SubB.size())
		Error("Diaps::GetRndDeg", "SubBIndex out of range");
	int RndDeg = Rnd.GetUniDev() * (SubB[SubBIndex].second - SubB[SubBIndex].first) + SubB[SubBIndex].first;
	return RndDeg;
}

// Priority = Sum_SubB: AvgDeg * Prob * Nodes
double Diaps::GetPriority(){
	double P = 0;
	if (SubB.size() != Prob.size())
		Error("Diaps::GetPriority", "SubB.size() != Prob.size()");
	for (size_t i = 0; i < SubB.size(); i++){
		double Deg = 0;
		for (size_t j = SubB[i].first; j <= SubB[i].second; j++)
			Deg += j;
		Deg /= SubB[i].second - SubB[i].first + 1;
		P += Deg * Prob[i] * abs(Nodes);
	}
	return P;
}

// add strat nodes
void Diaps::AddStratNodes(int S){
	if (StratNodes + S > Nodes)
		Error("Diaps::AddStratNodes", "StratNodes + S > Nodes");
	StratNodes += S;
}
Diaps::~Diaps(void)
{
}
