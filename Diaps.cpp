#include "StdAfx.h"
#include "Diaps.h"
#include "Error.h"
#include <iostream>

Diaps::Diaps(int I, pair<int, int> B, int BL) {
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
	BaseLen = BL;
	Len = Borders.second - Borders.first + 1;
	SetSubB();
	StratNodes = 0;
	ClusterClr();
}

void Diaps::SetNodes(int N){
	if (N == 0)
		Error("Diaps::SetNodes", "N == 0");
	Nodes = N;
}

// set subborders
void Diaps::SetSubB(){
	int &L = Borders.first, &R = Borders.second;
	double NodesPerDiap = static_cast<double>(Len) / BaseLen;
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
	int &L = Borders.first, &R = Borders.second;
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
	for (size_t i = 0; i < P.size(); i++){
		if (P[i] < 0 || P[i] > 1)
			Error("Diaps::SetProb", "P[i] < 0 || P[i] > 1");
		Prob.push_back(P[i]);
	}
}

int Diaps::AddStrat(int I, int N){
	/*if (I == Index)
		Error("Diaps::AddStrat", "Strategy (x, (x,y) is not allowed");*/
	/*if (Nodes > 0 && I > Index)
		Error("Diaps::AddStrat", "Strategy is not allowed");
	if (Nodes < 0 && I < Index)
		Error("Diaps::AddStrat", "Strategy is not allowed");*/
	if (N == 0)
		Error("Diaps::AddStrat", "Empty strategy");
	int NodesToStrat = StratNodes + N <= abs(Nodes) ? N : abs(Nodes) - StratNodes;
	if (NodesToStrat == 0)
		Error("Diaps::AddStrat", "Empty strategy");
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
	int RndDeg = Rnd.GetUniDev() * (SubB[SubBIndex].second - SubB[SubBIndex].first) + SubB[SubBIndex].first + 0.5;

	// DEBUG

	if (RndDeg < Borders.first || RndDeg > Borders.second)
		Error("Diaps::GetRndDeg", "RndDeg is out of range");

	// END DEBUG

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

// get index of neighbour diapason
int Diaps::GetNeighb(){
	if (Strat.size() == 0)
		Error("Diaps::GetNeighb", "There are no strategies available");
	return Strat[0].first;
}

// get number of nodes assigned to neighbour diapason
int Diaps::GetNeighbNodes(){
	if (Strat.size() == 0)
		Error("Diaps::GetNeighbNodes", "There are no strategies available");
	return Strat[0].second;
}

// check if cluster is full
bool Diaps::IsClusterFull(){
	int& ReqDeg = Cluster.first.first,
		&InitDeg = Cluster.first.second.second,
		&FreeNodes = Cluster.second.first;
	if (InitDeg + FreeNodes >= ReqDeg)
		return true;
	return false;
}

// delete first node of the cluster's list
void Diaps::DelClusterFirst(){
	if (GetClusterSize() == 0)
		Error("Diaps::DelClusterFirst", "Cluster size = 0");
	Cluster.second.second.erase(Cluster.second.second.begin());
	Cluster.second.first = -1;
	if (Cluster.second.second.size() == 0)
		ClusterClr();
}

// add node to cluster
void Diaps::AddToCluster(int Node, bool HasEdge){
	Cluster.second.second.push_back(Node);
	if (!HasEdge)
		Cluster.second.first++;
}

// clear the cluster
void Diaps::ClusterClr(){
	Cluster.first.first = -1;
	Cluster.first.second.first = -1;
	Cluster.first.second.second = -1;
	Cluster.second.first = -1;
	Cluster.second.second.clear();
}

// decrease nodes count for strategy
bool Diaps::DecreaseStratN(){
	if (Strat.size() == 0)
		Error("Diaps::DecreaseStratN", "Strat.size() == 0");
	Strat[0].second--;
	if (Strat[0].second == 0){
		Strat.erase(Strat.begin());
		if (Strat.size() == 0)
			ClusterClr();
		return true;
	}
	return false;
}

// reset cluster from first node in the list
void Diaps::ResetCluster(int ReqDeg, int CInitDeg, int TargNCount){
	if (Cluster.second.second.size() == 0)
		Error("Diaps::ResetCluster", "Cluster.size() == 0");
	Cluster.first.first = ReqDeg;
	Cluster.first.second.first = Cluster.second.second[0];
	Cluster.first.second.second = CInitDeg;
	Cluster.second.first = TargNCount;
	Cluster.second.second.erase(Cluster.second.second.begin());
}

// print node info
void Diaps::PrintInfo(ofstream& F){
	F << "Index: " << Index << "[" << Borders.first << ";" << Borders.second << "]" << " To add nodes: " << Nodes << endl;
	F << "Subborders: ";
	for (size_t i = 0; i < SubB.size(); i++) 
		F << "[" << SubB[i].first << ";" << SubB[i].second << "]" << " ";
	F << endl;
	F << "Prob: ";
	for (size_t i = 0; i < Prob.size(); i++) 
		F << Prob[i] << " ";
	F << endl;
}
// print strategies
void Diaps::PrintStrategies(ofstream& F){
	F << "Strat nodes (final): " << StratNodes << " Nodes remained: " << abs(Nodes) - StratNodes << endl;
	for (size_t i = 0; i < Strat.size(); i++) 
		F << "(" << Strat[i].first << "," << Strat[i].second << ")" << " ";
	F << endl;
}

// print cluster info
void Diaps::PrintClusterInfo(){
	std::cout << "Node: " << Cluster.first.second.first << " init deg: " << Cluster.first.second.second <<
		" req deg: " << Cluster.first.first << " cluster size: " << Cluster.second.second.size() <<
		" free nodes: " << Cluster.second.first << " nodes in the cluster: ";
	for (size_t i = 0; i < GetClusterSize(); i++)
		std::cout << Cluster.second.second[i] << " ";
	std::cout << endl;
}

Diaps::~Diaps(void)
{
}
