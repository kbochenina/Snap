#include "StdAfx.h"
#include "Diaps.h"
#include "Error.h"
#include <iostream>

Diaps::Diaps(int I, pair<int, int> B, int BL, double MK, double Prev, double W) : BaseDiap(I, B, BL, MK, Prev, W) {
	StratNodes = 0;
	ClusterClr();
}

void Diaps::SetNodes(int N){
	if (N == 0)
		Error("Diaps::SetNodes", "N == 0");
	Nodes = N;
}

void Diaps::SetSubB(const TFltPrV& KronDeg){
	if (Prob.size() != BaseLen)
		Error("Diaps::SetSubB", "Wrong size of Prob");
	if (SubB.size() != 0)
		SubB.clear();
	int& L = Borders.first, &R = Borders.second;
	int DiapNCount = 0;
	
	vector<int> KronNodes;
	for (int i = 0; i < Len; i++)
		KronNodes.push_back(0);

	for (size_t i = 0; i < KronDeg.Len(); i++){
		int Deg = KronDeg[i].Val1;
		if (Deg >= L && Deg <= R){
			int NCount = KronDeg[i].Val2;
			KronNodes[Deg-L] = NCount;
			DiapNCount += NCount;
		}
	}

	if (DiapNCount==0){
		SubB.push_back(make_pair(Borders.first, Borders.second));
		for(size_t i = 1; i < BaseLen; i++)
			SubB.push_back(make_pair(-1,-1));
		return;
	}


	if (KronNodes.size() != Len)
		Error("Diaps::SetSubB", "Wrong size of KronNodes");

	int KronNodesInd = 0;
	int DiapBegin = L;
	bool EndReached = false;
	int ZeroDiapsCount = 0;
	for (int i = 0; i < BaseLen; i++){
		double SubDNCount = NParts[i] * DiapNCount;
		int AccNodes = 0;
		int DiapEnd = 0;
		
		if (SubDNCount == 0 || EndReached){
			++ZeroDiapsCount;
			continue;
		}

		while (1){
			if (KronNodesInd == Len-1) {
				DiapEnd = L + KronNodesInd;
				break;
			}
			AccNodes += KronNodes[KronNodesInd];
			// if next diapason will cause "overflow"
			int NextAccNodes = AccNodes + KronNodes[KronNodesInd+1];
			if (NextAccNodes >= SubDNCount){
				// check if current or next degree is closer to desired right border
				bool IsCurr = (abs(AccNodes-SubDNCount) < abs(NextAccNodes-SubDNCount)) ? true : false;
				if (IsCurr)
					DiapEnd = L + KronNodesInd;
				else{
					DiapEnd = L + KronNodesInd + 1;
					KronNodesInd++;
				}
				KronNodesInd++;
				break;
			}
			KronNodesInd++;
		};
				
		if (DiapBegin < L && DiapEnd > R)
			Error("Diaps::SetSubB", "Attempt to add wrong subborders");
		SubB.push_back(make_pair(DiapBegin, DiapEnd));
		if (DiapEnd == R) EndReached = true;
		DiapBegin = DiapEnd + 1;
	}

	for (int i = 0; i < ZeroDiapsCount; ++i)
		SubB.push_back(make_pair(-1,-1));

	// if the end of the diapason is not equal to actual end, correct it
	int i = 0;
	for (; i < BaseLen; ++i){
		if (SubB[i].first == -1) break;
	}
	// i-1 is the index of last actual diapason
	if (SubB[i-1].second != R)
		SubB[i-1].second = R;

	if (SubB.size() != BaseLen)
		Error("Diaps::SetSubB", "Wrong SubB size");
}


void Diaps::SetProb(vector<double> P){
	/*if (P.size() != SubB.size())
	Error("Diaps::SetProb", "P.size() != SubB.size()");*/
	for (size_t i = 0; i < P.size(); i++){
		if (P[i] < 0 || P[i] > 1)
			Error("Diaps::SetProb", "P[i] < 0 || P[i] > 1");
		Prob.push_back(P[i]);
	}
}

void Diaps::SetNParts(vector<double> NP){
	for (size_t i = 0; i < NP.size(); i++){
		if (NP[i] < 0 || NP[i] > 1)
			Error("Diaps::SetProb", "NP[i] < 0 || NP[i] > 1");
		NParts.push_back(NP[i]);
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

void Diaps::ResetCumProb(){
	if (SubB.size() != BaseLen || Prob.size() != BaseLen)
		Error("Diaps::ResetCumProb", "Inconsistent size of SubB and/or Prob");
	double CumProbPrev = Prob[0];
	for (size_t i = 1; i < BaseLen; ++i){
		if (SubB[i].first == -1 && CumProbPrev != 1){
			// if there is only first diapason
			if (CumProbPrev == 0 && i == 1){
				for (size_t j = 0; j < BaseLen; j++)
					Prob[j] = 1;
				break;
			}
			/*if (CumProbPrev == 0)
				Error("Diaps::ResetCumProb()", "CumProbPrev == 0");*/
			double Least = 1 - CumProbPrev;
			for (size_t j = 0; j < i; ++j){
				if (CumProbPrev == 0)
					Prob[j] = (j+1) / static_cast<double>(i-1); 
				else 
					Prob[j] += Prob[j]/CumProbPrev * Least;
				if (Prob[j]>1)
					Prob[j] = 1;
			}
			for (size_t j = i; j < BaseLen; ++j)
				Prob[j] = 1;
			break;
		}
		else
			CumProbPrev = Prob[i];
	}
	TestProb();
}


void Diaps::TestProb(){
	for (size_t i = 0; i < BaseLen; ++i){
		if (Prob[i] == 1 && i != BaseLen-1){
			if (SubB[i].first == -1)
				Error("Diaps::TestProb()", "Inconsistent probability value");
			// all next probs should be = 1, 
			for (size_t j = i+1; j < BaseLen; ++j)
				if (Prob[j] != 1)
					Error("Diaps::TestProb()", "Inconsistent probability value");
			break;

		}
	}
}
// get random degree with account of probabilities
int Diaps::GetRndDeg(TRnd& Rnd){
	double RndVal = Rnd.GetUniDev();
	int SubBIndex = 0;
	while (Prob[SubBIndex] <= RndVal){
		SubBIndex++;
		if (SubBIndex == Prob.size())
			Error("Diaps::GetRndDeg", "SubBIndex is out of range");
	}
	if (SubBIndex > SubB.size())
		Error("Diaps::GetRndDeg", "SubBIndex is out of range");
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
	for (size_t i = 0; i < SubB.size(); ++i){
		double Deg = 0;
		// for empty diapasones
		if (SubB[i].first == -1) continue;
		for (size_t j = SubB[i].first; j <= SubB[i].second; ++j)
			Deg += j;
		Deg /= SubB[i].second - SubB[i].first + 1;
		if (i==0)
			P += Deg * Prob[i] * abs(Nodes);
		else P += Deg * (Prob[i]-Prob[i-1]) * abs(Nodes);
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
