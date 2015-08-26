#pragma once
// to work with diapasons
class Diaps
{
	// index of diapason
	int Index;
	// base length of diapason
	int BaseLen;
	// real length of diapason
	int Len;
	// borders of diapason
	pair<int, int> Borders;
	// subborders according to the relation of BaseLen to Len
	vector<pair<int, int>> SubB;
	// probabilities of subdiapasons
	vector<double> Prob;
	// strategy (the most important strategy is at the end of the vector)
	// (index of another diapasone; nodes count to be removed to another diapasone)
	vector<pair<int, int>> Strat;
	// number of nodes which are used in strategies
	int StratNodes;
	// references
	int& L, &R;
	// nodes to add to this diapason (if negative, nodes should be deleted)
	int Nodes;
	// set subborders
	void SetSubB();
	// test
	void TestSubB();
public:
	Diaps(int I, pair<int, int> B, int BL);
	void SetNodes(int N);
	void SetProb(vector<double> P);
	int Length() {return Len;}
	int GetSubBIndex(int Deg);
	int GetL(){return L;}
	int GetR(){return R;}
	int GetIndex() {return Index;}
	int GetNodes(){return Nodes;}
	int GetStratNodes(){return StratNodes;}
	// returning value: number of nodes actually added to a strategy
	int AddStrat(int I, int N);
	// get random degree with account of probabilities
	int GetRndDeg(TRnd& Rnd);
	double GetPriority();
	// are strategies found completely
	bool IsStratFound() {return abs(Nodes) == StratNodes; }
	// get nodes count available for strategies
	int GetFreeNodes() {if (Nodes >=0) return Nodes - StratNodes; return Nodes + StratNodes;}
	// add strat nodes
	void AddStratNodes(int S);
	// print node info
	void PrintInfo(ofstream& F);
	// print strategies
	void PrintStrategies(ofstream& F);
	~Diaps(void);
};

