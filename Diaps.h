#pragma once
// to work with diapasons
class Diaps : public BaseDiap
{
	// strategy (the most important strategy is at the beginning of the vector)
	// (index of another diapasone; nodes count to be removed to another diapasone)
	vector<pair<int, int>> Strat;
	// number of nodes which are used in strategies
	int StratNodes;
	// nodes to add to this diapason (if negative, nodes should be deleted)
	int Nodes;
	// clusters: ((ReqDeg, <NodeToCluster, InitialDeg>), (AppropriateNodesCount,vector<Nodes indexes>)) for Strat[0] (if it exists) for DMinus
	pair<pair<int, pair<int,int>>, pair<int,vector<int>>> Cluster;
	// test probability
	void TestProb();
public:
	Diaps(int I, pair<int, int> B, int BL, double MK, double Prev, double W);
	void SetSubB(const TFltPrV& KronDeg);
	void SetProb(vector<double> P);
	void ResetCumProb();
	void SetNParts(vector<double> NP);
	void SetNodes(int N);
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
	// check if diapason has strategies
	bool HasStrat() { if (Strat.size() != 0) return true; return false; }
	// decrease nodes count for strategy
	bool DecreaseStratN();
	// return "sign" of diapason
	bool Sign() { if (Nodes > 0) return true; return false; }
	// get index of neighbour diapason
	int GetNeighb();
	// get number of nodes assigned to neighbour diapason
	int GetNeighbNodes();
	// check if cluster is empty
	bool IsClusterEmpty() { if (Cluster.first.second.first == -1) return true; return false; }
	// get cluster required degree
	int GetClusterReqDeg() {return Cluster.first.first;}
	// get cluster initial degree
	int GetClusterInitDeg() {return Cluster.first.second.second;}
	// get cluster node index
	int GetClusterNode() {return Cluster.first.second.first;}
	// get target nodes count
	int GetTargNCount() {return Cluster.second.first;}
	// get cluster size
	size_t GetClusterSize() {return Cluster.second.second.size();}
	// delete first node of the cluster's list
	void DelClusterFirst();
	// get cluster
	vector<int>& GetCluster() {return Cluster.second.second;}
	// check if cluster is full
	bool IsClusterFull();
	// clear the cluster
	void ClusterClr();
	// set cluster degree
	void SetCluster(int ReqDeg, int Node, int InitDeg) {Cluster.first.first = ReqDeg; Cluster.first.second.first = Node; Cluster.first.second.second = InitDeg; Cluster.second.first = 0;}
	// reset cluster from first node in the list
	void ResetCluster(int ReqDeg, int CInitDeg, int TargNCount);
	// decrease CInitDeg by 1
	void DecreaseCInitDeg() {Cluster.first.second.second--;}
	// add node to cluster
	void AddToCluster(int Node, bool HasEdge);
	// print diapason info
	virtual void PrintInfo(ofstream& F);
	// print strategies
	void PrintStrategies(ofstream& F);
	// print cluster info
	void PrintClusterInfo();
	~Diaps(void);
};

