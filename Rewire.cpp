#include "stdafx.h"

bool GetDiap(const vector<Diap>& Diaps, const TFlt& DegPart, bool& DiapSign, TInt& DiapIndex){
	DiapIndex = 0;
	for (auto DiapIt = Diaps.begin(); DiapIt != Diaps.end(); DiapIt++){
		if (DiapIt->first.Val1 > DegPart) 
			return false;
		else if (DiapIt->first.Val1 <= DegPart && DiapIt->first.Val2 >= DegPart){
			DiapSign = DiapIt->first.Val2 > 0 ? true : false;
			return true;
		}
		DiapIndex++;
	}
	return false;
}

int GetRandDeg(TRnd& Rnd, const Diap& Borders, const int DegMin, const int DegMax){
	TInt DiapBegin = (DegMax - DegMin) * Borders.first.Val1 + DegMin,
		DiapEnd = (DegMax - DegMin) * Borders.first.Val2 + DegMin;
	return Rnd.GetUniDev() * (DiapEnd - DiapBegin) + DiapBegin;
}

// get appropriate number of nodes for each diapasone and its average degree
void GetDiaps(vector<Diaps>& DPlus, vector<Diaps>& DMinus, const vector<Diap>& SmoothedDiaps, const TFltPrV& KronDeg, const TInt& DegMin, const TInt& DegMax, const vector<int>& Prev){
	for (auto DiapIt = SmoothedDiaps.begin(); DiapIt != SmoothedDiaps.end(); DiapIt++){
		int Index = DiapIt - SmoothedDiaps.begin();
				

		int DiapBegin = static_cast<int>((DegMax - DegMin + 1) * DiapIt->first.Val1 + DegMin + 0.5),
			DiapEnd = static_cast<int>((DegMax - DegMin + 1) * DiapIt->first.Val2 + DegMin + 0.5);
		// check
		if (find(Prev.begin(), Prev.end(), DiapIt-SmoothedDiaps.begin()) != Prev.end()){
			// if in the sample this diapason starts after previous
			DiapBegin = (DegMax - DegMin + 1) * (DiapIt-1)->first.Val2 + DegMin + 0.5 + 1;
		}

		//// HACK!
		//if (Index == 1){
		//	DiapBegin = 3;
		//	DiapEnd = 11;
		//}

		//printf("DegBegin: %d DegEnd: %d\n", DiapBegin, DiapEnd);
		pair<int,int> Borders(DiapBegin, DiapEnd);
		int BaseLen = DiapIt->second.Len();

		// DEBUG (FOR MODEL_SIZE == KRON_SIZE)

		/*if (DiapIt->second.Len() != DiapEnd - DiapBegin + 1)
			Error("GetDiaps", "RelDiff count is not equal to nodes count");*/

		// DEBUG END

		Diaps NewDiap(Index, Borders, BaseLen);

		int NodesCount = 0;
		double NodesToAdd = 0, RelDiffSum = 0;
		int KronInd = 0, Count = 0, KDeg = 0;

		for (size_t Deg = DegMin; Deg <= DegMax; Deg++){
			if (KronInd < KronDeg.Len())
				KDeg = KronDeg[KronInd].Val1;
			if (KDeg == Deg){
				Count = KronDeg[KronInd].Val2;
				KronInd++;
			}
			else
				Count = 0;
			if (Deg > DiapEnd) break;
			if (Deg < DiapBegin) continue;
			int SubBIndex = NewDiap.GetSubBIndex(Deg);
			if (SubBIndex < 0 || SubBIndex > DiapIt->second.Len() - 1)
				Error("GetDiaps", "Incorrect index of subdiapason");
			// add nodes with account of relative difference
			double RelDiff = DiapIt->second[SubBIndex];
			if (abs(RelDiff) != 1 && RelDiff != 0) {
				printf("%3.2f\n", log10(static_cast<double>(Count+1)));
				NodesToAdd += pow(10, log10(static_cast<double>(Count+1)) * RelDiff) - Count - 1;
			}
			else 
				if (RelDiff == 1.0){ 
					if (Count == 0) 
						NodesToAdd += 1;
				}
			else if (RelDiff == -1)
				NodesToAdd -= Count;
			NodesCount += Count;
		}
		//printf("\n");
		int RoundNodesToAdd = static_cast<int>(NodesToAdd + 0.5);
		if (RoundNodesToAdd != 0){
			NewDiap.SetNodes(RoundNodesToAdd);
			// calculate sum of RelDiff for diapason
			for (size_t i = 0; i < DiapIt->second.Len(); i++)
				RelDiffSum += abs(DiapIt->second[i]);
			// calculate accumulated probabilities of subdiapasons
			vector<double> Prob;
			double RelDiffAcc = 0;
			for (size_t i = 0; i < DiapIt->second.Len(); i++){
				double RelDiff = abs(DiapIt->second[i]);
				//printf("%3.2f ", RelDiff);
				RelDiffAcc += RelDiff;
				double P = RelDiffAcc / RelDiffSum;
				if (RelDiffSum == 0)
					Error("GetDiaps", "Wrong value of P");
				if (P < 0 || P > 1)
					Error("GetDiaps", "Wrong value of P");
				Prob.push_back(P);
			}
			//printf("\n");
			NewDiap.SetProb(Prob);
			//printf("DiapBegin: %d DiapEnd: %d Nodes count: %d Nodes to add: %d Res: %d\n", NewDiap.GetL(),  NewDiap.GetR(), NodesCount, 
			//NewDiap.GetNodes(), NodesCount + NewDiap.GetNodes());
			if (NewDiap.GetNodes() > 0) DPlus.push_back(NewDiap);
			else DMinus.push_back(NewDiap);
		}
	}
}

bool PrVecComp(const pair<pair<bool,int>, double>& Pr1, const pair<pair<bool,int>, double>& Pr2){
	return Pr1.second > Pr2.second;
}

void GetPriorities(vector<pair<bool,int>>&Pr, vector<Diaps>& DPlus, vector<Diaps>& DMinus){
	vector<pair<pair<bool,int>, double>> PrVec;
	for (size_t i = 0; i < DPlus.size(); i++)
		PrVec.push_back(make_pair(make_pair(true, i), DPlus[i].GetPriority()));
	for (size_t i = 0; i < DMinus.size(); i++)
		PrVec.push_back(make_pair(make_pair(false, i), DMinus[i].GetPriority()));

	sort(PrVec.begin(), PrVec.end(), PrVecComp);
	for (size_t i = 0; i < PrVec.size(); i++)
		Pr.push_back(PrVec[i].first);
}

// get rewire stratergies
int GetRewireStrategies(vector<Diaps>& DPlus, vector<Diaps>& DMinus){
	
	TRnd Rnd;

	// vector of priorities
	// (true-DPlus, false-DMinus, Index)
	vector<pair<bool,int>> Pr;
	GetPriorities(Pr, DPlus, DMinus);
	
	int EdgesToAdd = 0;
	
	for (size_t i = 0; i < Pr.size(); i++){
		bool CurrTrue = Pr[i].first == true;
		int CurrLocalInd = Pr[i].second;
		Diaps& Curr =  CurrTrue ? DPlus[CurrLocalInd] : DMinus[CurrLocalInd];
		if (Curr.IsStratFound()) continue;
		for (size_t j = 0; j < Pr.size(); j++){
			bool NeighbTrue = Pr[j].first == true;
			if (CurrTrue == NeighbTrue) continue;
			int NeighbLocalInd = Pr[j].second;
			Diaps& Neighb =  NeighbTrue ? DPlus[NeighbLocalInd] : DMinus[NeighbLocalInd];
			if (Neighb.IsStratFound()) continue;
			
			bool CurrIndHighest = Curr.GetIndex() > Neighb.GetIndex() ? true : false;
			int LessNodes = CurrIndHighest ? Neighb.GetFreeNodes() : Curr.GetFreeNodes();
			int BiggNodes = CurrIndHighest ? Curr.GetFreeNodes() : Neighb.GetFreeNodes();
			int LessDeg, BiggDeg;

			// if we can delete edges for nodes with biggest degree to obtain nodes with smallest degree
			if (LessNodes > 0 && BiggNodes < 0){
				int Nodes = abs(BiggNodes) >= LessNodes ? LessNodes : abs(BiggNodes);
				if (CurrIndHighest){
					Curr.AddStrat(NeighbLocalInd, Nodes);
					Neighb.AddStratNodes(Nodes);
				}
				else {
					Neighb.AddStrat(CurrLocalInd, Nodes);
					Curr.AddStratNodes(Nodes);
				}

				int ApproxEdges = 0;
				for (int i = 0; i < Nodes; i++){
					BiggDeg = CurrIndHighest ? Curr.GetRndDeg(Rnd) : Neighb.GetRndDeg(Rnd);
					LessDeg = CurrIndHighest ? Neighb.GetRndDeg(Rnd) : Curr.GetRndDeg(Rnd);
					ApproxEdges += BiggDeg - LessDeg;
				}
				EdgesToAdd += ApproxEdges * 2;
			}
			// if we can add edges to nodes with smaller degree to obtain nodes with higher degree
			else if (LessNodes < 0 && BiggNodes > 0) {
				int Nodes = 0, ClusteredEdges = 0, Ind = 0;
				TIntPrV AddEdges; 
				while (1)
				{
					BiggDeg = CurrIndHighest ? Curr.GetRndDeg(Rnd) : Neighb.GetRndDeg(Rnd);
					LessDeg = CurrIndHighest ? Neighb.GetRndDeg(Rnd) : Curr.GetRndDeg(Rnd);
					Nodes += 1;

					for (size_t d = 0; d < AddEdges.Len(); d++){
						// if period is finished, delete it
						if (AddEdges[d].Val2 < Ind)
							AddEdges.Del(d);
					}
					int AddECount = AddEdges.Len();
					if (LessDeg + AddECount < BiggDeg){
						int E = BiggDeg - LessDeg - AddECount;
						ClusteredEdges += E;
						// from the next iteration to the (next + edges - 1) iteration inclusively
						TIntPr Diap(Ind + 1, Ind + E);
						AddEdges.Add(Diap);
					}
					Ind++;
					if (Nodes  >= abs(LessNodes) || Ind >= BiggNodes)
						break;
				}
								
				if (Ind > 0){
					if (CurrIndHighest){
						Neighb.AddStrat(CurrLocalInd, Nodes);
						Curr.AddStratNodes(Ind);
					}
					else {
						Curr.AddStrat(NeighbLocalInd, Nodes);
						Neighb.AddStratNodes(Ind);
					}
					EdgesToAdd -= 2 * ClusteredEdges;
				}
			}
		}
	}

	printf("Edges to add: %d\n", EdgesToAdd);
	//PrintDegDistr(DiapNodes, "DiapNodesAfter.tab");
	return EdgesToAdd;
}

// get diap index
bool GetDiap(int Deg, vector<Diaps>& DPlus, vector<Diaps>& DMinus, int& DiapIndex, bool& IsDPlus){
	// DPlus and DMinus are sorted in increasing order
	for (size_t i = 0; i < DPlus.size(); i++){
		if (DPlus[i].IsDegInDiap(Deg)){
			DiapIndex = i;
			IsDPlus = true;
			return true;
		}
	}
	for (size_t i = 0; i < DMinus.size(); i++){
		if (DMinus[i].IsDegInDiap(Deg)){
			DiapIndex = i;
			IsDPlus = false;
			return true;
		}
	}
	return false;
}

// get random number from diap
int GetRndDeg(TRnd& Rnd, const TIntPr& Borders){
	return Rnd.GetUniDev()*(Borders.Val2 - Borders.Val1) + Borders.Val1;
}

// get random number from diap with account of possibilities of subdiaps
int GetRndDeg(TRnd&Rnd, const TIntPr& Borders, const TFltV& DiapPossib){
	double RndVal = Rnd.GetUniDev();
	int DiapNum = 0;
	while (RndVal > DiapPossib[DiapNum]) ++DiapNum;
	int DiapSize = DiapPossib.Len();
	int RealDiapSize = Borders.GetVal2() - Borders.GetVal1() + 1;
	int NodesPerRealDiap = static_cast<int>(RealDiapSize / DiapSize);
	int DiapBegin = Borders.GetVal1() + DiapNum * NodesPerRealDiap;
	int DiapEnd = DiapBegin + NodesPerRealDiap - 1;
	if (DiapBegin < Borders.GetVal1() || DiapEnd > Borders.GetVal2())
		Error("GetRndDeg", "Wrong borders");
	TIntPr NewBorders(DiapBegin, DiapEnd);
	return GetRndDeg(Rnd, NewBorders);
}
	
vector<int>& FindCluster(ClusterMap& Clusters, int DiapIndex, int ReqDiap, int& ReqDeg){
	map<pair<int,int>, vector<int>>& ClustersDiap = Clusters[DiapIndex];
	for (map<pair<int,int>, vector<int>>::iterator CIt = ClustersDiap.begin(); CIt != ClustersDiap.end(); CIt++){
		if (CIt->first.first == ReqDiap){
			ReqDeg = CIt->first.second;
			return CIt->second;
		}
	}
	// create cluster and add reference to its vector of nodes to Cluster
	return ClustersDiap.insert(ClustersDiap.begin(), make_pair(make_pair(ReqDiap, ReqDeg), vector<int>()))->second;
}

void SetReqDeg(ClusterMap& Clusters, int DiapIndex, int ReqDiap, int ReqDeg, vector<int>& ClusterOld){
	map<pair<int,int>, vector<int>>& ClustersDiap = Clusters[DiapIndex];
	pair<int, int> KeyPair;
	vector <int> Cluster;
	map<pair<int,int>, vector<int>>::iterator EraseIt = ClustersDiap.end();
	for (map<pair<int,int>, vector<int>>::iterator CIt = ClustersDiap.begin(); CIt != ClustersDiap.end(); CIt++){
		if (CIt->first.first == ReqDiap){
			//CIt->first.second = ReqDeg;
			KeyPair = CIt->first;
			if (KeyPair.second == ReqDeg) return;
			KeyPair.second = ReqDeg;
			for (size_t i = 0; i < CIt->second.size(); i++) 
				Cluster.push_back(CIt->second[i]);
			EraseIt = CIt;
		}
	}
	// erase current 
	if (EraseIt == ClustersDiap.end())
		Error("SetReqDeg", "Cannot find diapasone");
	if (Cluster.size() == 0){
		system("pause");
	}
	ClustersDiap.insert(ClustersDiap.begin(), make_pair(KeyPair, Cluster));
	ClusterOld = Cluster;
	ClustersDiap.erase(EraseIt);
	
}

// add random edge
bool AddRndEdge(TRnd& Rnd, PNGraph&Kron, int Node, int DegMax){
	int Attempts = Kron->GetEdges();
	while (Attempts > 0){
		int Neighb = Rnd.GetUniDev() * (Kron->GetNodes() - 1);
		if (Kron->IsEdge(Node, Neighb) || Node == Neighb || Kron->GetNI(Neighb).GetInDeg() == DegMax){
			Attempts--; continue;
		}
		Kron->AddEdge(Node, Neighb);
		Kron->AddEdge(Neighb, Node);
		return true;
	}
	return false;
}

void Rewire(PNGraph& Kron, vector<Diaps>& DPlus, vector<Diaps>& DMinus, const int DegMin, const int DegMax){
	TRnd Rnd;
	int Add = 0, Del = 0;
	int BasicEdgesCount = Kron->GetEdges();
	double SecAdd = 0, SecDel = 0, TimeToFindDel = 0;
	TExeTm execTime;

	for (auto NodeIt = Kron->BegNI(); NodeIt != Kron->EndNI(); NodeIt++){
		//cout << "Add " << Add << " Del " << Del << endl;
		bool CanAdd = true, CanDel = true;
		// TEST CONDITION
		if (abs(Add-Del)/static_cast<double>(BasicEdgesCount) > 0.01){ 
			if (Add > Del) CanAdd = false;
			else CanDel = false;
		}
		
		// END TEST

		int Node = NodeIt.GetId();
		int Deg = NodeIt.GetInDeg();
		int DiapIndex = 0; 
		bool IsDPlus;
		bool IsDegInDiap = GetDiap(Deg, DPlus, DMinus, DiapIndex, IsDPlus);
		if (!IsDegInDiap) continue;
		Diaps& D = IsDPlus ? DPlus[DiapIndex] : DMinus[DiapIndex];
		if (!D.HasStrat() || IsDPlus) continue;
		int NeighbInd = D.GetNeighb();
		if (NeighbInd < 0 || NeighbInd > DPlus.size()-1)
			Error("Rewire", "Wrong neighbour index");
		Diaps& N = DPlus[NeighbInd];

		execTime.Tick();

		// if nodes are deleted from current diapason
		// and they are clustered to form a new value
		// and CanAdd == true
		if (IsDPlus == false && N.GetL() > Deg && CanAdd){
			
			// TEST
			if (N.GetL() - Deg > abs(D.GetNodes())){
				for (size_t i = 0; i < N.GetL() - Deg; i++){
					AddRndEdge(Rnd, Kron, Node, DegMax);
					Add+=2;
				}
				D.DecreaseStratN();
				continue;
			}
			// END TEST 

			int ReqDeg = 0, CNode = Node, CInitDeg = Deg;

			// cluster is empty if it has no CNode (value = -1)  
			if (D.IsClusterEmpty()){
				ReqDeg = N.GetRndDeg(Rnd);
				D.SetCluster(ReqDeg, CNode, CInitDeg);
			}
			else{
				ReqDeg = D.GetClusterReqDeg();
				CNode = D.GetClusterNode();
				CInitDeg = D.GetClusterInitDeg();
				bool HasDeg = Kron->IsEdge(CNode, Node);
				D.AddToCluster(Node, HasDeg);
			}
			
			//D.PrintClusterInfo();

			// while cluster is full
			while ( D.IsClusterFull() && !D.IsClusterEmpty() ){
				vector<int>& Cluster = D.GetCluster();
				vector<int> TargNodes;
				int TargNCount = D.GetClusterReqDeg() - D.GetClusterInitDeg();
				int Left = TargNCount;
				for (size_t i = 0; i < Cluster.size(); i++){
					if (!Kron->IsEdge(CNode, Cluster[i])){
						TargNodes.push_back(Cluster[i]);
						--Left;
						if (Left == 0)
							break;
					}
				}
				
				if (TargNodes.size() != TargNCount)
					Error("Rewire", "Incorrect size of TargNodes");
				
				//printf("Init degree: %d\n", Kron->GetNI(CNode).GetInDeg());

				for (size_t i = 0; i < TargNCount; i++){
					
					// DEBUG

					if (Kron->IsEdge(CNode, TargNodes[i]) || Kron->IsEdge(TargNodes[i], CNode))
						Error("Rewire", "Attempt to add existing edge");
					
					// DEBUG END

					Kron->AddEdge(CNode, TargNodes[i]);
					Kron->AddEdge(TargNodes[i], CNode);

					// DEBUG

					if (!Kron->IsEdge(CNode, TargNodes[i]))
						Error("Rewire", "Edge addition failed");

					// DEBUG END

					Add+=2;
				}
				
				// DEBUG

				int Edges = Kron->GetEdges();
				if (Edges != BasicEdgesCount + Add - Del)
					Error("Rewire", "Basic edges count != Edges + Add - Del");

				int NewCDeg = Kron->GetNI(CNode).GetInDeg();
				if (NewCDeg != ReqDeg)
					Error("Rewire", "Wrong degree after clustering");

				// DEBUG END

				ReqDeg = DPlus[D.GetNeighb()].GetRndDeg(Rnd);
				
				// while CInitDeg >= ReqDeg
				do{
					bool StratFin = D.DecreaseStratN();
					// if all strategies are completed
					if (StratFin && !D.HasStrat())	
						break;
					CNode = Cluster[0];
					CInitDeg = Kron->GetNI(CNode).GetInDeg();
					if (CInitDeg >= ReqDeg){
						D.DelClusterFirst();
					}
					else
						break;
				}
				while (!D.IsClusterEmpty());

				if (!D.IsClusterEmpty()){
					TargNCount = 0;
					// calculate TargNCount (indexing from 1!)
					for (size_t i = 1; i < Cluster.size(); i++){
						if (!Kron->IsEdge(CNode, Cluster[i])){
							TargNCount++;
						}
					}
					D.ResetCluster(ReqDeg, CInitDeg, TargNCount);
				}
			}
			continue;
		}

		SecAdd += execTime.GetSecs();

		if (CanDel == false || N.GetL() > Deg)
			continue;
		execTime.Tick();

		int ReqDeg = N.GetRndDeg(Rnd);
		//printf("Init degree: %d\n", Deg);
		int ToDel = Deg - ReqDeg;
		
		// DEBUG

		if (ToDel <= 0)
			Error("Rewire", "ToDel is non-positive");

		// END DEBUG
		TExeTm TimeToFind;
		int Attempts = 2 * Deg, EdgesDel = 0;
		while (Attempts != 0 && EdgesDel != ToDel){
			
			int NbInd = Rnd.GetUniDev() * (Deg - EdgesDel);
			int NNode = Kron->GetNI(Node).GetNbrNId(NbInd);
			int NDeg = Kron->GetNI(NNode).GetOutDeg();
			
			if (NDeg == DegMin ||
				NDeg == DegMax){
				Attempts--;
				continue;
			}
			
			

			Kron->DelEdge(Node, NNode, false);
			
			TimeToFind.Tick();

			// if the node deleted was the head of the cluster, we should decrease its CInitDeg by 1
			bool IsDegInDiap = GetDiap(NDeg, DPlus, DMinus, DiapIndex, IsDPlus);
			if (IsDegInDiap && !IsDPlus){
				Diaps& RDMin = DMinus[DiapIndex];
				if (RDMin.GetClusterNode() == NNode){
					RDMin.DecreaseCInitDeg();
					printf("Deg: %d CInitDeg:%d\n", Kron->GetNI(NNode).GetInDeg(), N.GetClusterInitDeg());
				}
			} 
			
			EdgesDel++;
			Del += 2;

			TimeToFindDel += TimeToFind.GetSecs();
			// DEBUG

			/*int Edges = Kron->GetEdges();
			if (Edges != BasicEdgesCount + Add - Del)
			Error("Rewire", "Basic edges count != Edges + Add - Del");*/

			

			// DEBUG END
		}
		//TimeToFindDel += TimeToFind.GetSecs();
		D.DecreaseStratN();
		SecDel += execTime.GetSecs();
	}
	
	cout <<  "Edges added " << Add << ", edges deleted " << Del << ", difference " << Add - Del << endl;
	cout <<  "Time of addition " << SecAdd << ", time of deletion " << SecDel <<  ", time to find del " << TimeToFindDel << endl;
	
}



// add missing or delete excess edges
void AddEdges(PNGraph&Kron, int Diff, int DegMin, int DegMax, int ModelEdges, vector<Diaps>& DMinus, vector<Diaps>& DPlus){
	TRnd Rnd;
	bool IsDPlus, IsInDiap;
	int DiapIndex;
	int E = 0;
	if (Diff < 0){
		while (E > Diff){
			int Node1 = Rnd.GetUniDev() * Kron->GetNodes(), Node2;
			int Deg1 = Kron->GetNI(Node1).GetOutDeg();
			IsInDiap = GetDiap(Deg1, DPlus, DMinus, DiapIndex, IsDPlus);
			if (Deg1 == DegMax )
				continue;
			while (1)
			{
				Node2 = Rnd.GetUniDev() * Kron->GetNodes();
				int Deg2 = Kron->GetNI(Node2).GetOutDeg();
				bool Node2DPlus;
				IsInDiap = GetDiap(Deg2, DPlus, DMinus, DiapIndex, Node2DPlus);
				if (Node1 == Node2 || Deg2 == DegMax || (IsInDiap && Node2DPlus == IsDPlus)) continue;
				break;
			}
			Kron->AddEdge(Node1, Node2);
			Kron->AddEdge(Node2, Node1);
			E -= 2;
		}
	}
	else {
		while (E < Diff){
			int Node1 = Rnd.GetUniDev() * Kron->GetNodes(), Node2;
			int Node1Deg = Kron->GetNI(Node1).GetOutDeg(), Node2Deg;
			if (Node1Deg == DegMin) // Node1Deg == DegMax || 
				continue;
			bool Node2Found = false;
			int Attempts = 0;
			while (Attempts != 2 * Node1Deg)
			{
				int Node2Ind = Node1Deg * Rnd.GetUniDev();
				Node2 = Kron->GetNI(Node1).GetNbrNId(Node2Ind);
				Node2Deg = Kron->GetNI(Node2).GetInDeg();
				if (Node1 == Node2 // || Node2Deg == DegMax 
					|| Node2Deg == DegMin) {
						Attempts++;
						continue;
				}
				Node2Found = true;
				break;
			}
			if (!Node2Found)
				continue;
			Kron->DelEdge(Node1, Node2, false);
			E+=2;
		}
	}

	// DEBUG

	if (Kron->GetEdges() - ModelEdges != 0) {
		TStr Msg = "Edges count of rewired graph != required edges count: ModelEdges = ";
		Msg += ModelEdges; Msg += ", KronEdges = "; Msg += Kron->GetEdges();
		//Error("AddEdges", Msg);
	}

	// END DEBUG
}



// rewire edges according to smoothed diaps
void Rewire(PNGraph& Kron, const vector<Diap>& SmoothedDiaps, const TIntPr& OutDegR, vector<int>& Prev, ofstream& TFile){
	TRnd Rnd;
	int ModelEdges = Kron->GetEdges();
	TFltPrV KronDeg;
	TExeTm execTime;
	TSnap::GetInDegCnt(Kron, KronDeg);
	TFile << "Time of getting degree sequence of Kronecker graph: " <<  execTime.GetSecs() << endl;
	execTime.Tick();
	//PrintDegDistr(KronDeg, "Kron.tab");
	KronDeg.Sort();
	const TInt& DegMin = OutDegR.Val1, &DegMax = OutDegR.Val2;
	vector<Diaps> DPlus, DMinus;
	GetDiaps(DPlus, DMinus, SmoothedDiaps, KronDeg, DegMin, DegMax, Prev);
	TFile << "Time of GetDiaps(): " <<  execTime.GetSecs() << endl;
	execTime.Tick();
	GetRewireStrategies(DPlus, DMinus);
	TFile << "Time of GetRewireStrategies(): " <<  execTime.GetSecs() << endl;
	execTime.Tick();
	//PrintDiapsInfo(DPlus, "DPlus.tab");
	//PrintDiapsInfo(DMinus, "DMinus.tab");
	Rewire(Kron, DPlus, DMinus, DegMin.Val, DegMax.Val);
	TFile << "Time of Rewire(): " <<  execTime.GetSecs() << endl;
	int Diff = Kron->GetEdges() - ModelEdges;
	cout << "Difference of edges: " << Kron->GetEdges() - ModelEdges << endl;
	AddEdges(Kron, Diff, DegMin, DegMax, ModelEdges, DMinus, DPlus);
	
	
}

// type 1: RD == 0
// type 2: RD == -1
// type 3: RD == 1
// type 4: RD > 1 
// type 5: 0 < RD < 1
int GetRelDiffType(double RD){
	if (RD == 0) return 1;
	if (RD == -1) return 2;
	if (RD == 1) return 3;
	if (RD > 1) return 4;
	return 5;
}

// get smoothed diapasons for scaling
void GetSmoothedDiaps(const TFltPrV& RelDiffNonCum, vector<Diap>& SmoothedDiaps, vector<int>& Prev){
	if (RelDiffNonCum.Len() == 0)
		Error("GetSmoothedDiaps", "Array size = 0");
	const int DegCount = RelDiffNonCum.Len(), 
		DegMin = RelDiffNonCum[0].Val1, 
		DegMax = RelDiffNonCum[DegCount-1].Val1, 
		DiffDegs = DegMax - DegMin + 1;
	TInt ToleranceVal = 4;
	TFltV DiapDevV;
	int PrevRDType = GetRelDiffType(RelDiffNonCum[0].Val2), RDType;
	int DiapBegin = 0, DiapEnd = 0; 
	int DiapIndex = 0, PrevDiapEnd = 0, PrevDeg = 0, Deg = 0;
	double Diff = 0;

	for (size_t i = 0; i < DegCount; i++){
			Deg = RelDiffNonCum[i].Val1;
			Diff = RelDiffNonCum[i].Val2;
			RDType = GetRelDiffType(Diff);
		
		// if there is no difference, diapason should be finalized
		//if (Diff == 0.0) Sign = DiapSign == true ? false : true;

		// if some degrees inside the diapasone are absent, add zero values of relative difference
		if (i != 0 && RDType == PrevRDType && Deg != PrevDeg + 1){
			printf("%d %d\n", PrevDeg, Deg);
			for (size_t i = 0; i < Deg-PrevDeg-1; i++)
				DiapDevV.Add(0);
			if (i == DegCount - 1){
				DiapEnd = i;
				DiapDevV.Add(Diff);
			}
		}

		// if sign was changed or it is last interval
		if (RDType != PrevRDType || i == DegCount - 1){
			int DegBegin = RelDiffNonCum[DiapBegin].Val1, 
				DegEnd = RelDiffNonCum[DiapEnd].Val1;
			//printf("DegBegin: %d DegEnd: %d\n", DegBegin, DegEnd);
			TInt DiapLength = DegEnd - DegBegin + 1;
			// if it is first interval or previous interval has enough length
			if (DiapIndex == 0 || DiapLength >= ToleranceVal){
				Diap NewDiap; 
				// [begin;end] as parts of [DegMin;DegMax]
				double ProbFirst = static_cast<double>(DegBegin - DegMin) / DiffDegs,
					ProbSecond = static_cast<double>(DegEnd - DegMin) / DiffDegs;
				int CalcDegBegin = static_cast<int>(ProbFirst * DiffDegs + DegMin + 0.5),
					CalcDegEnd = static_cast<int>(ProbSecond * DiffDegs + DegMin + 0.5);
				if (CalcDegBegin != DegBegin ||CalcDegEnd  != DegEnd){
					Error("GetSmoothedDiaps", "Diapason borders are violated");
				}
				NewDiap.first.Val1 = ProbFirst;
				NewDiap.first.Val2 = ProbSecond;

				printf("DiapBegin: %d DiapEnd: %d\n", DegBegin, DegEnd);

				if (DiapDevV.Len() != DiapLength){

					Error("GetSmoothedDiaps", "Inconsistent size of DiapDevV vector");
				}

				for (size_t i = 0; i < DiapDevV.Len(); i++)
					NewDiap.second.Add(DiapDevV[i]);
				SmoothedDiaps.push_back(NewDiap);

				// remember if previous interval is the neighbour
				if (PrevDiapEnd == DiapBegin - 1 && 
					RelDiffNonCum[PrevDiapEnd].Val1 + 1 == DegBegin)
					Prev.push_back(DiapIndex);

				PrevDiapEnd = DiapEnd;
				DiapBegin = i; DiapEnd = i; 
				// add value of RelDiffNonCum
				DiapDevV.Clr(); 
				DiapDevV.Add(Diff); 
				
				// for the last diapason
				if (i == DegCount - 1){
					NewDiap.first.Val1 = static_cast<double>(Deg - DegMin) / DiffDegs;
					NewDiap.first.Val2 = NewDiap.first.Val1;
					printf("DiapBegin: %d DiapEnd: %d\n", Deg, Deg);
					DiapDevV.Clr();
					DiapDevV.Add(Diff);
					NewDiap.second.Clr();
					NewDiap.second.Add(DiapDevV[0]);
					SmoothedDiaps.push_back(NewDiap);
				}

				PrevRDType = RDType;
				DiapIndex++;
			}
			// if DiapLength < ToleranceVal && Diap is finished
			else {
				DiapBegin = i;
				DiapEnd = i;
				DiapDevV.Clr(); 
				DiapDevV.Add(Diff); 
				PrevRDType = RDType;
			}
		}
		else {
			DiapEnd = i;
			DiapDevV.Add(Diff);
		}
		PrevDeg = Deg;
	}
}


// get relative differences of degrees
void GetRelativeDiff(const TFltPrV& MDeg, const TFltPrV& KronDeg, TFltPrV&  RelDiffV, bool NonCum){
	/*PrintDegDistr(MDeg, "ModelBasic.tab");
	PrintDegDistr(KronDeg, "KronBasic.tab");*/
	int MDegCount = MDeg.Len(), KronDegCount = KronDeg.Len();
	int MinDegModel = static_cast<int>(MDeg[0].Val1), MaxDegModel = static_cast<int>(MDeg[MDegCount-1].Val1),
		MinDegKron = static_cast<int>(KronDeg[0].Val1), MaxDegKron = static_cast<int>(KronDeg[KronDegCount-1].Val1);
	int MinDeg = MinDegModel < MinDegKron ? MinDegModel : MinDegKron,
		MaxDeg = MaxDegModel > MaxDegKron ? MaxDegModel : MaxDegKron;
	double CurrDeg = MinDeg, MInd = 0, KronInd = 0;
	if (NonCum){
		while (1){
			double MDegVal, KronDegVal, MDegCount, KronDegCount;
			if (MInd < MDeg.Len()){
				MDegVal = MDeg[MInd].Val1; MDegCount = MDeg[MInd].Val2;
			}
			else {
				MDegVal = MaxDeg + 1; MDegCount = 0;
			}
			if (KronInd < KronDeg.Len()) {
				KronDegVal = KronDeg[KronInd].Val1; KronDegCount = floor(KronDeg[KronInd].Val2 + 0.5);
			}
			else {
				KronDegVal = MaxDeg + 1; KronDegCount = 0;
			}
			bool MLessDeg = MDegVal < KronDegVal ? true : false;
			CurrDeg = MLessDeg ? MDegVal : KronDegVal;
			if (CurrDeg > MaxDeg) break;
			double RelDiff;
			if (MDegVal == CurrDeg && KronDegVal == CurrDeg){
				//RelDiff = (MDegCount - KronDegCount) / KronDegCount;
				if (MDegCount == KronDegCount)
					RelDiff = 0;
				else{ 
					// if averaged KronDegCount < 0.5
					if (KronDegCount == 0)
						RelDiff = 1;
					else 
						RelDiff = log10(MDegCount+1) / log10(KronDegCount+1);
				}
				
				//RelDiff = KronDegCount / MDegCount;
				MInd++; KronInd++;
			}
			else if (MLessDeg){
				RelDiff = 1;
				MInd++;
			}
			else {
				//RelDiff = KronDegCount;
				RelDiff = -1;
				KronInd++;
			}
			TFltPr RelDiffPr(CurrDeg, RelDiff);
			//printf("MDegCount=%3.2f KronDegCount=%3.2f RelDiff = %3.2f\n", MDegCount, KronDegCount, RelDiff);
			RelDiffV.Add(RelDiffPr);
		}
	}
}

double GetAvgDeviation(const TFltPrV& ModelDegCnt, const TFltPrV& KronDegCnt){
	double Dev = 0.0;
	bool SeqViewed = false;
	double MaxModelDeg = ModelDegCnt[ModelDegCnt.Len()-1].Val1, MaxKronDeg = KronDegCnt[KronDegCnt.Len()-1].Val1,
		MaxModelCount = ModelDegCnt[ModelDegCnt.Len()-1].Val2, MaxKronCount = KronDegCnt[KronDegCnt.Len()-1].Val2;
	if (MaxModelDeg == MaxKronDeg) 
		return (pow((MaxModelCount-MaxKronCount)/MaxModelCount,2));
	else {
		return (pow((MaxModelDeg-MaxKronDeg)/MaxModelDeg,2) + 1);
	}

	/*int IndexModel = 0, IndexKron = 0;
	while (1){
		if (IndexModel == ModelDegCnt.Len() && IndexKron == KronDegCnt.Len())
			break;
		double ModelDeg, KronDeg;
		if (IndexModel == ModelDegCnt.Len()) ModelDeg = 0;
		else ModelDeg = ModelDegCnt[IndexModel].Val1;
		if (IndexKron == KronDegCnt.Len()) KronDeg = 0;
		else KronDeg = KronDegCnt[IndexKron].Val1;

		if (ModelDeg < KronDeg ){
			Dev += 1 * (ModelDeg/MaxModelDeg); IndexModel++;
		}
		else if (ModelDeg == KronDeg){
			double ModelCount = ModelDegCnt[IndexModel].Val2, KronCount = KronDegCnt[IndexKron].Val2;
			Dev += pow((ModelCount-KronCount)/ModelCount, 2) * (ModelDeg/MaxModelDeg); IndexModel++; IndexKron++;
		}
		else {
			Dev += 1 * (ModelDeg/MaxModelDeg); IndexKron++;
		}
	}
	return Dev;*/
}