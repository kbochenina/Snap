#include "stdafx.h"

ofstream TFile;

void DefaultSettings(vector<TStr>& out){
	out.clear();
	printf("Applying default settings...\n");
	TStr s("1 ");
	for (int i = 0; i < NFUNC; i++){
		out.push_back(s);
	}
}

void ReadParameters(TStr settingsFN, vector<TStr>& out){
	ifstream f;
	f.open(settingsFN.CStr());
	if (!f.is_open())
	{
		while (1){
			printf("Error while opening file %s with command line settings. Apply default settings? (y/n)", settingsFN.CStr());
			char ch;
			scanf("%c",&ch);
			if (ch == 'Y' || ch == 'y'){
				DefaultSettings(out);
				break;
			}
			else if (ch == 'N' || ch == 'n'){
				printf("\nProgram terminated");
				exit(1);
			}
		}
	}
	string insteadOfName = "1 ";
	for (int i = 0; i < NFUNC; i++){
		string s;
		getline(f,s);
		if ( s!= FUNCNAMES[i]){
			printf("Wrong syntax in settings file. ");
			DefaultSettings(out);
			break;
		}
		bool isComment = true;
		while (isComment){
			getline(f,s);
			if (s.find_first_of("//")!=0)
				isComment = false;
		}
		
		if (s == "default")
			s = insteadOfName;
		else s = insteadOfName + s;
		TStr ts(s.c_str());
		out.push_back(ts);
	}
	f.close();
}

int GraphGen(const TStr args, PNGraph &GD){
	Env = TEnv(args, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Graph generation. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
	const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "output.txt", "Output graph filename");
	const TStr Plot = Env.GetIfArgPrefixStr("-g:", "e", "Which generator to use:"
	"\n\tf: Complete graph. Required parameters: n (number of nodes)"
	"\n\ts: Star graph. Required parameters: n (number of nodes)"
	"\n\t2: 2D Grid. Required parameters: n (number of rows), m (number of columns)"
	"\n\te: Erdos-Renyi (G_nm). Required parameters: n (number of nodes), m (number of edges)"
	"\n\tk: Random k-regular graph. Required parameters: n (number of nodes), k (degree of every node)"
	"\n\tb: Albert-Barabasi Preferential Attachment. Required parameters: n (number of nodes), k (edges created by each new node)"
	"\n\tp: Random Power-Law graph. Required parameters: n (number of nodes), p (power-law degree exponent)"
	"\n\tc: Copying model by Kleinberg et al. Required parameters: n (number of nodes), p (copying probability Beta)"
	"\n\tw: Small-world model. Required parameters: n (number of nodes), k (each node is connected to k nearest neighbors in ring topology), p (rewiring probability)\n"
	);
	const int N = Env.GetIfArgPrefixInt("-n:", 1000, "Number of nodes");
	const int M = Env.GetIfArgPrefixInt("-m:", 5000, "Number of edges");
	const double P = Env.GetIfArgPrefixFlt("-p:", 0.1, "Probability/Degree-exponent");
	const int K = Env.GetIfArgPrefixInt("-k:", 3, "Degree");
  
	if (Env.IsEndOfRun()) { return 0; }
	TInt::Rnd.PutSeed(0); // initialize random seed
	printf("Generating...\n");

	TStr DescStr;
	PUNGraph G;
	if (Plot == "f") {
	G = TSnap::GenFull<PUNGraph>(N);
	DescStr = TStr::Fmt("Undirected complete graph.");
	} else
	if (Plot == "s") {
	G = TSnap::GenStar<PUNGraph>(N, false);
	DescStr = TStr::Fmt("Undirected star graph (1 center node connected to all other nodes).");
	} else
	if (Plot == "2") {
	G = TSnap::GenGrid<PUNGraph>(N, M, false);
	DescStr = TStr::Fmt("Undirected 2D grid of %d rows and %d columns.", N, M);
	} else
	if (Plot == "e") {
	G = TSnap::GenRndGnm<PUNGraph>(N, M, false);
	DescStr = TStr::Fmt("Undirected Erdos-Renyi random graph.");
	} else
	if (Plot == "k") {
	G = TSnap::GenRndDegK(N, K);
	DescStr = TStr::Fmt("Undirected k-regular random graph (every node has degree K).");
	} else
	if (Plot == "b") {
	G = TSnap::GenPrefAttach(N, K);
	DescStr = TStr::Fmt("Undirected Albert-Barabasi Preferential Attachment graph (each new node creades k preferentially attached edges).");
	} else
	if (Plot == "p") {
	G = TSnap::GenRndPowerLaw(N, P, true);
	DescStr = TStr::Fmt("Random Graph with Power-Law degree distribution with exponent P.");
	} else
	if (Plot == "c") {
	G = TSnap::ConvertGraph<PUNGraph>(TSnap::GenCopyModel(N, P));
	DescStr = TStr::Fmt("Copying model by Kleinberg et al. Node u comes, selects a random v, and with prob P it links to v, with 1-P links u links to neighbor of v. Power-law degree slope is 1/(1-P).");
	} else
	if (Plot == "w") {
		G = TSnap::GenSmallWorld(N, K, P);
		DescStr = TStr::Fmt("Watts-Strogatz Small-world model. Every node links to K other nodes.");
	}
	printf("done.\n");
	// save graph
	//TSnap::SaveEdgeList(G, OutFNm, DescStr);
	GD = TSnap::ConvertGraph<PNGraph>(G);
	

	return 0;
	
	Catch
		printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
		return 0;
	}


int InitKronecker(const TStr args, PNGraph &GD, TKronMtx& FitMtx){
	Env = TEnv(args, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Kronecker graphs. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
	Env = TEnv(args, TNotify::StdNotify);
	//const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../as20graph.txt", "Input graph file (single directed edge per line)");
	TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "initMatrix.txt", "Output file prefix");
	const TInt NZero = Env.GetIfArgPrefixInt("-n0:", 2, "Innitiator matrix size");
	const TStr InitMtx = Env.GetIfArgPrefixStr("-m:", "0.9 0.7; 0.5 0.2", "Init Gradient Descent Matrix (R=random)").GetLc();
	const TStr Perm = Env.GetIfArgPrefixStr("-p:", "d", "Initial node permutation: d:Degree, r:Random, o:Order").GetLc();
	const TInt GradIter = Env.GetIfArgPrefixInt("-gi:", 10, "Gradient descent iterations");
	const TFlt LrnRate = Env.GetIfArgPrefixFlt("-l:", 1e-5, "Learning rate");
	const TFlt MnStep = Env.GetIfArgPrefixFlt("-mns:", 0.005, "Minimum gradient step");
	const TFlt MxStep = Env.GetIfArgPrefixFlt("-mxs:", 0.05, "Maximum gradient step");
	const TInt WarmUp =  Env.GetIfArgPrefixInt("-w:", 10000, "Samples to warm up");
	const TInt NSamples = Env.GetIfArgPrefixInt("-s:", 100000, "Samples per gradient estimation");
	//const TInt GradType = Env.GetIfArgPrefixInt("-gt:", 1, "1:Grad1, 2:Grad2");
	const bool ScaleInitMtx = Env.GetIfArgPrefixBool("-sim:", true, "Scale the initiator to match the number of edges");
	const TFlt PermSwapNodeProb = Env.GetIfArgPrefixFlt("-nsp:", 1.0, "Probability of using NodeSwap (vs. EdgeSwap) MCMC proposal distribution");
	//if (OutFNm.Empty()) { OutFNm = TStr::Fmt("%s-fit%d", InFNm.GetFMid().CStr(), NZero()); }
	printf("%s\n", OutFNm.CStr());
	// load graph
	/*PNGraph G;
	if (InFNm.GetFExt().GetLc()==".ungraph") {
	TFIn FIn(InFNm);  G=TSnap::ConvertGraph<PNGraph>(TUNGraph::Load(FIn), true); }
	else if (InFNm.GetFExt().GetLc()==".ngraph") {
	TFIn FIn(InFNm);  G=TNGraph::Load(FIn); }
	else {
	G = TSnap::LoadEdgeList<PNGraph>(InFNm, 0, 1);
	}*/
	// convert from undirected to directed via doubling of edges
	/*PNGraph GD;
	GD=TSnap::ConvertGraph<PNGraph>(G, true);*/
	// fit
	TKronMtx InitKronMtx = InitMtx=="r" ? TKronMtx::GetRndMtx(NZero, 0.1) : TKronMtx::GetMtx(InitMtx);
	InitKronMtx.Dump("INIT PARAM", true);
	TKroneckerLL KronLL(GD, InitKronMtx, PermSwapNodeProb);
	if (ScaleInitMtx) {
		InitKronMtx.SetForEdges(GD->GetNodes(), GD->GetEdges()); }
	KronLL.InitLL(GD, InitKronMtx);
	InitKronMtx.Dump("SCALED PARAM", true);
	KronLL.SetPerm(Perm.GetCh(0));
	double LogLike = 0;
	//if (GradType == 1) {
	LogLike = KronLL.GradDescent(GradIter, LrnRate, MnStep, MxStep, WarmUp, NSamples);
	//} else if (GradType == 2) {
	//  LogLike = KronLL.GradDescent2(GradIter, LrnRate, MnStep, MxStep, WarmUp, NSamples); }
	//else{ Fail; }
	//const TKronMtx& FitMtx = KronLL.GetProbMtx();
	FitMtx = KronLL.GetProbMtx();
	FILE *F = fopen(OutFNm.CStr(), "w");
//	fprintf(F, "Input\t%s\n", InFNm.CStr());
	TStrV ParamV; Env.GetCmLn().SplitOnAllCh(' ', ParamV);
	fprintf(F, "Command line options\n");
	for (int i = 0; i < ParamV.Len(); i++) {
		fprintf(F, "\t%s\n", ParamV[i].CStr()+(ParamV[i][0]=='-'?1:0)); }
	fprintf(F, "Loglikelihood\t%10.2f\n", LogLike);
	fprintf(F, "Absolute error (based on expected number of edges)\t%f\n", KronLL.GetAbsErr());
	fprintf(F, "RunTime\t%g\n", ExeTm.GetSecs());
	fprintf(F, "Estimated initiator\t%s\n", FitMtx.GetMtxStr().CStr());
	fclose(F);

	Catch
		printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
	return 0;
}

void RemoveZeroDegreeNodes(PNGraph& out){
	TRnd rnd;
	int nodesCount = out->GetNodes();
	for (int i = 0; i < nodesCount; i++){
		if (out->GetNI(i).GetInDeg() == 0){
			double val = rnd.GetUniDev();
			int nodeId = static_cast<int>(val * nodesCount);
			out->AddEdge(nodeId, i);
		}
	}
}

int KroneckerGen(const TInt NIter, const TKronMtx& FitMtx, PNGraph& out, const TStr& OutFNm){
	Env.PrepArgs(TStr::Fmt("Kronecker graphs. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
	const TKronMtx& SeedMtx = FitMtx;
	printf("\n*** Seed matrix:\n");
	SeedMtx.Dump();
	printf("\n*** Kronecker:\n");
	// slow but exact O(n^2) algorightm
	//PNGraph Graph = TKronMtx::GenKronecker(SeedMtx, NIter, true, Seed); 
	// fast O(e) approximate algorithm
	out = TKronMtx::GenFastKronecker(SeedMtx, NIter, true, 0);

	RemoveZeroDegreeNodes(out);

	// save edge list
	TSnap::SaveEdgeList(out, OutFNm, TStr::Fmt("Kronecker Graph: seed matrix [%s]", FitMtx.GetMtxStr().CStr()));
	Catch
		printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
	return 0;
}

void ReadPNGraphFromFile(const TStr args, PNGraph& G){
	Try
	Env = TEnv(args, TNotify::StdNotify);
	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../as20graph.txt", "Input graph file (single directed edge per line)");
	// load graph
	G = TSnap::LoadEdgeList<PNGraph>(InFNm, 0, 1);
	Catch
}



int FindVal1Elem(const TFltPrV& vec, const TFlt& elem, bool& isExact){
	for (int i = 0; i < vec.Len(); i++){
		if (vec[i].Val1.Val == elem){
			isExact = true;
			return i;
		}
		if (vec[i].Val1.Val > elem){
			return i-1;
		}
	}
}

void GetMtxFromSepLine(const TStr& line, const TStr& separator, TFltV& matrix){
	TStrV strVals;
	line.SplitOnAllAnyCh(separator, strVals);
	for (int i = 0; i < strVals.Len(); i++) 
		matrix.Add(strVals[i].GetFlt());
}

// get degrees from current and add it to degrees
void AddDegreeStat(const TFltPrV& current, TFltPrV& degrees, TIntPrV& samples){
	for (int j = 0; j < current.Len(); j++){
		const TFltPr& elem = current[j];
		const double& deg = elem.Val1.Val, &nodesCount = elem.Val2.Val;
		bool wasFound = false;
		// silly search
		for (int k = 0; k < degrees.Len(); k++){
			if (degrees[k].Val1.Val == deg){
				degrees[k].Val2.Val += nodesCount;
				samples[k].Val2.Val++;
				wasFound = true; break;
			}
		}
		if (!wasFound){
			TFlt d(deg), n(nodesCount);
			TFltPr val(d,n);
			degrees.Add(val);
			TInt di(static_cast<int>(deg));
			TIntPr valI(di, 1);
			samples.Add(valI);
		}
	}
}

// get graph and add in and out degrees to cumulative vectors
void AddDegreesStat(TFltPrV& deg, TIntPrV& samples, const PNGraph& G, bool isIn){
	TFltPrV current;
	if (isIn)
		TSnap::GetInDegCnt(G, current);
	else 
		TSnap::GetOutDegCnt(G, current);
	AddDegreeStat(current, deg, samples);
}

void GetAvgDegreeStat (TFltPrV& deg, const TIntPrV& samples){
	for (int i = 0; i < deg.Len(); i++)
		deg[i].Val2.Val /= static_cast<double>(samples[i].Val2.Val);
}

void GetAvgDegreeStat (TFltPrV& deg, const TInt& NKron){
	for (int i = 0; i < deg.Len(); i++)
		deg[i].Val2.Val /= NKron;
}

void GetPoints(const TFlt& maxDegLog, const TFlt& minDegLog, const int& NInt, const TFltPrV& base, TFltPrV& points){
	int beginIndex = 0;
	// ignore nodes with zero degree (for Kronecker graphs)
	/*if (base[0].Val1.Val != 0)
		points.Add(base[beginIndex]);
	else {
		points.Add(base[++beginIndex]);
	}*/
	points.Add(base[beginIndex]);
	TFlt baseMaxDeg = base[base.Len()-1].Val1.Val,
		baseMinDeg = base[beginIndex].Val1.Val;
	for (int i = beginIndex + 1; i < NInt; i++){
		// deg - degree to be found in base
		TFlt degRound (pow (10, minDegLog.Val + i * (maxDegLog.Val - minDegLog.Val) / NInt));
		TInt degInt(static_cast<int>(degRound.Val));
		TFlt deg(degInt);
		// if deg < baseMinDeg (for cases when baseMinDeg > minDeg)
		if (deg.Val <= baseMinDeg)
			continue;
		// if deg > baseMaxDeg, add last point and finish
		if (deg.Val >= baseMaxDeg){
			points.Add(base[base.Len()-1]);
			break;
		}
		// we have two cases: when we can find an exact value of deg, or when we have not this value
		bool isExact = false;
		int index = FindVal1Elem(base, deg, isExact);
		if (isExact){
			points.Add(base[index]);
		}
		else 
		{
			TFltPr x;
			x.Val1.Val = deg;
			x.Val2.Val = ( base[index].Val2.Val + base [index + 1].Val2.Val ) / 2;
			points.Add(x);
		}
	}
}

void PrintDegDistr(const TFltPrV& distr, const TStr& OutFNm){
	FILE *F = stdout;
	if (! OutFNm.Empty()) F = fopen(OutFNm.CStr(), "wt");
	fprintf(F, "\n");
	fprintf(F, "  Degree\tVal\n");
	for (int i = 0; i < distr.Len(); i++){
		fprintf(F, "%f\t%f\n", distr[i].Val1.Val, distr[i].Val2.Val);
	}
	if (! OutFNm.Empty()) { fclose(F); }
}

void PrintDegDistr(const TIntPrV& distr, const TStr& OutFNm){
	FILE *F = stdout;
	if (! OutFNm.Empty()) F = fopen(OutFNm.CStr(), "wt");
	fprintf(F, "\n");
	fprintf(F, "  Degree\tVal\n");
	for (int i = 0; i < distr.Len(); i++){
		fprintf(F, "%d\t%d\n", distr[i].Val1.Val, distr[i].Val2.Val);
	}
	if (! OutFNm.Empty()) { fclose(F); }
}

void GetNodesEdgesCountFromDegDistr(const TFltPrV& deg, int& nodes, int& edges){
	double nodesD = 0, edgesD = 0;
	for (int i = 0; i < deg.Len(); i++){
		nodesD += deg[i].Val2.Val;
		edgesD += deg[i].Val1.Val * deg[i].Val2.Val;
	}
	nodes = static_cast<int>(nodesD);
	edges = static_cast<int>(edgesD);
//	edges /= 2; as Deg = inDeg + outDeg
}

void GetNodesEdgesCountFromAccDegDistr(const TFltPrV& deg, int& nodes, int& edges){
	double nodesD = deg[0].Val2.Val, edgesD = 0;
	for (int i = 0; i < deg.Len(); i++){
		if (i == deg.Len()-1)
			edgesD += deg[i].Val1.Val * deg[i].Val2.Val;
		else edgesD += deg[i].Val1.Val * (deg[i].Val2.Val - deg[i+1].Val2.Val);
	}
	nodes = static_cast<int>(nodesD);
	edges = static_cast<int>(edgesD);
	//	edges /= 2; as Deg = inDeg + outDeg
}

void SaveDegree(const TFltPrV& deg, const TStr& n, bool isIn, bool isCum, bool calcCum = true){
	TFltPrV d(deg);
	d.Sort();
	int nodes, edges;
	GetNodesEdgesCountFromDegDistr(d, nodes, edges);
	TSnap::PlotDegDistr(d, nodes, edges, n, n, isCum, false, isIn, calcCum);
}

void SaveAndPlot(const PNGraph& G, const TStr& name, bool isCum){
    TFltPrV in, out;
	TSnap::GetInDegCnt(G, in);
	TSnap::GetOutDegCnt(G, out);
	int nodes = G->GetNodes(), edges = G->GetEdges();
	/*TSnap::PlotInDegDistr(G, name, name, isCum, false);
	TSnap::PlotOutDegDistr(G, name, name, isCum, false);*/
	TSnap::PlotDegDistr(in, nodes, edges, name, name, isCum, false, true);
	TSnap::PlotDegDistr(out, nodes, edges, name, name, isCum, false, false);
}

void GetMinMaxLogDegree(const vector<TFltPrV>& distr, TFlt& minLog, TFlt& maxLog){
	TFltV vecMin, vecMax;
	for (size_t i = 0; i < distr.size(); i++){
		const TFltPrV &vec = distr[i]; 
		vecMin.Add(vec[0].Val1);
		vecMax.Add(vec[vec.Len()-1].Val1);
	}
	vecMin.Sort();
	vecMax.Sort();
	minLog = log10(vecMin[0]);
	maxLog= log10(vecMax[vecMax.Len()-1]);
}

void ExpBinning(const TFltPrV& deg, TFltPrV& degSparse, const int& BinRadix){
	PrintDegDistr(deg, "degtest");
	TFlt maxDeg(deg[deg.Len()-1].Val1.Val);
	bool maxPowerReached = false;
	// idx - index of border, previdx - index of previous border
	int power = 0, previdx = 0, idx = 0, binSize = 0;
	bool isExact;
	while (!maxPowerReached){
		double binBorder = pow(static_cast<double>(BinRadix), power++);
		if (power == 1){
			// if there are nodes with degree 1
			idx = FindVal1Elem(deg, 1, isExact);
			if (isExact){
				TFltPr val(1, deg[idx].Val2.Val);
				degSparse.Add(val);
				previdx = idx;
			}
		}
		else {
			if (binBorder >= maxDeg){
				// when last element of deg was previous bin border
				if (previdx == deg.Len() - 1)
					break;
				// if we have another elements
				binBorder = maxDeg;
				maxPowerReached = true;
			}
			// find next element
			idx = FindVal1Elem(deg, binBorder, isExact);
			// if bin size == 0
			if (previdx + 1 == idx && !isExact)
				continue;
			if (!isExact)
				idx = idx - 1;
			double sum = 0.0;
			binSize = idx - previdx;
			for (int i = previdx + 1; i <= idx; i++){
				sum += deg[i].Val2.Val;
			}
			sum /= binSize;
			double avgDeg = (binBorder + static_cast<double>(binBorder) / BinRadix) / 2.0;
			TFltPr val(avgDeg, sum);
			degSparse.Add(val);
			previdx = idx;
		}
	}
}

void GetCumDistr(const TFltPrV& nonCum, TFltPrV& res){
	for (int i = nonCum.Len() - 1; i >=0; i--){
		TFlt count;
		if (i == nonCum.Len() - 1)
			count = nonCum[i].Val2.Val;
		else
			count = nonCum[i].Val2.Val + res[res.Len()-1].Val2.Val;
		TFltPr val(nonCum[i].Val1, count);
		res.Add(val);
	}
	res.Sort();
}

void SaveSparse(const TFltPrV& G, const int& BinRadix, bool isIn, const TStr&name, bool isCum){
	TFltPrV deg(G), degSparse;
	if (isCum){
		deg.Clr();
		GetCumDistr(G, deg);
	}
	//GetPoints(maxLog, minLog, NInt, deg, degSparse);
	ExpBinning(deg, degSparse, BinRadix);

	/*if (isCum){
	PrintDegDistr(G, "G.Tab");
	PrintDegDistr(deg, "degTest.Tab");
	PrintDegDistr(degSparse, "degSparseTest.Tab");
	}*/
	//printf("%s: Nodes %d, edges %d\n", name.CStr(), nodes, edges);
	SaveDegree(degSparse, name, isIn, isCum, false);
}

void GetPowerLawDistrib(TIntV& DegSeqV, const int& NodesCount, const double& Gamma){
	// part of vertices with degree deg is proportional to pow (c * deg, -Gamma)
	int NodesAssigned = 0;
	// c - normalizing coefficient
	double degmax = pow (NodesCount, 1.00 / Gamma), c = 0.0;
	int degsum = 0;
	int degmaxInt = static_cast<int>(degmax);
	double sub = degmax - degmaxInt;
	if (sub != 0) --degmax;
	c = static_cast<double>(NodesCount) / pow(degmax, Gamma);
	double Z = 1.0;
	for (int i = 2; i <= degmax; i++)
		Z += 1.00 / pow (i, Gamma);
	for (int i = degmax; i >= 1; i--){
		double n = ( 1.00 / pow (i, Gamma) ) * c * NodesCount / Z;
		int nodes = static_cast<int>(n + 0.5);
		
		NodesAssigned += nodes;
		if (NodesAssigned > NodesCount)
		{
			int s = NodesAssigned - NodesCount;
			nodes -= s;
			NodesAssigned -= s;
		}
		for (int j = 0; j < nodes; j++){
			DegSeqV.Add(i);
			degsum += i;
		}
	}

	if (degsum % 2 != 0){
		// if we have non-even sum of degrees, we add 1 edge to one 1-degree node
		int idx = DegSeqV.SearchForw(1);
		DegSeqV[idx]++;
	}

	if (NodesAssigned != NodesCount){
		printf("Nodes assigned: %d, nodes count: %d", NodesAssigned, NodesCount);
	}
}

void GenRandomMtx(const int& MtxRndSize, TKronMtx& FitMtx){
	FitMtx.SetRndMtx(MtxRndSize);
}

void GenNewMtx(PNGraph& model, const TStr& args, TKronMtx& FitMtx){
	TExeTm execTime;
	InitKronecker(args, model, FitMtx);
	TFile << "Time of creation of init matrix: " <<  execTime.GetTmStr() << endl;
}

void ReadMtx(const TStr& Mtx, TKronMtx& FitMtx){
	TFltV matrix;
	GetMtxFromSepLine(Mtx, ";", matrix);
	FitMtx.GenMtx(matrix.Len() / 2);
	FitMtx.SetMtx(matrix);
}

void GetCumulativeDistrib(const PNGraph& model, TFltPrV& distr, bool isIn){
	if (isIn)
		TSnap::GetInDegCnt(model, distr);
	else 
		TSnap::GetOutDegCnt(model, distr);
	for (int i = distr.Len() - 2; i >= 0; i-- )
		distr[i].Val2.Val += distr[i+1].Val2.Val;
}

void ConvertToNonCum(TFltPrV& distr){
	for (int i = 0; i < distr.Len() - 1; i++ )
		distr[i].Val2.Val -= distr[i+1].Val2.Val;
}



// get model graph according to args
void GetModel(const TStr& args, PNGraph& G, const TStr& name, const TStr& Plt){
	Env = TEnv(args, TNotify::StdNotify);
	const TStr Gen = Env.GetIfArgPrefixStr("-p:", "gen", "How to get model graph: read (read from file, -i: file name); gen (use generator); deg (create with power-law degree distribution)");
	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "", "Input graph file (single directed edge per line)");
	const TInt NodesCount = Env.GetIfArgPrefixInt("-n:", 1024, "Nodes count");
	const TFlt Gamma = Env.GetIfArgPrefixFlt("-gamma:", 2.0, "Gamma");
	TExeTm execTime;
	if (Gen == "gen")
		GraphGen(args, G);
	else if (Gen == "read")
		ReadPNGraphFromFile(InFNm, G);
	else if (Gen == "deg"){
		TIntV DegSeqV;
		GetPowerLawDistrib(DegSeqV, NodesCount, Gamma);
		PUNGraph GU = TSnap::GenDegSeq(DegSeqV);
		G = TSnap::ConvertGraph<PNGraph>(GU);
	}
	if (Plt == "cum" || Plt == "all")
		SaveAndPlot(G, name.CStr(), true);
	if (Plt == "noncum" || Plt == "all")
		SaveAndPlot(G, name.CStr(), false);
	TFile << "Time of getting model: " <<  execTime.GetTmStr() << endl;
}

// read or get random mtx
bool GetMtx(const TStr& MtxArgs, TKronMtx& FitMtxModel){
	Env = TEnv(MtxArgs, TNotify::StdNotify);
	// how to generate initiator matrix
	const TStr Mtx = Env.GetIfArgPrefixStr("-m:", "random", "Init Kronecker matrix");
	// if matrix will be generated, its size is an argument of KRONFIT cmd line
	const TInt MtxRndSize = Env.GetIfArgPrefixInt("-rs:", 2, "Size of randomized Kronecker matrix");
	// get Kronecker init matrix
	if (Mtx == "create") return false;
	if (Mtx == "random")
		GenRandomMtx(MtxRndSize, FitMtxModel);
	else 
		ReadMtx(Mtx, FitMtxModel);
	return true;
}

void GenKron(const TStr& args, const TKronMtx& FitMtx, TFltPrV& inDegAvgKronM, TFltPrV& outDegAvgKronM){
	Env = TEnv(args, TNotify::StdNotify);
	TExeTm execTime;
	// number of Kronecker graphs to generate
	const TInt NKron = Env.GetIfArgPrefixInt("-n:", 1, "Number of generated Kronecker graphs");
	// if IsSampled == true, during averaging each value is divided on number of samples having this value,
	// otherwise value is divided on NKron
	const TStr IsSampled = Env.GetIfArgPrefixStr("-s:", "true", "Averaging by number of samples");
	// iterations of Kronecker product
	const TInt NIter = Env.GetIfArgPrefixInt("-i:", 10, "Iterations of Kronecker product");
	// output file name
	const TStr OutFnm = Env.GetIfArgPrefixStr("-o:", "kronGen.txt", "Output file name (default: krongen.txt)");
	// Kronecker model of graph
	PNGraph kron;
	TIntPrV samplesIn, samplesOut;
	double sec = 0.0;
	for (int i = 0; i < NKron; i++){
		execTime.Tick();
		KroneckerGen(NIter, FitMtx, kron, OutFnm);
		sec += execTime.GetSecInt();
		printf("Nodes count: %d, nodes with non-zero degree: %d\n", kron->GetNodes(), TSnap::CntNonZNodes(kron));
		AddDegreesStat(inDegAvgKronM, samplesIn, kron, true);
		AddDegreesStat(outDegAvgKronM, samplesOut, kron, false);
	}
	sec /= NKron;
	if (IsSampled == "true"){
		GetAvgDegreeStat(inDegAvgKronM, samplesIn);
		GetAvgDegreeStat(outDegAvgKronM, samplesOut);
	}
	else {
		GetAvgDegreeStat(inDegAvgKronM, NKron);
		GetAvgDegreeStat(outDegAvgKronM, NKron);
	}
	inDegAvgKronM.Sort();
	outDegAvgKronM.Sort();
	TFile << "Average time of generation of Kronecker product: " <<  execTime.GetTmStrFromSec(sec) << endl;
}


// plot all points without binning
void PlotPoints(const TFltPrV& inFirst, const TFltPrV& outFirst, const TFltPrV& inSecond, const TFltPrV& outSecond, const TStr& Plt){
		if (Plt == "cum" || Plt == "all"){
			SaveDegree(inFirst, "kronModel", true, true);
			SaveDegree(outFirst, "kronModel", false, true);
			SaveDegree(inSecond, "kronSmall", true, true);
			SaveDegree(outSecond, "kronSmall", false, true);
		}
		if (Plt == "noncum" || Plt == "all"){
			SaveDegree(inFirst, "kronModel", true, false);
			SaveDegree(outFirst, "kronModel", false, false);
			SaveDegree(inSecond, "kronSmall", true, false);
			SaveDegree(outSecond, "kronSmall", false, false);
		}
}


void PlotSparse(const vector<TFltPrV>& distr, const TStrV& names, bool isIn, const TStr& Plt, const TInt& BinRadix){
	TFlt minLog, maxLog;
	for (size_t i = 0; i < distr.size(); i++){
		if (Plt == "cum" || Plt == "all"){
			SaveSparse(distr[i], BinRadix, isIn, names[i], true);
		}
		if (Plt == "noncum" || Plt == "all")
			SaveSparse(distr[i], BinRadix, isIn, names[i], false);
	}
}

ofstream OpenFile(const TStr& fileName)
{
	Try
	ofstream f(fileName.CStr());
	if (f.is_open())
	return f;		
	IAssert(1);
	Catch
}

void GetGraphs(vector <TStr>& parameters)
{
	GetModel(parameters[GRAPHGEN], G, parameters[NAME], parameters[PLT]);
	// generate Kronecker initiator matrix using big graph
	TKronMtx FitMtxM;
	if (!GetMtx(parameters[MTXGEN], FitMtxM))
		GenNewMtx(G, parameters[KRONFIT], FitMtxM);

	// in and out average degree distribution for kronM (non-accumulated)
	TFltPrV inDegAvgKronM, outDegAvgKronM;
	GenKron(parameters[KRONGEN], FitMtxM, inDegAvgKronM, outDegAvgKronM);
}

// generates Kronecker model using configuration model of small model network
// and compare it to big network
void KroneckerByConf(vector<TStr> commandLineArgs){
	Try
	Env = TEnv(commandLineArgs[KRONTEST], TNotify::StdNotify);
	// type of plots
	const TStr Plt = Env.GetIfArgPrefixStr("-plt:", "noncum", "Type of plots (cum, noncum, all)");
	// full - all points of distrib will be plotted; expbin - exponential binning
	const TStr PType = Env.GetIfArgPrefixStr("-ptype:", "all", "How to plot (full, expbin, all)");
	// if there is a need for a plot of the small model
	const TStr PlotMS = Env.GetIfArgPrefixStr("-ms:", "false", "Plot of small model is required: true, false");
	// radix of binning
	const TInt BinRadix = Env.GetIfArgPrefixInt("-bin:", 2, "Radix for exponential binning");
	// time estimates file name
	const TStr TimeFile = Env.GetIfArgPrefixStr("-ot:", "time.dat", "Name of output file with time estimates");
	// generation of big model and its Kronecker product is required
	const TStr ModelGen = Env.GetIfArgPrefixStr("-mg:", "true", "Generation of big model and its Kronecker product is required");

	TFile = OpenFile(TimeFile.CStr());

	if (ModelGen == "true")
	{
		TFile << "Model graph" << endl;
		vector <TStr> parameters;
		copy(commandLineArgs.begin() + 1, commandLineArgs.begin() + NPARCOPY, parameters);
		parameters.push_back(PType); parameters.push_back(Plt); parameters.push_back("model");
		GetGraphs(parameters);
	}

	PNGraph model;
	
	GetModel(commandLineArgs[GRAPHGEN_M], model, "model", Plt);
	// generate Kronecker initiator matrix using big graph
	TKronMtx FitMtxM;
	if (!GetMtx(commandLineArgs[MTXGEN_M], FitMtxM))
		GenNewMtx(model, commandLineArgs[KRONFIT_M], FitMtxM);
	
	
	// in and out average degree distribution for kronM (non-accumulated)
	TFltPrV inDegAvgKronM, outDegAvgKronM;
	GenKron(commandLineArgs[KRONGEN_M], FitMtxM, inDegAvgKronM, outDegAvgKronM);
	
	// generate small graph with the same degree distribution
	PNGraph modelS;
	GetModel(commandLineArgs[GRAPHGEN_MS], modelS, "modelSmall", Plt);
	// use configuration model
	/*PUNGraph modelU = TSnap::GenConfModel(TSnap::ConvertGraph<PUNGraph>(modelSmall));
	modelSmall = TSnap::ConvertGraph<PNGraph>(modelU);*/
		
	// generate Kronecker initiator matrix using modelSmall
	TKronMtx FitMtxMS;
	if (!GetMtx(commandLineArgs[MTXGEN_MS], FitMtxMS))
		GenNewMtx(modelS, commandLineArgs[KRONFIT_MS], FitMtxMS);

	// create NKron graphs
	TFltPrV inDegAvgKronMS, outDegAvgKronMS;
	GenKron(commandLineArgs[KRONGEN_MS], FitMtxMS, inDegAvgKronMS, outDegAvgKronMS);
	
	if ( PType == "full" || PType == "all" )
		PlotPoints(inDegAvgKronM, outDegAvgKronM, inDegAvgKronMS, outDegAvgKronMS, Plt);

	if ( PType == "exp" || PType == "all"){
		TFltPrV mDegIn, mDegOut, msDegIn, msDegOut;
		TSnap::GetInDegCnt(model, mDegIn);
		TSnap::GetOutDegCnt(model, mDegOut);
		TSnap::GetInDegCnt(modelS, msDegIn);
		TSnap::GetOutDegCnt(modelS, msDegOut);
		vector<TFltPrV> distrIn, distrOut;
		TStrV names;
		distrIn.push_back(mDegIn); distrOut.push_back(mDegOut); names.Add("modelSparse");
		if (PlotMS == "true"){
			distrIn.push_back(msDegIn); distrOut.push_back(msDegOut); names.Add("modelSmallSparse");
		}
		distrIn.push_back(inDegAvgKronM); distrOut.push_back(outDegAvgKronM); names.Add("kronModelSparse");
		distrIn.push_back(inDegAvgKronMS); distrOut.push_back(outDegAvgKronMS); names.Add("kronSmallSparse");
		PlotSparse(distrIn, names, true, Plt, BinRadix);
		PlotSparse(distrOut, names, false, Plt, BinRadix);
	}

	TFile.close();

	Catch
}

