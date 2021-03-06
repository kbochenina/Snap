#include "stdafx.h"

ofstream TFile;


int BasicGraphGen(const TStr args, PNGraph &GD){
	Env = TEnv(args, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Graph generation. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
	const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "output.txt", "Output graph filename");
	const TStr Plot = Env.GetIfArgPrefixStr("-gt:", "e", "Which generator to use:"
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
	cout << "n0 = " << NZero << endl;
	// fit
	FILE *F = fopen(OutFNm.CStr(), "w");
	TKronMtx InitKronMtx = InitMtx=="r" ? TKronMtx::GetRndMtx(NZero, 0.1) : TKronMtx::GetMtx(InitMtx);
	InitKronMtx.Dump("INIT PARAM", true);
	TKroneckerLL KronLL(GD, InitKronMtx, PermSwapNodeProb);
	fprintf(F, "INIT PARAM\t%s, MTX SUM %f\n", InitKronMtx.GetMtxStr().CStr(), InitKronMtx.GetMtxSum());
	if (ScaleInitMtx) {
		InitKronMtx.SetForEdges(GD->GetNodes(), GD->GetEdges()); }
	KronLL.InitLL(GD, InitKronMtx);
	InitKronMtx.Dump("SCALED PARAM", true);
	fprintf(F, "SCALED PARAM\t%s, MTX SUM %f\n", InitKronMtx.GetMtxStr().CStr(), InitKronMtx.GetMtxSum());
	KronLL.SetPerm(Perm.GetCh(0));
	double LogLike = 0;
	//if (GradType == 1) {
	LogLike = KronLL.GradDescent(GradIter, LrnRate, MnStep, MxStep, WarmUp, NSamples);
	//} else if (GradType == 2) {
	//  LogLike = KronLL.GradDescent2(GradIter, LrnRate, MnStep, MxStep, WarmUp, NSamples); }
	//else{ Fail; }
	//const TKronMtx& FitMtx = KronLL.GetProbMtx();
	FitMtx = KronLL.GetProbMtx();

//	fprintf(F, "Input\t%s\n", InFNm.CStr());
	TStrV ParamV; Env.GetCmLn().SplitOnAllCh(' ', ParamV);
	fprintf(F, "Command line options\n");
	for (int i = 0; i < ParamV.Len(); i++) {
		fprintf(F, "\t%s\n", ParamV[i].CStr()+(ParamV[i][0]=='-'?1:0)); }
	fprintf(F, "Loglikelihood\t%10.2f\n", LogLike);
	fprintf(F, "Absolute error (based on expected number of edges)\t%f\n", KronLL.GetAbsErr());
	fprintf(F, "RunTime\t%g\n", ExeTm.GetSecs());
	fprintf(F, "Estimated initiator\t%s, mtx sum %f\n", FitMtx.GetMtxStr().CStr(), FitMtx.GetMtxSum());
	/*if (ScaleInitMtx) {
		FitMtx.SetForEdgesNoCut(GD->GetNodes(), GD->GetEdges()); }
	fprintf(F, "Scaled initiator\t%s, mtx sum %f\n", FitMtx.GetMtxStr().CStr(), FitMtx.GetMtxSum());
	FitMtx.Normalize();
	fprintf(F, "Normalized initiator\t%s, mtx sum %f\n", FitMtx.GetMtxStr().CStr(), FitMtx.GetMtxSum());*/
	fclose(F);

	Catch
		printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
	return 0;
}



void KroneckerGen(PNGraph& out, const TKronMtx& FitMtx, const TInt NIter, const TStr& IsDir, const TIntPr& InDegR, const TIntPr& OutDegR, double NoiseCoeff){
	//Env.PrepArgs(TStr::Fmt("Kronecker graphs. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
	const TKronMtx& SeedMtx = FitMtx;
	/*printf("\n*** Seed matrix:\n");
	SeedMtx.Dump();*/
	//printf("\n*** Kronecker:\n");
	bool Dir = IsDir == "true" ? true: false;
	if (Dir){
		//printf("Directed graph with restrictions. Required functional is under construction. Undirected graph will be generated instead\n");
		Dir = false;
	}
	
	// if we have constraints on degrees, run corresponding algorithm
	if (InDegR.Val1 != 0 || OutDegR.Val1 != 0 || InDegR.Val2 != INT_MAX || InDegR.Val2 != INT_MAX) {
		if (Dir){
			//printf("Directed graph with restrictions. Required functional is under construction. Undirected graph will be generated instead\n");
			Dir = false;
		}
		// we want to generate graph without restrictions on degrees
		TIntPr DegR(0, INT_MAX);
		TKronMtx::GenFastKronecker(out, SeedMtx, NIter, Dir, DegR, DegR, NoiseCoeff);
	}
	else {
		TKronMtx::GenFastKronecker(out, SeedMtx, NIter, static_cast<int>(pow(SeedMtx.GetMtxSum(), NIter)), Dir);
		// slow and exact version
		//TKronMtx::GenKronecker(out, SeedMtx, NIter, Dir);
	}
	//cout << "Kronecker edges: " << out->GetEdges() << endl;
	//PrintNodeDegrees(out, SeedMtx, NIter);
	TKronMtx::RemoveZeroDegreeNodes(out, SeedMtx, NIter, InDegR.Val1, InDegR.Val2);
	//printf("             %d edges [%s]\n",out->GetEdges(), ExeTm.GetTmStr());
	// save edge list
	//TSnap::SaveEdgeList(out, OutFNm, TStr::Fmt("Kronecker Graph: seed matrix [%s]", FitMtx.GetMtxStr().CStr()));
	Catch
		//printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
}


void GetPowerLawDistrib(TIntV& DegSeqV, const int& NodesCount, const double& Gamma){
	TExeTm execTime;
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
	TFile << "Time of getting power law distrib: " <<  execTime.GetTmStr() << endl;
}

void GenRandomMtx(const int& MtxRndSize, TKronMtx& FitMtx){
	FitMtx.SetRndMtx(MtxRndSize);
}

void GenNewMtx(PNGraph& model, const TStr& args, TKronMtx& FitMtx){
	TExeTm execTime;
	InitKronecker(args, model, FitMtx);
	TFile << "Time of creation of init matrix: " <<  execTime.GetTmStr() << endl;
}



// get model graph according to args
void GetModel(const TStr& Args, PNGraph& G){
	Env = TEnv(Args, TNotify::NullNotify);
	const TStr Gen = Env.GetIfArgPrefixStr("-g:", "gen", "How to get model graph: read, gen, deg, genpy");
	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "", "Input graph file (single directed edge per line)");
	
	TExeTm execTime;
	if (Gen == "gen")
		BasicGraphGen(Args, G);
	else if (Gen == "read")
		ReadPNGraphFromFile(InFNm, G);
	else if (Gen == "deg"){
		const TInt NodesCount = Env.GetIfArgPrefixInt("-n:", 1024, "Nodes count");
		const TFlt Gamma = Env.GetIfArgPrefixFlt("-gamma:", 2.0, "Gamma");
		TIntV DegSeqV;
		GetPowerLawDistrib(DegSeqV, NodesCount, Gamma);
		TExeTm T;
		PUNGraph GU = TSnap::GenDegSeq(DegSeqV);
		TFile << "Time of execution of configuration model: " <<  T.GetTmStr() << endl;
		G = TSnap::ConvertGraph<PNGraph>(GU);
	}
	if (Gen == "genpy")
	{
		PUNGraph GU;
		GenPy(GU, TFile, Args);	
		G = TSnap::ConvertGraph<PNGraph>(GU);
	}
	TFile << "Time of getting model: " <<  execTime.GetTmStr() << endl;
	//TFile << "Model graph: " << G->GetNodes() << " nodes, " << G->GetEdges() << " edges\n";
	/*TIntV DegV;
	TSnap::GetDegSeqV(G, DegV);
	execTime.Tick();
	PUNGraph Conf = TSnap::GenConfModel(DegV);
	TFile << "Time of getting configuration model: " <<  execTime.GetTmStr() << endl;
	cout << "Undirected configuration model: " << Conf->GetNodes() << " nodes, " << Conf->GetEdges() << " edges\n";*/
	//PNGraph ConfD = TSnap::ConvertGraph<PNGraph>(Conf);
	//SaveAndPlot(ConfD, "conf", false);
	//TFile << "Clustering coefficient of configuration model: " << TSnap::GetClustCf(ConfD) << endl;
	//TSnap::PlotClustCf(ConfD,"conf");
}

// read or get random mtx
bool GetMtx(const TStr& MtxArgs, TKronMtx& FitMtxModel){
	Env = TEnv(MtxArgs, TNotify::StdNotify);
	// how to generate initiator matrix
	const TStr Mtx = Env.GetIfArgPrefixStr("-m:", "random", "Init Kronecker matrix");
	// if matrix will be generated, its size is an argument of KRONFIT cmd line
	const TInt MtxSize = Env.GetIfArgPrefixInt("-rs:", 2, "Size of randomized Kronecker matrix");
	// get Kronecker init matrix
	if (Mtx == "create") return false;
	if (Mtx == "random")
		GenRandomMtx(MtxSize, FitMtxModel);
	else 
		ReadMtx(Mtx, MtxSize, FitMtxModel);
	return true;
}



void GenKron(const TStr& Args, TKronMtx& FitMtx, TFltPrV& KronDegAvgIn, TFltPrV& KronDegAvgOut, vector<BaseDiap>& BaseDiaps, PNGraph& Kron){
	Env = TEnv(Args, TNotify::NullNotify);
	TExeTm ExecTime;
	// number of Kronecker graphs to generate
	const TInt NKron = Env.GetIfArgPrefixInt("-n:", 1, "Number of generated Kronecker graphs");
	// if IsSampled == true, during averaging each value is divided on number of samples having this value, otherwise value is divided on NKron
	const TStr IsSampled = Env.GetIfArgPrefixStr("-s:", "true", "Averaging by number of samples");
	// iterations of Kronecker product
	const TInt NIter = Env.GetIfArgPrefixInt("-i:", 10, "Iterations of Kronecker product");
	// output file name
	const TStr OutFnm = Env.GetIfArgPrefixStr("-o:", "kronGen.txt", "Output file name (default: krongen.txt)");
	// output file name
	const TStr ScaleMtx = Env.GetIfArgPrefixStr("-scalemtx:", "true", "Scale init matrix to match number of edges");
	// is graph directed?
	TStr IsDir = Env.GetIfArgPrefixStr("-isdir:", "false", "Produce directed graph (true, false)");
	// restrictions to in- and out- degrees count for 1 vertex
	const TInt InMin = Env.GetIfArgPrefixInt("-inmin:", numeric_limits<int>::lowest(), "In-degree minimum");
	const TInt InMax = Env.GetIfArgPrefixInt("-inmax:", INT_MAX, "In-degree maximum");
	const TInt OutMin = Env.GetIfArgPrefixInt("-outmin:", numeric_limits<int>::lowest(), "Out-degree minimum");
	const TInt OutMax = Env.GetIfArgPrefixInt("-outmax:", INT_MAX, "Out-degree maximum");
	// part of maximum possible noise to use in initiator matrix
	const double NoiseCoeff = Env.GetIfArgPrefixFlt("-noise:", 0, "NoiseCoeff");
	// scaling factor for maximum expected degree (REMOVE)
	const double MinScale = Env.GetIfArgPrefixFlt("-minscale:", 0.1, "MinScale");
	const TInt ModelNodes = Env.GetIfArgPrefixInt("-nodes:", 0, "Model nodes count");
	const TInt ModelEdges = Env.GetIfArgPrefixInt("-e:", 0, "Model edges count");
	const TInt MaxModelInDeg = Env.GetIfArgPrefixInt("-modelinmax:", INT_MAX, "Max model in-degree");
	const TInt MaxModelOutDeg = Env.GetIfArgPrefixInt("-modeloutmax:", INT_MAX, "Max model out-degree");
	const TInt MinModelInDeg = Env.GetIfArgPrefixInt("-modelinmin:", 0, "Min model in-degree");
	const TInt MinModelOutDeg = Env.GetIfArgPrefixInt("-modeloutmin:", 0, "Min model out-degree");

	TFlt ExpectedNodes = FitMtx.GetNodes(NIter), ExpectedEdges = FitMtx.GetEdges(NIter);
	TFile << "Model nodes: " << ModelNodes << ", expected model edges: " << ModelEdges << endl;
	TFile << "Kronecker nodes: " << ExpectedNodes << ", expected Kronecker edges: " << ExpectedEdges << endl;
	TFile << "Maximum degree in model graph: in " << MaxModelInDeg << ",  out " << MaxModelOutDeg << endl;

	
	//const TIntPr InDegR(0, 1000000); const TIntPr OutDegR(0, 1000000);
	const TIntPr InDegR(MinModelInDeg, MaxModelInDeg); const TIntPr OutDegR(MinModelOutDeg, MaxModelOutDeg);
	if (InDegR.Val1 * ExpectedNodes > ExpectedEdges || OutDegR.Val1 * ExpectedNodes > ExpectedEdges || InDegR.Val2 < InDegR.Val1 || OutDegR.Val2 < OutDegR.Val1) 
		Error("GenKron", "Constraints do not match to the number of edges");
	if (IsDir == "false" && (InDegR.Val1 != OutDegR.Val1 || InDegR.Val2 != OutDegR.Val2))
		Error("GenKron", "InDegR and OutDegR should be the same for undirected graph");

	// Kronecker  graph (UNCOMMENT FOR REAL USE)
	//PNGraph Kron;
	TIntPrV SamplesIn, SamplesOut;
	double Sec = 0.0;
	int AvgMaxOutDeg = 0, AvgMaxInDeg = 0, MinMaxOutDeg = 0, MaxMaxOutDeg = 0, MinMaxInDeg = 0, MaxMaxInDeg = 0;

	for (int i = 0; i < NKron; i++){
		ExecTime.Tick();
		KroneckerGen(Kron, FitMtx, NIter, IsDir, InDegR, OutDegR, NoiseCoeff);
		Rewire(Kron, BaseDiaps, OutDegR, TFile);
		//TKronMtx::RemoveZeroDegreeNodes(Kron, FitMtx, NIter, InDegR.Val1, InDegR.Val2);
		/*if (IsDir == "false" && !CheckReciprocity(Kron)){
			Error("GenKron", "Violation of reciprocity for undirected graph");
		}*/
		Sec += ExecTime.GetSecs();
		printf("Calculating maximum degree...\n");
		int MaxOutDeg = GetMaxMinDeg(Kron, IsDir, "false", "true"), MaxInDeg = GetMaxMinDeg(Kron, IsDir, "true", "true");
		CompareDeg(i, MaxOutDeg, MinMaxOutDeg, MaxMaxOutDeg, AvgMaxOutDeg);
		CompareDeg(i, MaxInDeg, MinMaxInDeg, MaxMaxInDeg, AvgMaxInDeg);

		//printf("Nodes count: %d, nodes with non-zero degree %d, edges count %d\n max deg = %d\n", kron->GetNodes(), TSnap::CntNonZNodes(kron), kron->GetEdges(), MaxDeg);
		if (i == NKron - 1){
			//TFile << "Clustering coefficient: " << TSnap::GetClustCf(kron) << endl;
			//TSnap::PlotClustCf(kron,"kronSingle");
			//TSnap::PlotHops(kron, "kronSingle");
			//PrintLargestEigenVal(kron, TFile, "Kron");
			if (IsDir == "false" && (MinMaxOutDeg != MinMaxInDeg || MaxMaxOutDeg != MaxMaxInDeg)){
				Error("GenKron", "Violation of reciprocity for undirected graph");
			}
			TFile << "Maximum output degree in kron graph: " << "from " << MinMaxOutDeg << " to " << MaxMaxOutDeg << " (average: " << (double)AvgMaxOutDeg / (double)NKron << ")" << endl;
			TFile << "Maximum input degree in kron graph: " << "from " << MinMaxInDeg << " to " << MaxMaxInDeg << " (average: " << (double)AvgMaxInDeg / (double)NKron << ")" << endl;
		}
		AddDegreesStat(KronDegAvgIn, SamplesIn, Kron, true);
		AddDegreesStat(KronDegAvgOut, SamplesOut, Kron, false);
	}
	Sec /= NKron;

	if (IsSampled == "true"){
		GetAvgDegreeStat(KronDegAvgIn, SamplesIn);
		GetAvgDegreeStat(KronDegAvgOut, SamplesOut);
	}
	else {
		GetAvgDegreeStat(KronDegAvgIn, NKron);
		GetAvgDegreeStat(KronDegAvgOut, NKron);
	}
	KronDegAvgIn.Sort();
	KronDegAvgOut.Sort();
	TFile << "Average time of generation of Kronecker product: " <<  Sec << endl;
}

void GetGraphFromAvgDistr(TFltPrV in_deg_avg_kron, PNGraph& t_pt)
{
	TIntV vec;
	for (size_t i = 0; i < in_deg_avg_kron.Len(); i++){
		// for all nodes with the same degree
		for (int j = 0; j < static_cast<int>(in_deg_avg_kron[i].Val2); j++)
			// add degree to vec
			vec.Add(static_cast<int>(in_deg_avg_kron[i].Val1));
	}
	TExeTm execTime;
	PUNGraph G = TSnap::GenConfModel(vec);
	TFile << "Time of configuration model: " <<  execTime.GetTmStr() << endl;
	t_pt = TSnap::ConvertGraph<PNGraph>(G);

}

// create string with parameters of model graph as an input to GenKron() function
TStr GetModelParamsStr(const PNGraph& G, const TStr& IsDir){
	TStr ModelParamsStr = " -nodes:"; 
	ModelParamsStr += to_string((long long)G->GetNodes()).c_str();
	ModelParamsStr += " -e:";
	ModelParamsStr += to_string((long long)G->GetEdges()).c_str();
	int MaxInDeg = GetMaxMinDeg(G, IsDir, "true", "true"); int MaxOutDeg = GetMaxMinDeg(G, IsDir, "false", "true");
	int	MinInDeg = GetMaxMinDeg(G, IsDir, "true", "false"); int MinOutDeg = GetMaxMinDeg(G, IsDir, "false", "false");
	ModelParamsStr += " -modelinmax:";
	ModelParamsStr += to_string((long long)MaxInDeg).c_str();
	ModelParamsStr += " -modeloutmax:";
	ModelParamsStr += to_string((long long)MaxOutDeg).c_str();
	ModelParamsStr += " -modelinmin:";
	ModelParamsStr += to_string((long long)MinInDeg).c_str();
	ModelParamsStr += " -modeloutmin:";
	ModelParamsStr += to_string((long long)MinOutDeg).c_str();
	return ModelParamsStr;
}

void GetGraphs(const vector <TStr>& Parameters, const TStr& ModelGen, const TStr&ModelPlt)
{
	PNGraph G;
	GetModel(Parameters[GRAPHGEN], G);
	
	const TStr& NEigenStr = Parameters[NEIGEN];
	int NEigen = NEigenStr.GetInt();
	if (NEigen != 0)
		PlotEigen(G, NEigenStr, "model");
	
	
	TFltPrV MDegIn, MDegOut;
	TSnap::GetInDegCnt(G, MDegIn);
	TSnap::GetOutDegCnt(G, MDegOut);

	PlotDegrees(Parameters, MDegIn, MDegOut, "model");
	PlotMetrics(Parameters, G, "model", TFile);
			
	if (ModelGen == "model+kron"){
		// generate (or read) Kronecker initiator matrix
		TKronMtx FitMtxM;
		if (!GetMtx(Parameters[MTXGEN], FitMtxM))
			GenNewMtx(G, Parameters[KRONFIT], FitMtxM);

		int ModelNodes = G->GetNodes(), ModelEdges = G->GetEdges();

		Env = TEnv(Parameters[KRONGEN], TNotify::NullNotify);
		TStr IsDir = Env.GetIfArgPrefixStr("-isdir:", "false", "Produce directed graph (true, false)");
		const TInt NIter = Env.GetIfArgPrefixInt("-i:", 1, "Number of iterations of Kronecker product");
		const TStr ScaleMtx = Env.GetIfArgPrefixStr("-scalemtx:", "false", "Scale initiator matrix (yes/no)");

		int MaxModelOutDeg = MDegOut[MDegOut.Len()-1].Val1;
		// scaling of initiator matrix
		// 3rd and 4th parameters - the maximum in/out degree of model graph
		ScaleFitMtx(ModelNodes, ModelEdges, MDegIn[MDegIn.Len()-1].Val1, MaxModelOutDeg, FitMtxM, NIter, IsDir, ScaleMtx, TFile);
		double ScalingCoeff = 0.0;
		if (ScaleMtx == "true"){
			TExeTm T;
			T.Tick();
			ScalingCoeff = GetScalingCoefficient(MDegIn, MDegOut, FitMtxM, NIter, IsDir);
			TFile << "Time of getting scaling coefficient: " << T.GetSecs() << endl;
			TFile << "Scailing coefficient: " << ScalingCoeff << endl;
			cout << "Scailing coefficient: " << ScalingCoeff << endl;
			system("pause");
			ScaleFitMtx(FitMtxM, NIter, ModelNodes, MaxModelOutDeg + MaxModelOutDeg * ScalingCoeff, IsDir);
		}

		// in and out average degrees of Kronecker graphs
		TFltPrV KronDegAvgIn, KronDegAvgOut;
		
		TStr KronParameters = GetModelParamsStr(G, IsDir);
		vector<BaseDiap> BaseDiaps; vector<int> Prev;

		PNGraph Kron;
		GenKron(Parameters[KRONGEN] + KronParameters, FitMtxM, KronDegAvgIn, KronDegAvgOut, BaseDiaps, Kron);

		PlotDegrees(Parameters, KronDegAvgIn, KronDegAvgOut, "kron");
		
		PNGraph  K;
		GetGraphFromAvgDistr(KronDegAvgOut, K);
		PlotMetrics(Parameters, K, "kron", TFile);
	}
}

// get FitMtx and scaling coefficient from small model
void GetFitMtxFromMS(TKronMtx& FitMtxM, TFlt& ScalingCoeff, vector<BaseDiap>& BaseDiaps, const vector<TStr>& Parameters){
	PNGraph G;
	GetModel(Parameters[GRAPHGEN], G);
	int ModelNodes = G->GetNodes(), ModelEdges = G->GetEdges();
	TFltPrV MDegIn, MDegOut;
	TSnap::GetInDegCnt(G, MDegIn);
	TSnap::GetOutDegCnt(G, MDegOut);
	const TFltPr MaxOut = MDegOut[MDegOut.Len()-1], MinOut = MDegOut[0];
	int MinDeg = MinOut.Val1, MaxDeg = MaxOut.Val1;
	TIntPr InDegR(MinDeg, MaxDeg), OutDegR(MinDeg, MaxDeg);
	// generate (or read) Kronecker initiator matrix
	if (!GetMtx(Parameters[MTXGEN], FitMtxM))
		GenNewMtx(G, Parameters[KRONFIT], FitMtxM);
	Env = TEnv(Parameters[KRONGEN], TNotify::NullNotify);
	TStr IsDir = Env.GetIfArgPrefixStr("-isdir:", "false", "Produce directed graph (true, false)");
	const TInt NIter = Env.GetIfArgPrefixInt("-i:", 1, "Number of iterations of Kronecker product");
	const TStr ScaleMtx = Env.GetIfArgPrefixStr("-scalemtx:", "false", "Scale initiator matrix (yes/no)");
	int MaxModelOutDeg = MDegOut[MDegOut.Len()-1].Val1;
	// scaling of initiator matrix
	// 3rd and 4th parameters - the maximum in/out degree of model graph
	ScaleFitMtx(ModelNodes, ModelEdges, MDegIn[MDegIn.Len()-1].Val1, MaxModelOutDeg, FitMtxM, NIter, IsDir, ScaleMtx, TFile);
	if (ScaleMtx == "true"){
		TKronMtx NewMtx(FitMtxM);
		NewMtx.Dump(TFile);
		TExeTm T;
		T.Tick();
		// if scaling coefficient is not given
		if (ScalingCoeff == -1.00)
			ScalingCoeff = GetScalingCoefficient(MDegIn, MDegOut, NewMtx, NIter, IsDir);
		TFile << "Time of getting scaling coefficient: " << T.GetSecs() << endl;
		TFile << "Scailing coefficient: " << ScalingCoeff << endl;
		cout << "Scailing coefficient: " << ScalingCoeff << endl;
		TFltPrV KronDeg;
		const TInt NKron = 15;
		ScaleFitMtx(NewMtx, NIter, ModelNodes, MaxModelOutDeg + MaxModelOutDeg * ScalingCoeff, IsDir);
		cout << "NewMtx:" << endl;
		NewMtx.Dump();
		GetAvgKronDeg(NewMtx, NIter, IsDir, NKron, OutDegR, KronDeg);
		TFltPrV RelDiffNonCum;
		//GetRelativeDiff(MDegOut, KronDeg, RelDiffNonCum);
		PrintDegDistr(KronDeg, "KronAvg.tab");
		//PrintDegDistr(MDegOut, "Model.tab");
		//PrintRelDiff(RelDiffNonCum, "RelDiff.tab");
		GetBaseDiaps(MDegIn, KronDeg, BaseDiaps);
		PrintBaseDiaps(BaseDiaps, "BaseDiaps.tab");
		/*PNGraph Kron; 
		KroneckerGen(Kron, NewMtx, NIter, IsDir, OutDegR, OutDegR, 0);
		Rewire(Kron, SmoothedDiaps, OutDegR);
		KronDeg.Clr();
		TSnap::GetInDegCnt(Kron, KronDeg);
		PlotPoints(KronDeg,	KronDeg, "KronScaled", "all");*/
		//system("pause");
	}
}

TFlt GetDegPart(const TInt& Deg, const TInt& DegMax, const TInt& DegMin){
	return (Deg - DegMin) / static_cast<double>(DegMax - DegMin);
}



// get average estimates of out-degree of NKron Kronecker graphs
void GetAvgKronDeg(const TKronMtx& NewMtx, const TInt& NIter, const TStr& IsDir, const TInt& NKron, const TIntPr& ModelDegR, TFltPrV& KronDeg){
	TIntPrV Samples;
	for (int i = 0; i < NKron; i++){
		PNGraph G;
		KroneckerGen(G, NewMtx, NIter, IsDir, ModelDegR, ModelDegR, 0);
		/*char Buf[15];
		itoa(i, Buf, 10);
		TStr Str(Buf);
		PrintDegDistr(G, "Kron" + Str);*/
		// ! out degrees
		AddDegreesStat(KronDeg, Samples, G, false);
	}
	GetAvgDegreeStat(KronDeg, NKron);
	KronDeg.Sort();
}

// generates Kronecker model using configuration model of small model network
// and compare it to big network
void KroneckerByConf(vector<TStr> CommandLineArgs){
	Try
	Env = TEnv(CommandLineArgs[KRONTEST], TNotify::NullNotify);
	// generation of big model and its Kronecker product is required
	const TStr ModelGen = Env.GetIfArgPrefixStr("-mgen:", "model", "Generation of big model and/or its Kronecker product (model, kron, model+kron)");
	// generation of big model and its Kronecker product is required
	const TStr ModelPlt = Env.GetIfArgPrefixStr("-mplt:", "model", "Plotting of big model and/or its Kronecker product (model, kron, model+kron)");
	// generation of big model and its Kronecker product is required
	const TStr MSGen = Env.GetIfArgPrefixStr("-msgen:", "model+kron", "Generation of small model and/or its Kronecker product (model, kron, model+kron)");
	// generation of big model and its Kronecker product is required
	const TStr MSPlt = Env.GetIfArgPrefixStr("-msplt:", "kron", "Plotting of small model and/or its Kronecker product (model, kron, model+kron)");
	// time estimates file name
	const TStr StatFile = Env.GetIfArgPrefixStr("-ot:", "stat.tab", "Name of output file with statistics");
	// is it required to generate big Kronecker graphs using FitMtx estimated on small graph (yes/no)
	const TStr FromMS = Env.GetIfArgPrefixStr("-fromms:", "yes", "Create big graph using FitMtx estimated for small graph (yes/none)");
	// scaling coefficient (if the parameter is given, scaling coefficient will not be calculated)
	TFlt ScalingCoeff = Env.GetIfArgPrefixFlt("-sc:", -1, "Scaling coefficient");

	CheckParams(ModelGen, ModelPlt);
	CheckParams(MSGen, MSPlt);

	TFile = OpenFile(StatFile.CStr());

	PyInit("PySettings.txt");

	if (FromMS == "yes"){
		TKronMtx FitMtx;
		vector <TStr> Parameters;
		GetParameters(CommandLineArgs, "Small", Parameters);
		vector<BaseDiap> BaseDiaps; 
		GetFitMtxFromMS(FitMtx, ScalingCoeff, BaseDiaps, Parameters);
		Parameters.clear();
		GetParameters(CommandLineArgs, "Model", Parameters);
		// generate big graph and plot its degrees
		PNGraph G;
		GetModel(Parameters[GRAPHGEN], G);
		Env = TEnv(Parameters[GRAPHGEN], TNotify::NullNotify);
		const TInt ModelNodes = Env.GetIfArgPrefixInt("-n:", 1024, "Model nodes count");
		Env = TEnv(Parameters[KRONGEN], TNotify::NullNotify);
		const TInt NIter = Env.GetIfArgPrefixInt("-i:", 1, "Number of iterations of Kronecker product");
		if (pow(FitMtx.GetDim(), (double)NIter) != ModelNodes)
			Error("KroneckerByConf", "ModelNodes != 2^NIter");

		TFltPrV MDegIn, MDegOut;
		TSnap::GetInDegCnt(G, MDegIn);
		TSnap::GetOutDegCnt(G, MDegOut);
		PlotDegrees(Parameters, MDegIn, MDegOut, "model");
		PlotMetrics(Parameters, G, "model", TFile);


		TStr IsDir = Env.GetIfArgPrefixStr("-isdir:", "false", "Produce directed graph (true, false)");
		
		// in and out average degrees of Kronecker graphs
		TFltPrV KronDegAvgIn, KronDegAvgOut;
		TStr KronParameters = GetModelParamsStr(G, IsDir);
		double MaxModelOutDeg = MDegOut[MDegOut.Len()-1].Val1;
		double RequiredDeg = MaxModelOutDeg + MaxModelOutDeg * ScalingCoeff;
		ScaleFitMtx(G->GetNodes(), G->GetEdges(), RequiredDeg, RequiredDeg, FitMtx, NIter, IsDir, "true", TFile);
		PNGraph Kron;
		GenKron(Parameters[KRONGEN] + KronParameters, FitMtx, KronDegAvgIn, KronDegAvgOut, BaseDiaps, Kron);
		//TestModelKronRelation(G, Kron);
		//PlotDegrees(Parameters, KronDegAvgIn, KronDegAvgOut, "kron");
		PlotPoints(KronDegAvgIn, KronDegAvgOut, "KronModel", "all");
		PNGraph  K;
		GetGraphFromAvgDistr(KronDegAvgOut, K);
		PlotMetrics(Parameters, K, "kron", TFile);
	}

	else {
		if (MSGen != "none"){
			vector <TStr> Parameters;
			GetParameters(CommandLineArgs, "Small", Parameters);
			GetGraphs(Parameters, MSGen, MSPlt);
		}

		if (ModelGen != "none")
		{
			vector <TStr> Parameters;
			GetParameters(CommandLineArgs, "Model", Parameters);
			GetGraphs(Parameters, ModelGen, ModelPlt);
		}
	}
	
	TFile.close();

	Py_Finalize();

	Catch
}

