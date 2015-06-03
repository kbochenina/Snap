#include "stdafx.h"
#include <iostream>
#include "GenPy.h"
#include "Eigen.h"
#include "Error.h"

ofstream TFile;


int GraphGen(const TStr args, PNGraph &GD){
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


int KroneckerGen(const TInt NIter, const TKronMtx& FitMtx, PNGraph& out, const TStr& OutFNm, const TIntPr& InDegR, const TIntPr& OutDegR, const TStr& IsDir, double ModelClustCf){
	Env.PrepArgs(TStr::Fmt("Kronecker graphs. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
	const TKronMtx& SeedMtx = FitMtx;
	printf("\n*** Seed matrix:\n");
	SeedMtx.Dump();
	printf("\n*** Kronecker:\n");
	// slow but exact O(n^2) algorightm
	//out = TKronMtx::GenKronecker(SeedMtx, NIter, false, 0); 
	// fast O(e) approximate algorithm
	// if we need to save clustering coefficient
	bool Dir = IsDir == "true" ? true: false;
	ModelClustCf = 0;
	// if we don't have constraints on degrees, run basic algorithm
	if (InDegR.Val1 == numeric_limits<int>::lowest() && InDegR.Val2 == INT_MAX && OutDegR.Val1 == numeric_limits<int>::lowest() && INT_MAX)
		out = TKronMtx::GenFastKronecker(SeedMtx, NIter, Dir, 0, ModelClustCf);
	else {
		TKronMtx::GenFastKronecker(SeedMtx, NIter, Dir, 0, InDegR, OutDegR, out, ModelClustCf);
	}

	//RemoveZeroDegreeNodes(out);

	// save edge list
	//TSnap::SaveEdgeList(out, OutFNm, TStr::Fmt("Kronecker Graph: seed matrix [%s]", FitMtx.GetMtxStr().CStr()));
	Catch
		printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
	
	return 0;
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
void GetModel(const TStr& args, PNGraph& G, const TStr& name, const TStr& Plt){
	Env = TEnv(args, TNotify::StdNotify);
	const TStr Gen = Env.GetIfArgPrefixStr("-g:", "gen", "How to get model graph: read, gen, deg, genpy");
	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "", "Input graph file (single directed edge per line)");
	
	TExeTm execTime;
	if (Gen == "gen")
		GraphGen(args, G);
	else if (Gen == "read")
		ReadPNGraphFromFile(InFNm, G);
	else if (Gen == "deg"){
		const TInt NodesCount = Env.GetIfArgPrefixInt("-n:", 1024, "Nodes count");
		const TFlt Gamma = Env.GetIfArgPrefixFlt("-gamma:", 2.0, "Gamma");
		TIntV DegSeqV;
		GetPowerLawDistrib(DegSeqV, NodesCount, Gamma);
		TExeTm t;
		PUNGraph GU = TSnap::GenDegSeq(DegSeqV);
		TFile << "Time of execution of configuration model: " <<  t.GetTmStr() << endl;
		G = TSnap::ConvertGraph<PNGraph>(GU);
	}
	if (Gen == "genpy")
	{
		PUNGraph GU;
		GenPy(GU, TFile, args);	
		G = TSnap::ConvertGraph<PNGraph>(GU);
	}
	if (Plt == "cum" || Plt == "all")
		SaveAndPlot(G, name.CStr(), true);
	if (Plt == "noncum" || Plt == "all")
		SaveAndPlot(G, name.CStr(), false);
	TFile << "Time of getting model: " <<  execTime.GetTmStr() << endl;
	TFile << "Model graph: " << G->GetNodes() << " nodes, " << G->GetEdges() << " edges\n";
	TIntV DegV;
	TSnap::GetDegSeqV(G, DegV);
	execTime.Tick();
	PUNGraph Conf = TSnap::GenConfModel(DegV);
	TFile << "Time of getting configuration model: " <<  execTime.GetTmStr() << endl;
	cout << "Undirected configuration model: " << Conf->GetNodes() << " nodes, " << Conf->GetEdges() << " edges\n";
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

int GetMaxDeg(const PNGraph& G)
{
	TIntPrV DegCnt;
	TSnap::GetDegCnt(TSnap::ConvertGraph<PUNGraph>(G), DegCnt);
	// sort in descending order
	DegCnt.Sort(false);
	return DegCnt[0].Val1;
}

int GetExpectedModelEdges(const PNGraph& G, const int k, const TStr& order){
	int expectedEdges = 0;
	map<int, int> DegAcc;
	TIntV DegSeq; 
	TSnap::GetDegSeqV(G, DegSeq);
	for (int i = 0; i < DegSeq.Len(); i++){
		if (order == "linear")
			expectedEdges += DegSeq[i] * k;
		else if (order == "square")
			expectedEdges += DegSeq[i] * k * k; 
	}
	// each edge is considered twice for both vertices
	expectedEdges /= 2;
	return expectedEdges;
}

void ScaleFitMtxForEdges(TKronMtx& FitMtx, const TInt& NIter, const int& ExpectedModelEdges){
	int ExpectedNodes = FitMtx.GetNodes(NIter);
	TFile << "Expected nodes: " << ExpectedNodes << " expected edges: " << FitMtx.GetEdges(NIter) << endl;
	double KronEdges = 0;
	while (abs (ExpectedModelEdges - KronEdges)  > 0.001 * ExpectedModelEdges){
		// after that there could be elements more that 1
		FitMtx.SetForEdgesNoCut(ExpectedNodes, ExpectedModelEdges);
		KronEdges = FitMtx.GetEdges(NIter);
	}
	TFile << "Scaled nodes " << FitMtx.GetNodes(NIter) << " scaled edges " << FitMtx.GetEdges(NIter) << endl;
	FitMtx.Dump(TFile);
	FitMtx.Normalize();
	cout << "Normalized matrix: \n";
	FitMtx.Dump(TFile);
}

void ScaleFitMtx(TKronMtx& FitMtx, const TInt& NIter, const int& InitModelNodes, const int& ExpectedModelEdges, const int& MaxModelDeg){
	ScaleFitMtxForEdges(FitMtx, NIter, ExpectedModelEdges);	
	TFile << "Before scaling " << endl;
	FitMtx.Dump(TFile);
	// check ceil()
	double ModelIter = ceil(log10((double)InitModelNodes) / log10((double)FitMtx.GetDim()));
	// rename function and variable
	int MinMaxDeg = FitMtx.GetMaxExpectedDeg(NIter);
	TFile << "Expected model maximum degree: " << MaxModelDeg << endl << "Expected Kronecker maximum degree: "<<  MinMaxDeg << endl;
	FitMtx.SetForMaxDeg(MaxModelDeg, ModelIter);
	TFile << "After scaling " << endl;
	FitMtx.Dump(TFile);
}

void GenKron(const TStr& args, TKronMtx& FitMtx, TFltPrV& inDegAvgKronM, TFltPrV& outDegAvgKronM, const PNGraph& G, const int NEigen, double ModelClustCf = 0.0){
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
	// output file name
	const TStr ScaleMtx = Env.GetIfArgPrefixStr("-scalemtx:", "true", "Scale init matrix to match number of edges");
	// output file name
	const TStr IsDir = Env.GetIfArgPrefixStr("-isdir:", "false", "Produce directed graph (true, false)");
	// restrictions to in- and out- degrees count for 1 vertex
	const TInt InMin = Env.GetIfArgPrefixInt("-inmin:", numeric_limits<int>::lowest(), "In-degree minimum");
	const TInt InMax = Env.GetIfArgPrefixInt("-inmax:", INT_MAX, "In-degree maximum");
	const TInt OutMin = Env.GetIfArgPrefixInt("-outmin:", numeric_limits<int>::lowest(), "Out-degree minimum");
	const TInt OutMax = Env.GetIfArgPrefixInt("-outmax:", INT_MAX, "Out-degree maximum");
	const TIntPr InDegR(InMin, InMax); const TIntPr OutDegR(OutMin, OutMax);
	
	float ModelNodes = G->GetNodes(), ModelEdges = G->GetEdges(), ExpectedNodes = FitMtx.GetNodes(NIter), ExpectedEdges = FitMtx.GetEdges(NIter);
	TFile << "Init model nodes: " << ModelNodes << ", init model edges: " << ModelEdges << endl;
	int ExpectedModelEdges = (ModelNodes == ExpectedNodes) ? ModelEdges : GetExpectedModelEdges(G, ExpectedNodes / ModelNodes, "linear");
	TFile << "Expected model nodes: " << ExpectedNodes << ", expected model edges: " << ExpectedModelEdges << endl;
	TFile << "Expected nodes: " << ExpectedNodes << ", expected edges: " << ExpectedEdges << endl;
	// check function
	int MaxModelDeg = GetMaxDeg(G);
	TFile << "Maximum degree in model graph: " << MaxModelDeg << endl;
	//TFile << "Maximum expected degree in kron graph: " << FitMtx.GetMinMaxPossibleDeg(NIter) << endl;

	if (ScaleMtx == "true"){
		ScaleFitMtx(FitMtx, NIter, ModelNodes, ExpectedModelEdges, MaxModelDeg);
		//ScaleFitMtxForEdges(FitMtx, NIter, ExpectedModelEdges);	
	}

	if (NEigen != 0)
		PrintEigen(FitMtx, NIter, NEigen);

	// Kronecker model of graph
	PNGraph kron;
	TIntPrV samplesIn, samplesOut;
	double sec = 0.0;
	
	int Nodes = FitMtx.GetNodes(NIter);
	int Edges = FitMtx.GetEdges(NIter);
	if (InDegR.Val1 * Nodes > Edges || OutDegR.Val1 * Nodes > Edges || InDegR.Val2 < InDegR.Val1 || OutDegR.Val2 < OutDegR.Val1) 
		Error("GenKron", "Constraints do not match to the number of edges");
	if (IsDir == "false" && (InDegR.Val1 != OutDegR.Val1 || InDegR.Val2 != OutDegR.Val2))
		Error("GenKron", "InDegR and OutDegR should be the same for undirected graph");
	if (IsDir == "false"){
		FitMtx.Dump(TFile);
		TFile << "Maximum expected degree in kron graph: " << FitMtx.GetMaxExpectedDeg(NIter) << endl;
	}
		

	TFltV KronEigen;
	int AverageMaxDeg = 0, MinMaxDeg = 0, MaxMaxDeg = 0;

	for (int i = 0; i < NKron; i++){
		execTime.Tick();
		// ModelClustCf!
		FitMtx.Dump();
		KroneckerGen(NIter, FitMtx, kron, OutFnm, InDegR, OutDegR, IsDir, ModelClustCf);
		sec += execTime.GetSecs();
		int MaxDeg = GetMaxDeg(kron);
		if (i == 0) MinMaxDeg = MaxDeg;
		else if (MaxDeg < MinMaxDeg) MinMaxDeg = MaxDeg;
		if (MaxDeg > MaxMaxDeg) MaxMaxDeg = MaxDeg;
		AverageMaxDeg += MaxDeg;
		printf("Nodes count: %d, nodes with non-zero degree %d, edges count %d\n max deg = %d\n", kron->GetNodes(), TSnap::CntNonZNodes(kron), kron->GetEdges(), MaxDeg);
		if (i == NKron - 1){
			//TFile << "Clustering coefficient: " << TSnap::GetClustCf(kron) << endl;
			//TSnap::PlotClustCf(kron,"kronSingle");
			//TSnap::PlotHops(kron, "kronSingle");
			//PrintLargestEigenVal(kron, TFile, "Kron");
			if (NEigen != 0){
				KronEigen.Clr();
				TSnap::PlotEigValRank(TSnap::ConvertGraph<PUNGraph>(kron), NEigen, "KronEigen", KronEigen);
				TFile << "Maximum eigenvalue in kron graph: " << KronEigen[0].Val << endl;
			}
			TFile << "Maximum degree in kron graph: " << "from " << MinMaxDeg << " to " << MaxMaxDeg << " (average: " << (double)AverageMaxDeg / (double)NKron << ")" << endl;
		}
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
	TFile << "Average time of generation of Kronecker product: " <<  sec << endl;
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

// Get graph size from GRAPHGEN string
TInt GetNFromGraphGenStr(const TStr& GraphGenStr){
	TInt PosBegin = GraphGenStr.SearchStr("-n:");
	TInt PosEnd = 0;
	while (PosEnd < PosBegin) 
		PosEnd = GraphGenStr.SearchStr(" ", PosEnd + 1);
	TStr NodesCountStr = GraphGenStr.GetSubStr(PosBegin + 3, PosEnd - 1);
	return NodesCountStr.GetInt();
}

// return GRAPHGEN string with NodesCount = NodesCount * N
TStr GetModifiedStr(const TStr& SrcStr, const TInt& N){
	TStr NewGen(SrcStr);
	TInt NodesCount = GetNFromGraphGenStr(NewGen);
	TInt NewNodesCount = N.Val * NodesCount.Val;
	NewGen.ChangeStr(NodesCount.GetStr(), NewNodesCount.GetStr());
	return NewGen;
}

// Get average scaling coefficient and vector of scaling coefficients from two samples of graphs
TFlt GetAvgScaleCf(const TFltV& ModelEigValV, const TFltV& G2EigValV, const TInt& N, TFltV& ScaleCfV){
	TFlt AvgScaleCf = 0.0;
	TInt NEigen = 0;
	for (int i = 0; i < ModelEigValV.Len(); i++){
		if (G2EigValV[i].Val < 0 || ModelEigValV[i].Val < 0) break;
		//TFlt ScaleRatio = (G2EigValV[i].Val / ModelEigValV[i].Val) / N;
		TFlt ScaleRatio = log10(G2EigValV[i].Val) / log10(ModelEigValV[i].Val);
		ScaleCfV.Add(ScaleRatio);
		AvgScaleCf += ScaleRatio;
		NEigen = NEigen + 1;
	}
	AvgScaleCf = AvgScaleCf / NEigen;

	//TFlt MaxEigenVal = ModelEigValV[0].Val;
	//TInt SampleCountSum = 0;
	//for (int i = 0; i < NEigen; i++){
	//	TInt SampleCount = ModelEigValV[i].Val / MaxEigenVal * NEigen;
	//	//printf("%f %d\n", ModelEigValV[i].Val / MaxEigenVal, SampleCount);
	//	AvgScaleCf += ScaleCfV[i] * SampleCount;
	//	SampleCountSum += SampleCount;
	//}
	//AvgScaleCf = AvgScaleCf / SampleCountSum;

	TFile << "AvgScaleCf: " << AvgScaleCf << endl;

	TFlt SD = 0.0;

	for (int i = 0; i < ScaleCfV.Len(); i++){
		TFile << "ScaleCfV[" << i + 1 << "] = " << ScaleCfV[i];
		TFile << " Deviation: " << sqrt(abs(pow(ScaleCfV[i].Val, 2.00) - pow(AvgScaleCf.Val, 2.00))) << endl;
		SD = SD + pow(ScaleCfV[i].Val - AvgScaleCf.Val, 2.00);
	}
	SD = sqrt(SD.Val / (NEigen.Val - 1));
	TFile << "SD: " << SD << endl;
	return AvgScaleCf;
}


void TestScalingEigen(const TInt& ScaleCount, const TInt& ScaleSize, const TFlt& AvgScaleCf, const TStr& GraphGenStr, const TFltV& ModelEigValV, const TInt& ModelEdges, const TFltV & ScaleCfV){
	// initial nodes count
	TInt InitN = GetNFromGraphGenStr(GraphGenStr);
	TInt NewN = InitN;
	TInt NEigen = ScaleCfV.Len();
	for (TInt i = 0; i < ScaleCount; i++){
		PNGraph G;
		NewN = NewN * ScaleSize;
		printf("NewN = %d\n", NewN);
		TInt SizeRatio = NewN / InitN;
		TStr NewStr = GetModifiedStr(GraphGenStr, SizeRatio);
		// produce and plot model graph
		GetModel(NewStr, G, "time.tab", "none");
		TFltV GEigValV;
		TSnap::PlotEigValRank(TSnap::ConvertGraph<PUNGraph>(G), NEigen, "ModelEigen" + i.GetStr(), GEigValV);
		//TFlt ScaleCoeff = pow(AvgScaleCf * ScaleSize, i + 1);
		//printf("ScaleCoeff: %f\n", ScaleCoeff.Val);
		TFlt ScaleCoeff = AvgScaleCf * SizeRatio / ScaleSize;
		TFltV ApproxEigValV;
		//TFlt ApproxE = ModelEdges * ScaleCoeff;
		TFlt ApproxE = pow(ModelEdges, ScaleCoeff);
		ofstream F = OpenFile("ApproxEigen" + i.GetStr() + ".tab");
		F << "#" << endl;
		F << " ApproxEigen. G(" << NewN << ", " << ApproxE.Val << "). ";
		printf("NEigen = %d\n", NEigen);
		for (int j = 0; j < NEigen; j++){
			//ScaleCoeff = pow(ScaleCfV[j] * ScaleSize, i + 1);
			//ScaleCoeff = ScaleCfV[j] * SizeRatio / ScaleSize;
			//TFlt ApproxEigenVal = ModelEigValV[j] * ScaleCoeff;
			TFlt ApproxEigenVal = ModelEigValV[j];
			printf("ApproxEigenVal = %f ScaleCfV[j] = %f SizeRatio / ScaleSize = %d\n", ApproxEigenVal, ScaleCfV[j], SizeRatio / ScaleSize);
			for (int i = 0; i < SizeRatio / ScaleSize; i++){
				ApproxEigenVal = pow(ApproxEigenVal, ScaleCfV[j]);
				printf("ApproxEigenVal = %f ScaleCfV[j] = %f SizeRatio / ScaleSize = %d\n", ApproxEigenVal, ScaleCfV[j], SizeRatio / ScaleSize);
			}
			if (j == 0){
				F << "Largest eig val = " << ApproxEigenVal.Val << endl << "#" << endl << "# Rank	Eigen value" << endl;
			}
			ApproxEigValV.Add(ApproxEigenVal);
			F << j + 1 << " " << ApproxEigenVal.Val << endl;
		}
		F.close();
	}
}

void GetGraphs(vector <TStr>& parameters, vector<TFltPrV>& distrIn, vector<TFltPrV>& distrOut, TStrV& names, const TStr& ModelGen, const TStr&ModelPlt)
{
	PNGraph G;
	const TStr& name = parameters[NAME];
	const TStr& Plt = parameters[PLT];
	const TStr& PType = parameters[PTYPE];
	const TStr& NEigenStr = parameters[NEIGEN];
	const TStr& ScaleSizeStr = parameters[SCALE_SIZE];
	const TStr& ScaleCountStr = parameters[SCALE_COUNT];

	const TInt& ScaleSize = ScaleSizeStr.GetInt();
	const TInt& ScaleCount = ScaleCountStr.GetInt();

	GetModel(parameters[GRAPHGEN], G, name, parameters[PLT]);
	int NEigen = NEigenStr.GetInt();
	if (NEigen != 0){
		TFltV ModelEigValV;
		TSnap::PlotEigValRank(TSnap::ConvertGraph<PUNGraph>(G), NEigenStr.GetInt(), "ModelEigen", ModelEigValV);
		TFile << "Maximum eigenvalue in model graph: " << ModelEigValV[0].Val << endl;
	}

	/*TStr ModifiedStr = GetModifiedStr(parameters[GRAPHGEN], ScaleSize);
	PNGraph G2;
	GetModel(ModifiedStr, G2, name, parameters[PLT]);
	TFltV G2EigValV;
	TSnap::PlotEigValRank(TSnap::ConvertGraph<PUNGraph>(G2), NEigenStr.GetInt(), "G2Eigen", G2EigValV);
	TFltV ScaleCfV;
	TFlt AvgScaleCf = GetAvgScaleCf(ModelEigValV, G2EigValV, ScaleSize, ScaleCfV);

	TestScalingEigen(ScaleCount, ScaleSize, AvgScaleCf, parameters[GRAPHGEN], ModelEigValV, G->GetEdges(), ScaleCfV);*/

	double ModelClustCf = 0.0;

	if ( PType == "exp" || PType == "all" )
	{
		TFltPrV mDegIn, mDegOut;
		TSnap::GetInDegCnt(G, mDegIn);
		TSnap::GetOutDegCnt(G, mDegOut);
		distrIn.push_back(mDegIn); distrOut.push_back(mDegOut); names.Add(name + "Sparse");
		//ModelClustCf = TSnap::GetClustCf(G);
		TExeTm execTime;
		//TFile << "Clustering coefficient: " << ModelClustCf << endl;
		//TSnap::PlotClustCf(G, name);
		//TSnap::PlotHops(G, name);
		TFile << "Time of calculating the metrics: " << execTime.GetTmStr() << endl;
	}

	
	if (ModelGen == "model+kron"){
		// generate Kronecker initiator matrix using big graph
		TKronMtx FitMtxM;
		if (!GetMtx(parameters[MTXGEN], FitMtxM))
			GenNewMtx(G, parameters[KRONFIT], FitMtxM);

		// in and out average degree distribution for kronM (non-accumulated)
		TFltPrV inDegAvgKron, outDegAvgKron;
		
		GenKron(parameters[KRONGEN], FitMtxM, inDegAvgKron, outDegAvgKron, G, NEigenStr.GetInt());
		if ( PType == "full" || PType == "all" ){
		PlotPoints(inDegAvgKron, outDegAvgKron, name, Plt);
		}
		if ( PType == "exp" || PType == "all"){
			distrIn.push_back(inDegAvgKron); distrOut.push_back(outDegAvgKron); names.Add("kron" + name + "Sparse");
		}
		//PNGraph  K;
		//TExeTm execTime;
		//GetGraphFromAvgDistr(outDegAvgKron, K);
		//TFile << "Clustering coefficient: " << TSnap::GetClustCf(K) << endl;
		//TSnap::PlotClustCf(K, "kron" + name);
		//TSnap::PlotHops(K, "kron" + name);
		//TFile << "Time of calculating the metrics: " << execTime.GetTmStr() << endl;
	}
	
	
}



// generates Kronecker model using configuration model of small model network
// and compare it to big network
void KroneckerByConf(vector<TStr> commandLineArgs){
	Try
	Env = TEnv(commandLineArgs[KRONTEST], TNotify::StdNotify);
	// type of plots
	const TStr Plt = Env.GetIfArgPrefixStr("-plt:", "all", "Type of plots (cum, noncum, all)");
	// full - all points of distrib will be plotted; expbin - exponential binning
	const TStr PType = Env.GetIfArgPrefixStr("-ptype:", "all", "How to plot (full, expbin, all)");
	// radix of binning
	const TInt BinRadix = Env.GetIfArgPrefixInt("-bin:", 2, "Radix for exponential binning");
	// time estimates file name
	const TStr TimeFile = Env.GetIfArgPrefixStr("-ot:", "time.tab", "Name of output file with time estimates");
	// generation of big model and its Kronecker product is required
	const TStr ModelGen = Env.GetIfArgPrefixStr("-mgen:", "model", "Generation of big model and/or its Kronecker product (model, kron, model+kron)");
	// generation of big model and its Kronecker product is required
	const TStr ModelPlt = Env.GetIfArgPrefixStr("-mplt:", "model", "Plotting of big model and/or its Kronecker product (model, kron, model+kron)");
	// generation of big model and its Kronecker product is required
	const TStr MSGen = Env.GetIfArgPrefixStr("-msgen:", "model+kron", "Generation of small model and/or its Kronecker product (model, kron, model+kron)");
	// generation of big model and its Kronecker product is required
	const TStr MSPlt = Env.GetIfArgPrefixStr("-msplt:", "kron", "Plotting of small model and/or its Kronecker product (model, kron, model+kron)");
	// number of eigenvalues to investigate
	const TStr NEigen = Env.GetIfArgPrefixStr("-neigen:", "10", "Number of eigenvalues");
	// scale size
	const TStr ScaleSize = Env.GetIfArgPrefixStr("-ss:", "2", "Scale size");
	// scale count
	const TStr ScaleCount = Env.GetIfArgPrefixStr("-sc:", "5", "Scale count");

	CheckParams(ModelGen, ModelPlt);
	CheckParams(MSGen, MSPlt);

	TFile = OpenFile(TimeFile.CStr());

	vector<TFltPrV> distrIn, distrOut;
	TStrV names;

	PyInit("PySettings.txt");

	if (ModelGen != "none")
	{
		TFile << "Kronecker" << endl;
		vector <TStr> parameters;
		for (size_t i = 1; i <= NPARCOPY; i++)
			parameters.push_back(commandLineArgs[i]);
		parameters.push_back(PType); parameters.push_back(Plt); parameters.push_back("Model"); parameters.push_back(NEigen); parameters.push_back(ScaleSize); parameters.push_back(ScaleCount);
		GetGraphs(parameters, distrIn, distrOut, names, ModelGen, ModelPlt);
	}
	
	if (MSGen != "none"){
		TFile << "Kronecker from reduced size" << endl;
		vector <TStr> parameters;
		for (size_t i = NPARCOPY + 1; i <= 2 * NPARCOPY; i++)
			parameters.push_back(commandLineArgs[i]);
		parameters.push_back(PType); parameters.push_back(Plt); parameters.push_back("Small"); parameters.push_back(NEigen); parameters.push_back(ScaleSize); parameters.push_back(ScaleCount);
		GetGraphs(parameters, distrIn, distrOut, names, MSGen, MSPlt);
	}

	//system("pause");
	PlotSparse(distrIn, names, true, Plt, BinRadix);
	PlotSparse(distrOut, names, false, Plt, BinRadix);

	TFile.close();

	Py_Finalize();

	Catch
}

