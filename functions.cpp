#include "stdafx.h"
#include <iostream>
#include "GenPy.h"

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

int KroneckerGen(const TInt NIter, const TKronMtx& FitMtx, PNGraph& out, const TStr& OutFNm, const TIntPr& InDegR = TIntPr(-1,-1), const TIntPr& OutDegR = TIntPr(-1,-1)){
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
	// if we don't have constraints on degrees, run basic algorithm
	if (InDegR.Val1 == numeric_limits<int>::lowest() && InDegR.Val2 == INT_MAX && OutDegR.Val1 == numeric_limits<int>::lowest() && INT_MAX)
		out = TKronMtx::GenFastKronecker(SeedMtx, NIter, true, 0);
	else {
		out = TKronMtx::GenFastKronecker(SeedMtx, NIter, true, 0, InDegR, OutDegR);
	}

	//RemoveZeroDegreeNodes(out);

	// save edge list
	TSnap::SaveEdgeList(out, OutFNm, TStr::Fmt("Kronecker Graph: seed matrix [%s]", FitMtx.GetMtxStr().CStr()));
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

void GenKron(const TStr& args, TKronMtx& FitMtx, TFltPrV& inDegAvgKronM, TFltPrV& outDegAvgKronM, int modelNodes = 0, int modelEdges = 0){
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
	// restrictions to in- and out- degrees count for 1 vertex
	const TInt InMin = Env.GetIfArgPrefixInt("-inmin:", numeric_limits<int>::lowest(), "In-degree minimum");
	const TInt InMax = Env.GetIfArgPrefixInt("-inmax:", INT_MAX, "In-degree maximum");
	const TInt OutMin = Env.GetIfArgPrefixInt("-outmin:", numeric_limits<int>::lowest(), "Out-degree minimum");
	const TInt OutMax = Env.GetIfArgPrefixInt("-outmax:", INT_MAX, "Out-degree maximum");
	const TIntPr InDegR(InMin, InMax); const TIntPr OutDegR(OutMin, OutMax);
	
	
	// Kronecker model of graph
	PNGraph kron;
	TIntPrV samplesIn, samplesOut;
	double sec = 0.0;
	
	if (ScaleMtx == "true")
	{
		double expectedNodes = pow(static_cast<double>(FitMtx.GetDim()), NIter);
		cout << "Model nodes " << modelNodes << " model edges " << modelEdges << endl;
		int expectedEdges = expectedNodes / modelNodes * modelEdges;
		cout << "Expected nodes " << expectedNodes << " Expected edges " << expectedEdges << endl;
		FitMtx.SetForEdges(expectedNodes, expectedEdges);
		cout << "Scaled nodes " << FitMtx.GetNodes(NIter) << " scaled edges " << FitMtx.GetEdges(NIter) << endl;
	}

	try {
		int Nodes = FitMtx.GetNodes(NIter);
		int Edges = FitMtx.GetEdges(NIter);
		if (InDegR.Val1 * Nodes > Edges || OutDegR.Val1 * Nodes > Edges || InDegR.Val2 < InDegR.Val1 || OutDegR.Val2 < OutDegR.Val1) 
			throw exception("FastKronecker() exception. Constraints do not match to the number of edges");
	}
	catch (const exception e) {
		cerr << e.what() << endl;
		system("pause");
		exit(0);
	}

	for (int i = 0; i < NKron; i++){
		execTime.Tick();
		KroneckerGen(NIter, FitMtx, kron, OutFnm, InDegR, OutDegR);
		sec += execTime.GetSecs();
		printf("Nodes count: %d, nodes with non-zero degree: %d\n", kron->GetNodes(), TSnap::CntNonZNodes(kron));
		if (i == NKron - 1){
			TFile << "Clustering coefficient: " << TSnap::GetClustCf(kron) << endl;
			TSnap::PlotClustCf(kron,"kronSingle");
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

void GetGraphs(vector <TStr>& parameters, vector<TFltPrV>& distrIn, vector<TFltPrV>& distrOut, TStrV& names, const TStr& ModelGen, const TStr&ModelPlt)
{
	PNGraph G;
	const TStr& name = parameters[NAME];
	const TStr& Plt = parameters[PLT];
	const TStr& PType = parameters[PTYPE];

	GetModel(parameters[GRAPHGEN], G, name, parameters[PLT]);

	if ( PType == "exp" || PType == "all" )
	{
		TFltPrV mDegIn, mDegOut;
		TSnap::GetInDegCnt(G, mDegIn);
		TSnap::GetOutDegCnt(G, mDegOut);
		distrIn.push_back(mDegIn); distrOut.push_back(mDegOut); names.Add(name + "Sparse");
		TExeTm execTime;
		TFile << "Clustering coefficient: " << TSnap::GetClustCf(G) << endl;
		TSnap::PlotClustCf(G, name);
		TSnap::PlotHops(G, name);
		TFile << "Time of calculating the metrics: " << execTime.GetTmStr() << endl;
	}

	
	if (ModelGen == "model+kron"){
		// generate Kronecker initiator matrix using big graph
		TKronMtx FitMtxM;
		if (!GetMtx(parameters[MTXGEN], FitMtxM))
			GenNewMtx(G, parameters[KRONFIT], FitMtxM);

		// in and out average degree distribution for kronM (non-accumulated)
		TFltPrV inDegAvgKron, outDegAvgKron;
		GenKron(parameters[KRONGEN], FitMtxM, inDegAvgKron, outDegAvgKron, G->GetNodes(), G->GetEdges());
		if ( PType == "full" || PType == "all" ){
		PlotPoints(inDegAvgKron, outDegAvgKron, name, Plt);
		}
		if ( PType == "exp" || PType == "all"){
			distrIn.push_back(inDegAvgKron); distrOut.push_back(outDegAvgKron); names.Add("kron" + name + "Sparse");
		}
		PNGraph  K;
		TExeTm execTime;
		GetGraphFromAvgDistr(inDegAvgKron, K);
		TFile << "Clustering coefficient: " << TSnap::GetClustCf(K) << endl;
		TSnap::PlotClustCf(K, "kron" + name);
		TSnap::PlotHops(K, "kron" + name);
		TFile << "Time of calculating the metrics: " << execTime.GetTmStr() << endl;
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

	CheckParams(ModelGen, ModelPlt);
	CheckParams(MSGen, MSPlt);

	TFile = OpenFile(TimeFile.CStr());

	vector<TFltPrV> distrIn, distrOut;
	TStrV names;

	PyInit();

	if (ModelGen != "none")
	{
		TFile << "Kronecker" << endl;
		vector <TStr> parameters;
		for (size_t i = 1; i <= NPARCOPY; i++)
			parameters.push_back(commandLineArgs[i]);
		parameters.push_back(PType); parameters.push_back(Plt); parameters.push_back("Model"); 
		GetGraphs(parameters, distrIn, distrOut, names, ModelGen, ModelPlt);
	}
	
	if (MSGen != "none"){
		TFile << "Kronecker from reduced size" << endl;
		vector <TStr> parameters;
		for (size_t i = NPARCOPY + 1; i <= 2 * NPARCOPY; i++)
			parameters.push_back(commandLineArgs[i]);
		parameters.push_back(PType); parameters.push_back(Plt); parameters.push_back("Small"); 
		GetGraphs(parameters, distrIn, distrOut, names, MSGen, MSPlt);
	}

	system("pause");
	PlotSparse(distrIn, names, true, Plt, BinRadix);
	PlotSparse(distrOut, names, false, Plt, BinRadix);

	TFile.close();

	Py_Finalize();

	Catch
}

