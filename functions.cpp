#include "stdafx.h"
#include <iostream>
#include "GenPy.h"
#include "Eigen.h"
#include "Error.h"
#include "Diaps.h"

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

void ScaleFitMtx(int ModelNodes, int ModelEdges, TFlt MaxDegInModel, TFlt MaxDegOutModel, TKronMtx& FitMtx, const TInt& NIter, TStr IsDir, const TStr& ScaleMtx){
	// check if ModelNodes and NIter are consistent
	if (pow(FitMtx.GetDim(), (double)NIter) != ModelNodes)
		Error("ScaleFitMtx", "ModelNodes is not equal to 2^NIter");
	// scale init matrix to match number of edges
	ScaleFitMtxForEdges(FitMtx, NIter, ModelEdges);	

	// equalizing B and C for undirected graphs
	if (IsDir == "false" && FitMtx.At(0,1) != FitMtx.At(1,0)){
		cout << "WARNING. B and C is not equal for an undirected graph. B and C will be eqialized" << endl;
		ScaleFitMtxForUnDir(FitMtx);
	}

	TFile << "Expected maximum degree in Kronecker graph: in " << FitMtx.GetMaxExpectedDeg(NIter, IsDir, true) << 
	",  out " << FitMtx.GetMaxExpectedDeg(NIter, IsDir, false) << endl;

	cout << "Diapason of degrees: from " << pow(FitMtx.GetMinPossibleDeg(), NIter) << " to " << pow(FitMtx.GetMaxPossibleDeg(), NIter) << endl
			<< "  degree: " << MaxDegOutModel << endl;

	// scale initiator matrix to match maximum degree
	if (ScaleMtx == "true"){
		if (IsDir == "true"){
			cout << "WARNING. Matrix will be scaled to match maximum output degree. Required functional is under development\n";
			IsDir = "false";
			ScaleFitMtxForUnDir(FitMtx);
		}
		ScaleFitMtx(FitMtx, NIter, ModelNodes, MaxDegOutModel, IsDir);
		FitMtx.Dump();
		TFile << "Expected maximum degree in Kronecker graph: in " << FitMtx.GetMaxExpectedDeg(NIter, IsDir, true) << 
			",  out " << FitMtx.GetMaxExpectedDeg(NIter, IsDir, false) << endl;
	}
}


void ScaleFitMtxForEdges(TKronMtx& FitMtx, const TInt& NIter, const int& ExpectedModelEdges){
	TFile << "Scaling initiator matrix for correct number of edges..." << endl;
	int ExpectedNodes = FitMtx.GetNodes(NIter);
	//TFile << "Expected nodes: " << ExpectedNodes << " expected edges: " << FitMtx.GetEdges(NIter) << endl;
	double KronEdges = 0;
	while (abs (ExpectedModelEdges - KronEdges)  > 0.001 * ExpectedModelEdges){
		// after that there could be elements more that 1
		FitMtx.SetForEdgesNoCut(ExpectedNodes, ExpectedModelEdges);
		KronEdges = FitMtx.GetEdges(NIter);
	}
	TFile << "Scaled Kronecker nodes " << FitMtx.GetNodes(NIter) << " scaled Kronecker edges " << FitMtx.GetEdges(NIter) << endl;
	//FitMtx.Dump(TFile);
	FitMtx.Normalize();
	/*TFile << "Normalized matrix: \n";
	FitMtx.Dump(TFile);*/
}

void ScaleFitMtxForUnDir(TKronMtx& FitMtx){
	TFile << "Equalizing B and C..." << endl;
	FitMtx.EqualizeBC();
}

void ScaleFitMtx(TKronMtx& FitMtx, const TInt& NIter, const int& InitModelNodes, const int& MaxModelDeg, const TStr& IsDir){
	//TFile << "Before scaling: " << endl;
	//FitMtx.Dump(TFile);
	// check ceil()
	double ModelIter = ceil(log10((double)InitModelNodes) / log10((double)FitMtx.GetDim()));
	// rename function and variable
	int MinMaxDeg = (FitMtx.GetMaxExpectedDeg(NIter.Val, IsDir, false) + 0.5);
	//TFile << "Maximum degree in model graph: " << MaxModelDeg << endl << "Expected Kronecker maximum degree: "<<  MinMaxDeg << endl;
	// ModelIter instead of NIter
	FitMtx.SetForMaxDeg(MaxModelDeg, ModelIter, "false", false);
	//TFile << "After scaling " << endl;
	//FitMtx.Dump(TFile);
}

double GetAvgDeviation(const TFltPrV& ModelDegCnt, const TFltPrV& KronDegCnt){
	double Dev = 0.0;
	bool SeqViewed = false;
	double MaxModelDeg = ModelDegCnt[ModelDegCnt.Len()-1].Val1, MaxKronDeg = KronDegCnt[KronDegCnt.Len()-1].Val1,
		MaxModelCount = ModelDegCnt[ModelDegCnt.Len()-1].Val2, MaxKronCount = KronDegCnt[KronDegCnt.Len()-1].Val2;
	if (MaxModelDeg == MaxKronDeg) 
		return (pow((MaxModelCount-MaxKronCount)/MaxModelCount,2));
	else {
		return (pow((MaxModelDeg-MaxKronDeg)/MaxModelDeg,2));
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

double GetBestCoeff(TFltPrV& ScalingResults){
	if (ScalingResults.Len() == 0) return 0;
	double ScalingCoeff = 0, Dev = ScalingResults[0].Val2;
	if (ScalingResults.Len() == 1) return ScalingResults[0].Val1;
	
	for (int i = 0; i < ScalingResults.Len()-1; i++){
		ScalingResults[i].Val2 = (ScalingResults[i].Val2 + ScalingResults[i+1].Val2) / 2;
	}

	for (int i = 0; i < ScalingResults.Len()-1; i++){
		if (ScalingResults[i].Val2 < Dev){
			Dev = ScalingResults[i].Val2;
			if (i == ScalingResults.Len() - 1)
				ScalingCoeff = ScalingResults[i].Val1;
			else 
				ScalingCoeff = (ScalingResults[i].Val1 + ScalingResults[i+1].Val1) / 2;
		}
	}
	return ScalingCoeff;
}

// estimate scaling coefficient
double GetScalingCoefficient(const TFltPrV& InDegCnt, const TFltPrV& OutDegCnt, const TKronMtx& FitMtxM, const TInt& NIter, const TStr& IsDir){
	// !!!
	//return 1.4;
	TKronMtx FitMtx(FitMtxM);
	double ScalingCoeff = 0;
	double ScalingStep = 0.2;
	bool DecFound = false;
	const int NKron = 40;
	const TFltPr MaxOut = OutDegCnt[OutDegCnt.Len()-1], MinOut = OutDegCnt[0];
	int MinDeg = MinOut.Val1, MaxDeg = MaxOut.Val1;
	TIntPr InDegR(MinDeg, MaxDeg), OutDegR(MinDeg, MaxDeg);
	//TIntPr InDegR(MinDeg, 1000000), OutDegR(MinDeg, 1000000);
	bool FirstTry = true; bool ToIncrease = true;
	double Dev = 0.0;
	TInt Count = 0;
	TFltPrV ScalingResults;
	vector<TFltPrV> KronDeg;

	while (!DecFound){
		if (!FirstTry){
			ScalingCoeff += ScalingStep;
			int DegToScale = (int) (MaxDeg + MaxDeg * ScalingCoeff);
			if (DegToScale <= 1 || !FitMtx.CanScaleToDeg(DegToScale, NIter)){ // !!!abs(ScalingCoeff) >= 1 ||
				return GetBestCoeff(ScalingResults);
			}
			// we use the assumption that ModelNodes = KronNodes
			ScaleFitMtx(FitMtx, NIter, (int)pow(2.00, NIter), (int) (MaxDeg + MaxDeg * ScalingCoeff), IsDir);
			/*TKronMtx NewMtx(FitMtx);
			ScaleFitMtx(NewMtx, NIter, (int)pow(2.00, NIter), (int) (MaxDeg + MaxDeg * ScalingCoeff), IsDir);*/
			FitMtx.Dump();
		}
		// get average statistics on degrees
		PNGraph Kron;
		TFltPrV KronDegAvgOut;
		// for the compatibility
		TIntPrV SamplesOut;
		for (int i = 0; i < NKron; i++){
			KroneckerGen(Kron, FitMtx, NIter, IsDir, InDegR, OutDegR, 0);
			AddDegreesStat(KronDegAvgOut, SamplesOut, Kron, false);
		}
		//GetAvgDegreeStat(KronDegAvgOut, NKron);
		for (int i = 0; i < KronDegAvgOut.Len(); i++) KronDegAvgOut[i].Val2 /= NKron;
		KronDegAvgOut.Sort();
		KronDeg.push_back(KronDegAvgOut);
		PlotPoints(KronDegAvgOut, KronDegAvgOut, "scale" + Count.GetStr(), "all");
		if (FirstTry){
			// if maximum kron degree > maximum model degree
			if (KronDegAvgOut[KronDegAvgOut.Len()-1] > MaxOut){
				ToIncrease = false;
				//!!!
				ScalingStep *= -1;
			}
			Dev = GetAvgDeviation(OutDegCnt, KronDegAvgOut);
			printf("Scaling coeff: %f Dev: %f\n", ScalingCoeff, Dev);
			//system("pause");
			FirstTry = false;
		}
		else {
			double CurrDev = GetAvgDeviation(OutDegCnt, KronDegAvgOut);
			TFltPr Res(ScalingCoeff, CurrDev);
			ScalingResults.Add(Res);
			printf("Scaling coeff: %f Dev: %f\n", ScalingCoeff, CurrDev);
			//system("pause");
			//if (CurrDev > Dev && Dev == 0){
			//	//ScalingCoeff -= ScalingStep;
			//	return GetBestCoeff(ScalingResults);
			//	DecFound = true;
			//}

			Dev = CurrDev;
		}
		Count++;
	}
	return ScalingCoeff;
}

void GenKron(const TStr& Args, TKronMtx& FitMtx, TFltPrV& KronDegAvgIn, TFltPrV& KronDegAvgOut, vector<Diap>& SmoothedDiaps, vector<int>& Prev){
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

	// Kronecker  graph
	PNGraph Kron;
	TIntPrV SamplesIn, SamplesOut;
	double Sec = 0.0;
	int AvgMaxOutDeg = 0, AvgMaxInDeg = 0, MinMaxOutDeg = 0, MaxMaxOutDeg = 0, MinMaxInDeg = 0, MaxMaxInDeg = 0;

	for (int i = 0; i < NKron; i++){
		ExecTime.Tick();
		KroneckerGen(Kron, FitMtx, NIter, IsDir, InDegR, OutDegR, NoiseCoeff);
		Rewire(Kron, SmoothedDiaps, OutDegR, Prev);
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
		ScaleFitMtx(ModelNodes, ModelEdges, MDegIn[MDegIn.Len()-1].Val1, MaxModelOutDeg, FitMtxM, NIter, IsDir, ScaleMtx);
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
		vector<Diap> SmoothedDiaps; vector<int> Prev;

		GenKron(Parameters[KRONGEN] + KronParameters, FitMtxM, KronDegAvgIn, KronDegAvgOut, SmoothedDiaps, Prev);

		PlotDegrees(Parameters, KronDegAvgIn, KronDegAvgOut, "kron");
		
		PNGraph  K;
		GetGraphFromAvgDistr(KronDegAvgOut, K);
		PlotMetrics(Parameters, K, "kron", TFile);
	}
}

// get FitMtx and scaling coefficient from small model
void GetFitMtxFromMS(TKronMtx& FitMtxM, TFlt& ScalingCoeff, vector<Diap>& SmoothedDiaps, const vector<TStr>& Parameters, vector<int>& Prev){
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
	ScaleFitMtx(ModelNodes, ModelEdges, MDegIn[MDegIn.Len()-1].Val1, MaxModelOutDeg, FitMtxM, NIter, IsDir, ScaleMtx);
	if (ScaleMtx == "true"){
		TKronMtx NewMtx(FitMtxM);
		NewMtx.Dump(TFile);
		TExeTm T;
		T.Tick();
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
		GetRelativeDiff(MDegOut, KronDeg, RelDiffNonCum);
		PrintDegDistr(KronDeg, "KronAvg.tab");
		//PrintDegDistr(MDegOut, "Model.tab");
		PrintRelDiff(RelDiffNonCum, "RelDiff.tab");
		GetSmoothedDiaps(RelDiffNonCum, SmoothedDiaps, Prev);
		PrintSmoothedDiaps(SmoothedDiaps, "SmoothedDiaps.tab");
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

		for (size_t i = 0; i < KronDeg.Len(); i++){
			int Deg = KronDeg[i].Val1, Count = KronDeg[i].Val2;
			if (Deg > DiapEnd) break;
			if (Deg < DiapBegin) continue;
			int SubBIndex = NewDiap.GetSubBIndex(Deg);
			if (SubBIndex < 0 || SubBIndex > DiapIt->second.Len() - 1)
				Error("GetDiaps", "Incorrect index of subdiapason");
			// add nodes with account of relative difference
			double RelDiff = DiapIt->second[SubBIndex];
			if (abs(RelDiff) != 1) NodesToAdd += RelDiff * Count;
			else if (RelDiff == 1.0) NodesToAdd += 1;
			else NodesToAdd -= Count;
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

	for (auto NodeIt = Kron->BegNI(); NodeIt != Kron->EndNI(); NodeIt++){
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

		// if nodes are deleted from current diapason
		// and they are clustered to form a new value
		if (IsDPlus == false && N.GetL() > Deg){
			
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

		int ReqDeg = N.GetRndDeg(Rnd);
		//printf("Init degree: %d\n", Deg);
		int ToDel = Deg - ReqDeg;

		// DEBUG

		if (ToDel <= 0)
			Error("Rewire", "ToDel is non-positive");

		// END DEBUG

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

			// DEBUG

			int Edges = Kron->GetEdges();
			if (Edges != BasicEdgesCount + Add - Del)
				Error("Rewire", "Basic edges count != Edges + Add - Del");

			// DEBUG END
		}

		D.DecreaseStratN();

	}
	
	cout <<  "Edges added " << Add << ", edges deleted " << Del << ", difference " << Add - Del << endl;

	
}

// ModelEdges - debug parameter
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
void Rewire(PNGraph& Kron, const vector<Diap>& SmoothedDiaps, const TIntPr& OutDegR, vector<int>& Prev){
	TRnd Rnd;
	int ModelEdges = Kron->GetEdges();
	TFltPrV KronDeg;
	TSnap::GetInDegCnt(Kron, KronDeg);
	PrintDegDistr(KronDeg, "Kron.tab");
	KronDeg.Sort();
	const TInt& DegMin = OutDegR.Val1, &DegMax = OutDegR.Val2;
	vector<Diaps> DPlus, DMinus;
	GetDiaps(DPlus, DMinus, SmoothedDiaps, KronDeg, DegMin, DegMax, Prev);
	GetRewireStrategies(DPlus, DMinus);
	PrintDiapsInfo(DPlus, "DPlus.tab");
	PrintDiapsInfo(DMinus, "DMinus.tab");
	Rewire(Kron, DPlus, DMinus, DegMin.Val, DegMax.Val);
	
	int Diff = Kron->GetEdges() - ModelEdges;
	cout << "Difference of edges: " << Kron->GetEdges() - ModelEdges << endl;
	AddEdges(Kron, Diff, DegMin, DegMax, ModelEdges, DMinus, DPlus);
	
	
}

// get smoothed diapasons for scaling
void GetSmoothedDiaps(const TFltPrV& RelDiffNonCum, vector<Diap>& SmoothedDiaps, vector<int>& Prev){
	if (RelDiffNonCum.Len() == 0)
		Error("GetSmoothedDiaps", "Array size = 0");
	const int DegCount = RelDiffNonCum.Len(), 
		DegMin = RelDiffNonCum[0].Val1, 
		DegMax = RelDiffNonCum[DegCount-1].Val1, 
		DiffDegs = DegMax - DegMin + 1;
	TInt ToleranceVal = 1;
	TFltV DiapDevV;
	bool DiapSign = RelDiffNonCum[0].Val2 > 0 ? true : false;
	int DiapBegin = 0, DiapEnd = 0; 
	int DiapIndex = 0, PrevDiapEnd = 0, PrevDeg = 0;

	for (size_t i = 0; i < DegCount; i++){
		int Deg = RelDiffNonCum[i].Val1;
		double Diff = RelDiffNonCum[i].Val2;
		bool Sign = Diff > 0 ? true : false;
		// if there is no difference, diapason should be finalized
		if (Diff == 0.0) Sign = DiapSign == true ? false : true;

		// if some degrees inside the diapasone are absent, add zero values of relative difference
		if (i != 0 && Sign == DiapSign && Deg != PrevDeg + 1){
			printf("%d %d\n", PrevDeg, Deg);
			for (size_t i = 0; i < Deg-PrevDeg-1; i++)
				DiapDevV.Add(0);
			if (i == DegCount - 1){
				DiapEnd = i;
				DiapDevV.Add(Diff);
			}
		}

		// if sign was changed or it is last interval
		if (Sign != DiapSign || i == DegCount - 1){
			int DegBegin = RelDiffNonCum[DiapBegin].Val1, 
				DegEnd = RelDiffNonCum[DiapEnd].Val1;
			printf("DegBegin: %d DegEnd: %d\n", DegBegin, DegEnd);
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
				DiapSign = Sign;
				DiapIndex++;
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
	TInt MDegCount = MDeg.Len(), KronDegCount = KronDeg.Len();
	TInt MinDegModel = static_cast<int>(MDeg[0].Val1), MaxDegModel = static_cast<int>(MDeg[MDegCount-1].Val1),
		MinDegKron = static_cast<int>(KronDeg[0].Val1), MaxDegKron = static_cast<int>(KronDeg[KronDegCount-1].Val1);
	TInt MinDeg = MinDegModel < MinDegKron ? MinDegModel : MinDegKron,
		MaxDeg = MaxDegModel > MaxDegKron ? MaxDegModel : MaxDegKron;
	TFlt CurrDeg = MinDeg, MInd = 0, KronInd = 0;
	if (NonCum){
		while (1){
			TFlt MDegVal, KronDegVal, MDegCount, KronDegCount;
			if (MInd < MDeg.Len()){
				MDegVal = MDeg[MInd].Val1; MDegCount = MDeg[MInd].Val2;
			}
			else {
				MDegVal = MaxDeg + 1; MDegCount = 0;
			}
			if (KronInd < KronDeg.Len()) {
				KronDegVal = KronDeg[KronInd].Val1; KronDegCount = KronDeg[KronInd].Val2;
			}
			else {
				KronDegVal = MaxDeg + 1; KronDegCount = 0;
			}
			bool MLessDeg = MDegVal < KronDegVal ? true : false;
			CurrDeg = MLessDeg ? MDegVal : KronDegVal;
			if (CurrDeg > MaxDeg) break;
			TFlt RelDiff;
			if (MDegVal == CurrDeg && KronDegVal == CurrDeg){
				RelDiff = (MDegCount - KronDegCount) / KronDegCount;
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
			RelDiffV.Add(RelDiffPr);
		}
	}
}


// get average estimates of out-degree of NKron Kronecker graphs
void GetAvgKronDeg(const TKronMtx& NewMtx, const TInt& NIter, const TStr& IsDir, const TInt& NKron, const TIntPr& ModelDegR, TFltPrV& KronDeg){
	TIntPrV Samples;
	for (int i = 0; i < NKron; i++){
		PNGraph G;
		KroneckerGen(G, NewMtx, NIter, IsDir, ModelDegR, ModelDegR, 0);
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

	CheckParams(ModelGen, ModelPlt);
	CheckParams(MSGen, MSPlt);

	TFile = OpenFile(StatFile.CStr());

	PyInit("PySettings.txt");

	if (FromMS == "yes"){
		TKronMtx FitMtx;
		TFlt ScalingCoeff = 0;
		vector <TStr> Parameters;
		GetParameters(CommandLineArgs, "Small", Parameters);
		vector<Diap> SmoothedDiaps; vector<int> Prev;
		GetFitMtxFromMS(FitMtx, ScalingCoeff, SmoothedDiaps, Parameters, Prev);
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
		ScaleFitMtx(G->GetNodes(), G->GetEdges(), RequiredDeg, RequiredDeg, FitMtx, NIter, IsDir, "true");
		GenKron(Parameters[KRONGEN] + KronParameters, FitMtx, KronDegAvgIn, KronDegAvgOut, SmoothedDiaps, Prev);
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

