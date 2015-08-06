#include "stdafx.h"
#include "Error.h"

void DefaultSettings(vector<TStr>& out){
	out.clear();
	printf("Applying default settings...\n");
	TStr s("1 ");
	for (int i = 0; i < NFUNC; i++){
		out.push_back(s);
	}
}

double PrintLargestEigenVal(const PNGraph& G, ofstream& F, const TStr& GName){
	TFltV EigenVals;
	TSnap::GetEigVals(TSnap::ConvertGraph<PUNGraph>(G), 1, EigenVals);
	double MaxVal = EigenVals[0].Val;
	F << GName.CStr() << " largest eigenvalue: " << EigenVals[0].Val << endl; 
	return MaxVal;
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

void CheckParams(const TStr& model_gen, const TStr& model_plt)
{
	try
	{
		if (model_gen != model_plt)
		{
			if (model_gen != "model+kron")
				throw 1;
			if (model_gen == "none" && model_plt != "none")
				throw 1;
			if (model_plt == "model+kron" && model_gen != "model+kron")
				throw 1;
		}
	}
	catch (int i){
		printf("Inconsistency in KRONTEST parameters\n");
		system("pause");
		exit(1);
	}
}


void GetMtxFromSepLine(const TStr& line, const TStr& separator, TFltV& matrix){
	TStrV strVals;
	line.SplitOnAllAnyCh(separator, strVals);
	for (int i = 0; i < strVals.Len(); i++) 
		matrix.Add(strVals[i].GetFlt());
}

void ReadPNGraphFromFile(const TStr args, PNGraph& G){
	Try
	Env = TEnv(args, TNotify::StdNotify);
	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../as20graph.txt", "Input graph file (single directed edge per line)");
	// load graph
	G = TSnap::LoadEdgeList<PNGraph>(InFNm, 0, 1);
	Catch
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

void PrintDegDistr(const PNGraph& G, const TStr& OutFNm){
	TFltPrV InDeg, OutDeg;
	TSnap::GetInDegCnt(G, InDeg);
	TSnap::GetInDegCnt(G, OutDeg);
	PrintDegDistr(InDeg, OutFNm + "In.tab");
	PrintDegDistr(OutDeg, OutFNm + "Out.tab");
}

void ReadMtx(const TStr& Mtx, const TInt& MtxSize, TKronMtx& FitMtx){
	TFltV matrix;
	GetMtxFromSepLine(Mtx, ";", matrix);
	FitMtx.GenMtx(matrix.Len() / MtxSize);
	FitMtx.SetMtx(matrix);
}

void MakeDatFile(const TStr& Name, const TStr& AddStr, const TStrV& ColumnNames, const vector<vector<double>>& Data, const int& Nodes, const int& Edges){
	if (Data.size() == 0)
		Error("MakeDatFile", "Data array is empty");
	int ColumnsCount = ColumnNames.Len();
	if (Data.size() != ColumnsCount)
		Error("MakeDatFile", "Data array is inconsistent with columns count");
	for (int i = 0; i < Data.size()-1; i++){
		if (Data[i].size() != Data[i+1].size())
			Error("MakeDatFile", "The sizes of data arrays don't match together");
	}
	
	time_t ltime;  time(&ltime);
    char* TimeStr = ctime(&ltime);  TimeStr[strlen(TimeStr) - 1] = 0;
	TStr FileName = Name;
	FileName.InsStr(FileName.Len(), ".tab");
	ofstream F(FileName.CStr());
	F << "#\n";
	F << Name.CStr() << ". ";
	if (Nodes != 0)
		F << "G(" << Nodes << "," << Edges << "). ";
	F << AddStr.CStr() << " (" << TimeStr << ") ";
	F << "\n#\n# ";

	for (int i = 0; i < ColumnsCount; i++)
		F << ColumnNames[i].CStr() << " ";
	F << "\n";
	int DataSize = Data[0].size();
	for (int i = 0; i < DataSize; i++){
		for (int j = 0; j < ColumnsCount; j++){
			F << Data[j][i] << "\t";
		}
		F << "\n";
	}

	F.close();
}

int GetDigits(int i){
	return ceil(log(static_cast<double>(i)))/log(2.00) + 1;
}

int GetWeight(int i){
	if (i == 0) return 0;
	else {
		double Weight = 0;
		for (int j = 0; j < GetDigits(i); j++)
			Weight += (i >> j) & 1;
		return Weight;
	}
}

// i is the count of (i=ival, j=jval) bits in i,j
int GetPow(int i, int j, int ival, int jval){
	int DigitsI = GetDigits(i), DigitsJ = GetDigits(j);
	int Digits = DigitsI > DigitsJ ? DigitsI : DigitsJ;
	int Res = 0;
	for (int k = 0; k < Digits; k++){
		int IBit = (i >> k) & 1, JBit = (j >> k) & 1;
		if (IBit == ival && JBit == jval) Res++; 
	}
	return Res;
}

double GetProb(const int &i, const int&j, const TKronMtx& FitMtx, const int& NIter){
	double Prob = 1.0;
	for (int k = 0; k < NIter; k++){
		int IBit = (i >> k) & 1, JBit = (j >> k) & 1;
		Prob *= FitMtx.At(IBit, JBit);
	}
	return Prob;
}

void PrintNodeDegrees(const PNGraph& G, const TKronMtx& FitMtx, const int& NIter){
	vector<vector<double>> Data;
	vector<double> NodeId; vector<double> Deg; vector<double> ExpectedDeg; vector<double> SumDeg;
	double A = FitMtx.At(0,0), B = FitMtx.At(0,1), C = FitMtx.At(1,0), D = FitMtx.At(1,1);
	for (int i = 0; i < G->GetNodes(); i++){
		NodeId.push_back(i);
		// get weight of node i (count of 1's in its binary representation)
		int Weight = GetWeight(i);
		Deg.push_back(G->GetNI(i).GetOutDeg());
		ExpectedDeg.push_back(pow(A+B, NIter-Weight)*pow(C+D, Weight));
		double SDeg = 0.0;
		for (int j = 0; j < G->GetNodes(); j++){
			double Prob = GetProb(i, j, FitMtx, NIter);
			//printf("%.2f ", Prob);
			SDeg += Prob;
			//int PowI = GetPow(i, j, 0, 1), PowJ = GetPow(i, j, 1, 1);
			//SDeg += pow(A, NIter - Weight - PowI) * pow(B, PowI) * pow(C, Weight - PowJ) * pow(D, PowJ);
		}
		//printf("\n");
		SumDeg.push_back(SDeg);		
	}
	Data.push_back(NodeId); Data.push_back(Deg); Data.push_back(ExpectedDeg); Data.push_back(SumDeg);
	TStrV ColumnNames; ColumnNames.Add("NodeId");ColumnNames.Add("Deg");ColumnNames.Add("E_N[Deg]");ColumnNames.Add("E_N[Deg](Sum)");
	MakeDatFile("DegTest", "Node degrees", ColumnNames, Data, G->GetNodes(), G->GetEdges());
}

void GetParameters(const vector<TStr>& CommandLineArgs, const TStr& Type, vector<TStr>& Parameters){
	Env = TEnv(CommandLineArgs[KRONTEST], TNotify::NullNotify);
	// type of plots
	const TStr Plt = Env.GetIfArgPrefixStr("-plt:", "all", "Type of plots (cum, noncum, all)");
	// full - all points of distrib will be plotted; expbin - exponential binning
	const TStr PType = Env.GetIfArgPrefixStr("-ptype:", "all", "How to plot (full, expbin, all)");
	// radix of binning
	const TInt BinRadix = Env.GetIfArgPrefixInt("-bin:", 2, "Radix for exponential binning");
	// number of eigenvalues to investigate
	const TStr NEigen = Env.GetIfArgPrefixStr("-neigen:", "10", "Number of eigenvalues");
	// plot hops (plot/none)
	const TStr Hops = Env.GetIfArgPrefixStr("-hops:", "none", "Plot hops (plot/none)");
	// clustering coefficient (yes - calculate general CC, yes+plot - caculate general CC + plot local CC, none)
	const TStr Clust = Env.GetIfArgPrefixStr("-clust:", "none", "Calculate/plot clustering coefficient (yes/yes+plot/none)");

	int BegIndex = 1;
	if (Type == "Small") BegIndex += NPARCOPY;
	for (size_t i = BegIndex; i < BegIndex + NPARCOPY; i++)
		Parameters.push_back(CommandLineArgs[i]);
	Parameters.push_back(PType); Parameters.push_back(Plt); Parameters.push_back(Type); Parameters.push_back(NEigen); 
	Parameters.push_back(BinRadix.GetStr()); Parameters.push_back(Hops); Parameters.push_back(Clust);
}