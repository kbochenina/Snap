#include "stdafx.h"

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