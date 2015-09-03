#include "stdafx.h"


void ScaleFitMtx(int ModelNodes, int ModelEdges, TFlt MaxDegInModel, TFlt MaxDegOutModel, TKronMtx& FitMtx, const TInt& NIter, TStr IsDir, const TStr& ScaleMtx, ofstream& TFile){
	// check if ModelNodes and NIter are consistent
	if (pow(FitMtx.GetDim(), (double)NIter) != ModelNodes)
		Error("ScaleFitMtx", "ModelNodes is not equal to 2^NIter");
	// scale init matrix to match number of edges
	ScaleFitMtxForEdges(FitMtx, NIter, ModelEdges, TFile);	

	// equalizing B and C for undirected graphs
	if (IsDir == "false" && FitMtx.At(0,1) != FitMtx.At(1,0)){
		cout << "WARNING. B and C is not equal for an undirected graph. B and C will be eqialized" << endl;
		ScaleFitMtxForUnDir(FitMtx, TFile);
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
			ScaleFitMtxForUnDir(FitMtx, TFile);
		}
		ScaleFitMtx(FitMtx, NIter, ModelNodes, MaxDegOutModel, IsDir);
		FitMtx.Dump();
		TFile << "Expected maximum degree in Kronecker graph: in " << FitMtx.GetMaxExpectedDeg(NIter, IsDir, true) << 
			",  out " << FitMtx.GetMaxExpectedDeg(NIter, IsDir, false) << endl;
	}
}


void ScaleFitMtxForEdges(TKronMtx& FitMtx, const TInt& NIter, const int& ExpectedModelEdges, ofstream& TFile){
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

void ScaleFitMtxForUnDir(TKronMtx& FitMtx, ofstream& TFile){
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
	return 0.3;
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