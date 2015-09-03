#include "stdafx.h"

int GetNodesCountInDiap(const TFltPrV& Deg, int L, int R){
	int N = 0;
	for (size_t i = 0; i < Deg.Len(); i++){
		if (Deg[i].Val1 < L)
			continue;
		if (Deg[i].Val1 > R)
			break;
		N += Deg[i].Val2;
	}
	return N;
}

// add to Ratio (DiapBegin, DiapEnd, log(MNodes+1)/log(KNodes+1))
void GetDiapsNodesKMRatio(const TFltPrV& MDeg, const TFltPrV& KDeg, Diaps& D, IntIntDblV& Ratio){
	int L = D.GetL(), R = D.GetR();
	int MN = GetNodesCountInDiap(MDeg, L, R), 
		KN = GetNodesCountInDiap(KDeg, L, R);
	double Val = (KN != 0) ? log10(MN + 1.00) / log10(KN + 1.00) : 0;
	Ratio.push_back(make_pair(make_pair(L, R), Val));
}

// add to Ratio (CurrBegin, CurrEnd, log(PrevNodes+1)/log(CurrNodesNodes+1))
void GetDiapNodesMMRatio(const TFltPrV& MDeg, Diaps& Prev, Diaps& Curr, IntIntDblV& Ratio){
	int LPrev = Prev.GetL(), RPrev = Prev.GetR(),
		L = Curr.GetL(), R = Curr.GetR();
	int PrevN = GetNodesCountInDiap(MDeg, LPrev, RPrev), 
		N = GetNodesCountInDiap(MDeg, L, R);
	double Val = (N != 0) ? log10(PrevN + 1.00) / log10(N + 1.00) : 0;
	Ratio.push_back(make_pair(make_pair(L, R), Val));
}

void TestModelKronRelation(const PNGraph& M, const PNGraph& K){
	TFltPrV MDeg, KDeg, RD;
	TSnap::GetInDegCnt(M, MDeg);
	TSnap::GetInDegCnt(K, KDeg);
	GetRelativeDiff(MDeg, KDeg, RD);
	vector<BaseDiap> SDiaps;
	vector<int> Prev;
	GetBaseDiaps(MDeg, KDeg, SDiaps);
	vector<Diaps> DPlus, DMinus;
	int DegMin = static_cast<int>(MDeg[0].Val1),
		DegMax = static_cast<int>(MDeg[MDeg.Len()-1].Val1);
	GetDiaps(DPlus, DMinus, SDiaps, KDeg, DegMin, DegMax);
	int DPlusInd = 0, DMinusInd = 0, DPlusLDeg = 0, DMinusLDeg = 0;
	int DiapInd = 0;
	// true - DPlus, false - DMinus, int - index in DPlus/DMinus
	pair<bool, int> PrevDiap;
	bool DPlusFin = false, DMinusFin = false;
	IntIntDblV MKRatio, MMRatio;
	while (!DPlusFin || !DMinusFin){
		if (!DPlusFin)
			DPlusLDeg = DPlus[DPlusInd].GetL();
		if (!DMinusFin)
			DMinusLDeg = DMinus[DMinusInd].GetL();

		if (DPlusFin || (!DMinusFin && DMinusLDeg < DPlusLDeg)){
			if (DiapInd != 0){
				if (PrevDiap.first)
					GetDiapNodesMMRatio(MDeg, DPlus[PrevDiap.second], DMinus[DMinusInd], MMRatio);
				else
					GetDiapNodesMMRatio(MDeg, DMinus[PrevDiap.second], DMinus[DMinusInd], MMRatio);
			}
			PrevDiap.first = false;
			PrevDiap.second = DMinusInd;
			GetDiapsNodesKMRatio(MDeg, KDeg, DMinus[DMinusInd++], MKRatio);
			if (DMinusInd == DMinus.size()) DMinusFin = true;
		}
		else if (DMinusFin || (!DPlusFin && DPlusLDeg < DMinusLDeg)){
			if (DiapInd != 0){
				if (PrevDiap.first)
					GetDiapNodesMMRatio(MDeg, DPlus[PrevDiap.second], DPlus[DPlusInd], MMRatio);
				else
					GetDiapNodesMMRatio(MDeg, DMinus[PrevDiap.second], DPlus[DPlusInd], MMRatio);
			}
			PrevDiap.first = true;
			PrevDiap.second = DPlusInd;
			GetDiapsNodesKMRatio(MDeg, KDeg, DPlus[DPlusInd++], MKRatio);
			if (DPlusInd == DPlus.size()) DPlusFin = true;
		}
		DiapInd++;
	}
	PrintIntIntDblV(MKRatio, "MKRatio.tab");
	PrintIntIntDblV(MMRatio, "MMRatio.tab");
}