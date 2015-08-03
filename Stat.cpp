#include "stdafx.h"

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



int GetMaxDeg(const PNGraph& G, const TStr& IsDir, const TStr& IsIn){
	TIntPrV DegCnt;
	if (IsDir == "false"){
		PUNGraph U = TSnap::ConvertGraph<PUNGraph>(G);
		if (IsIn == "true")
			TSnap::GetInDegCnt(U, DegCnt);
		else 
			TSnap::GetOutDegCnt(U, DegCnt);
	}
	else{
		if (IsIn == "true")
			TSnap::GetInDegCnt(G, DegCnt);
		else
			TSnap::GetOutDegCnt(G, DegCnt);
	}
	// sort in descending order
	DegCnt.Sort(false);
	return DegCnt[0].Val1;
}

