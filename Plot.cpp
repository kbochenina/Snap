#include "stdafx.h"

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
	return false;
}
// CHECK
void ExpBinning(const TFltPrV& deg, TFltPrV& degSparse, const int& BinRadix){
	TFlt maxDeg(deg[deg.Len()-1].Val1.Val), minDeg(deg[0].Val1.Val);
	bool maxPowerReached = false;
	// idx - index of border, previdx - index of previous border
	int power = 0, previdx = 0, idx, binSize;
	TFltPr val;
	double binBorder = 0.0;
	while (binBorder <= minDeg)
		binBorder = pow(static_cast<double>(BinRadix), power++);

	TFltPr v(minDeg, deg[0].Val2.Val);
	degSparse.Add(v);

	bool isExact = false;
	while (!maxPowerReached){
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
		// if prevBinBorder was the smallest degree, it can be more than binBorder / BinRadix
		double SumBinBorder = previdx > 0 ? binBorder + static_cast<double>(binBorder) / BinRadix : binBorder + static_cast<double>(minDeg); 
		double avgDeg = SumBinBorder / 2.0;
		val.Val1 = avgDeg; val.Val2 = sum;
		degSparse.Add(val);
		previdx = idx;
		binBorder = pow(static_cast<double>(BinRadix), power++);
	}
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

void SaveDegree(const TFltPrV& deg, const TStr& n, bool isIn, bool isCum, bool calcCum){
	TFltPrV d(deg);
	d.Sort();
	int nodes, edges;
	GetNodesEdgesCountFromDegDistr(d, nodes, edges);
	TSnap::PlotDegDistr(d, nodes, edges, n, n, isCum, false, isIn, calcCum);
}

// plot all points without binning
void PlotPoints(const TFltPrV& in, const TFltPrV& out, const TStr& name, const TStr& Plt){
		if (Plt == "cum" || Plt == "all"){
			SaveDegree(in, "kron" + name, true, true);
			SaveDegree(out, "kron" + name, false, true);
		}
		if (Plt == "noncum" || Plt == "all"){
			SaveDegree(in, "kron" + name, true, false);
			SaveDegree(out, "kron" + name, false, false);
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

void PlotSparse(const vector<TFltPrV>& distr, const TStrV& names, bool isIn, const TStr& Plt, const TInt& BinRadix){
	for (size_t i = 0; i < distr.size(); i++){
		if (Plt == "cum" || Plt == "all"){
			SaveSparse(distr[i], BinRadix, isIn, names[i], true);
		}
		if (Plt == "noncum" || Plt == "all")
			SaveSparse(distr[i], BinRadix, isIn, names[i], false);
	}
}