#include "stdafx.h"

void ExpBinning(const TFltPrV& deg, TFltPrV& degSparse, const int& BinRadix){
	TFlt maxDeg(deg[deg.Len()-1].Val1.Val), minDeg(deg[0].Val1.Val);
	bool maxPowerReached = false;
	// idx - index of border, previdx - index of previous border
	int power = 0, previdx = 0, idx = 0, binSize = 0;
	double binBorder = 0.0;
	while (binBorder <= minDeg)
		binBorder = pow(static_cast<double>(BinRadix), power++);
	bool isExact = false;
	while (!maxPowerReached){
		if (power == 1){
			// if there are nodes with degree 1
			idx = FindVal1Elem(deg, 1, isExact);
			if (isExact){
				TFltPr val(1, deg[idx].Val2.Val);
				degSparse.Add(val);
				previdx = idx;
			}
		}
		else {
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
			double avgDeg = (binBorder + static_cast<double>(binBorder) / BinRadix) / 2.0;
			TFltPr val(avgDeg, sum);
			degSparse.Add(val);
			previdx = idx;
		}
		binBorder = pow(static_cast<double>(BinRadix), power++);
	}
}
