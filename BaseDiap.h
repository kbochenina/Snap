#pragma once

// MN - Model Nodes
// KN - Kron Nodes

class BaseDiap
{
protected:
	// index of diapason
	int Index;
	// base length of diapason
	int BaseLen;
	// real length of diapason
	int Len;
	// log10(MN + 1)/log10(KN + 1) for the whole diapason
	double MKRatio;
	// log10(PrevKN + 1)/log10(CurrKN+1)
	double PrevCurrKRatio;
	// borders of diapason
	pair<int, int> Borders;
	// subborders according to the relation of BaseLen to Len
	vector<pair<int, int>> SubB;
	// probabilities of subdiapasons
	vector<double> Prob;
	// test
	void TestSubB();
public:
	BaseDiap(int I, pair<int, int> B, int BL, double MK, double Prev);
	int Length() {return Len;}
	int GetSubBIndex(int Deg);
	int GetL(){return Borders.first;}
	int GetR(){return Borders.second;}
	int GetIndex() {return Index;}
	double GetMKRatio() {return MKRatio;}
	double GetPrevCurrKRatio() {return PrevCurrKRatio;}
	// set subborders
	void SetSubB(vector<pair<int,int>> SB, vector<double> P);
	// check if degree belongs to diapasone
	bool IsDegInDiap(int Deg) { if (Deg >= Borders.first && Deg <= Borders.second) return true; return false; }
	// print diapason info
	virtual void PrintInfo(ofstream& F);
	~BaseDiap(void);
};

