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
	// weight of diapason (NDiap/NCount)
	double Weight;
	// borders of diapason
	pair<int, int> Borders;
	// subborders according to the relation of BaseLen to Len
	vector<pair<int, int>> SubB;
	// cumulative probabilities to add nodes to subdiapasons
	vector<double> Prob;
	// part of nodes of diapason in subdiapasons
	vector<double> NParts;
	// test
	void TestSubB();
public:
	BaseDiap(int I, pair<int, int> B, int BL, double MK, double Prev, double W);
	int Length() const {return Len;}
	int BaseLength() const {return BaseLen;}
	int GetSubBIndex (int Deg);
	int GetL() const {return Borders.first;}
	int GetR() const {return Borders.second;}
	double GetWeight() const {return Weight;}
	int GetIndex() {return Index;}
	void GetProb(vector<double>& P) const;
	void GetNParts(vector<double>& NP) const;
	double GetMKRatio() const {return MKRatio;}
	double GetPrevCurrKRatio() const {return PrevCurrKRatio;}
	// set subborders
	void SetSubB(vector<pair<int,int>> SB, vector<double> P, vector<double> NP);
	// check if degree belongs to diapasone
	bool IsDegInDiap(int Deg) { if (Deg >= Borders.first && Deg <= Borders.second) return true; return false; }
	// print diapason info
	virtual void PrintInfo(ofstream& F);
	~BaseDiap(void);
};

