#include "StdAfx.h"
#include "BaseDiap.h"


BaseDiap::BaseDiap(int I, pair<int, int> B, int BL, double MK, double Prev) {
	if (I < 0)
		Error("BaseDiap::Diaps", "I < 0");
	if (MK < 0)
		Error("BaseDiap::Diaps", "MKRatio < 0");
	if (Prev < 0)
		Error("BaseDiap::Diaps", "PrevCurrKRatio < 0");
	if (B.first < 0 || B.second < 0)
		Error("BaseDiap::Diaps", "Negative border");
	if (B.first > B.second)
		Error("BaseDiap::Diaps", "B.first > B.second");
	if (BL > B.second - B.first + 1)
		Error("BaseDiap::Diaps", "Base length is more that real length");
	Index = I;
	Borders.first = B.first;
	Borders.second = B.second;
	BaseLen = BL;
	Len = Borders.second - Borders.first + 1;
	MKRatio = MK;
	PrevCurrKRatio = Prev;
}

// set subborders
void BaseDiap::SetSubB(vector<pair<int,int>> SB, vector<double> P, vector<double> NP){
	if (SB.size() != BaseLen)
		Error("BaseDiap::SetSubB", "SB.size() != BaseLen");
	if (P.size() != BaseLen)
		Error("BaseDiap::SetSubB", "P.size() != BaseLen");
	if (NP.size() != BaseLen)
		Error("BaseDiap::SetSubB", "NP.size() != BaseLen");
	for (size_t i = 0; i < BaseLen; i++){
		SubB.push_back(make_pair(SB[i].first, SB[i].second));
		Prob.push_back(P[i]);
		NParts.push_back(NP[i]);
	}
	TestSubB();
}

// test
void BaseDiap::TestSubB(){
	int Count = 0;
	int &L = Borders.first, &R = Borders.second;
	if (SubB.size() != BaseLen)
		Error("BaseDiap::TestSubB", "Subborders count != BaseLen");
	for (auto it = SubB.begin(); it != SubB.end(); it++){
		if (it - SubB.begin() == 0 && it->first != L)
			Error("BaseDiap::TestSubB", "Left border of 1st subdiapason != L");
		if (it - SubB.begin() == SubB.size() && it->second != R)
			Error("BaseDiap::TestSubB", "Left border of 1st subdiapason != L");
		if (it - SubB.begin() != 0 && it->first != (it-1)->second + 1)
			Error("BaseDiap::TestSubB", "Left border of current subdiapason != right border of previous diapason + 1");
		if (it->first > it->second)
			Error("BaseDiap::TestSubB", "B.first > B.second");
		Count += it->second - it->first + 1;
	}
	if (Count != Len)
		Error("BaseDiap::TestSubB", "Count != Len");
}

int BaseDiap::GetSubBIndex(int Deg){
	if (SubB.size() == 0)
		Error("BaseDiap::GetSubBIndex", "SubB size == 0");
	if (Deg < SubB[0].first || Deg > SubB[SubB.size()-1].second)
		Error("BaseDiap::GetSubBIndex", "Deg is out of range");
	for (int i = 0; i < SubB.size(); i++){
		if (Deg >= SubB[i].first && Deg <= SubB[i].second)
			return i;
	}
	Error("BaseDiap::GetSubBIndex", "Deg not found");
}

void BaseDiap::GetProb(vector<double>& P) const{
	P.clear();
	if (Prob.size() != BaseLen)
		Error("BaseDiap::GetProb", "Prob is not set");
	for (int i = 0; i < BaseLen; i++)
		P.push_back(Prob[i]);
}

void BaseDiap::GetNParts(vector<double>& NP) const{
	NP.clear();
	if (NParts.size() != BaseLen)
		Error("BaseDiap::GetProb", "Prob is not set");
	for (int i = 0; i < BaseLen; i++)
		NP.push_back(NParts[i]);
}

// print node info
void BaseDiap::PrintInfo(ofstream& F){
	F << "Index: " << Index << "[" << Borders.first << ";" << Borders.second << "]" <<  endl;
	F << "MKRatio: " << MKRatio << endl;
	F << "PrevCurrKRatio: " << PrevCurrKRatio << endl;
	F << "Subborders: " << endl;
	for (size_t i = 0; i < SubB.size(); i++) 
		F << "[" << SubB[i].first << ";" << SubB[i].second << "]" << " Prob: " << Prob[i] << " Np: " << NParts[i] << endl;
	F << endl;
}

BaseDiap::~BaseDiap(void)
{
}
