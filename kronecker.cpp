#include "stdafx.h"
#include "kronecker.h"
#include "Error.h"
#include "Eigen.h"

void PrintDeg(const TFltPrV& distr, const TStr& OutFNm){
	FILE *F = stdout;
	if (! OutFNm.Empty()) F = fopen(OutFNm.CStr(), "wt");
	fprintf(F, "\n");
	fprintf(F, "  Degree\tVal\n");
	for (int i = 0; i < distr.Len(); i++){
		fprintf(F, "%f\t%f\n", distr[i].Val1.Val, distr[i].Val2.Val);
	}
	if (! OutFNm.Empty()) { fclose(F); }
}

void PrintDeg(const PNGraph& G, const TStr& OutFNm){
	TFltPrV InDeg, OutDeg;
	TSnap::GetInDegCnt(G, InDeg);
	TSnap::GetOutDegCnt(G, OutDeg);
	PrintDeg(InDeg, OutFNm + "In.tab");
	PrintDeg(OutDeg, OutFNm + "Out.tab");
}


/////////////////////////////////////////////////
// Kronecker Graphs
const double TKronMtx::NInf = -DBL_MAX;
TRnd TKronMtx::Rnd = TRnd(0);

TKronMtx::TKronMtx(const TFltV& SeedMatrix) : SeedMtx(SeedMatrix) {
  MtxDim = (int) sqrt((double)SeedMatrix.Len());
  IAssert(MtxDim*MtxDim == SeedMtx.Len());
}

void TKronMtx::SaveTxt(const TStr& OutFNm) const {
  FILE *F = fopen(OutFNm.CStr(), "wt");
  for (int i = 0; i < GetDim(); i++) {
    for (int j = 0; j < GetDim(); j++) {
      if (j > 0) fprintf(F, "\t");
      fprintf(F, "%f", At(i,j)); }
    fprintf(F, "\n");
  }
  fclose(F);
}

TKronMtx& TKronMtx::operator = (const TKronMtx& Kronecker) {
  if (this != &Kronecker){
    MtxDim=Kronecker.MtxDim;
    SeedMtx=Kronecker.SeedMtx;
  }
  return *this;
}

bool TKronMtx::IsProbMtx() const {
  for (int i = 0; i < Len(); i++) {
    if (At(i) < 0.0 || At(i) > 1.0) return false;
  }
  return true;
}

void TKronMtx::SetRndMtx(const int& PrmMtxDim, const double& MinProb) {
  MtxDim = PrmMtxDim;
  SeedMtx.Gen(MtxDim*MtxDim);
  for (int p = 0; p < SeedMtx.Len(); p++) {
    do {
      SeedMtx[p] = TKronMtx::Rnd.GetUniDev();
    } while (SeedMtx[p] < MinProb);
  }
}

void TKronMtx::SetEpsMtx(const double& Eps1, const double& Eps0, const int& Eps1Val, const int& Eps0Val) {
  for (int i = 0; i < Len(); i++) {
    double& Val = At(i);
    if (Val == Eps1Val) Val = double(Eps1);
    else if (Val == Eps0Val) Val = double(Eps0);
  }
}

// scales parameter values to allow Edges
void TKronMtx::SetForEdges(const int& Nodes, const int& Edges) {
  const int KronIter = GetKronIter(Nodes);
  const double EZero = pow((double) Edges, 1.0/double(KronIter));
  const double Factor = EZero / GetMtxSum();
  for (int i = 0; i < Len(); i++) {
    At(i) *= Factor;
    if (At(i) > 1)
    {
	    At(i) = 1;
    }
  }
}

// scales parameter values to allow Edges - without cutting
void TKronMtx::SetForEdgesNoCut(const int& Nodes, const int& Edges) {
  const int KronIter = GetKronIter(Nodes);
  const double EZero = pow((double) Edges, 1.0/double(KronIter));
  const double Factor = EZero / GetMtxSum();
  for (int i = 0; i < Len(); i++) {
    At(i) *= Factor;
   }
}

double TKronMtx::GetMax() const {
	double Val = 0.0;
	if(At(0,0) > At(0,1))
		Val = At(0,0);
	else Val = At(0,1);
	if (Val < At(1,0))
		Val = At(1,0);
	if (Val < At(1,1))
		Val = At(1,1);
	return Val;
}

double TKronMtx::GetMaxExpectedDeg(){
	double A = At(0,0), B = At(0,1), C = At(1, 0), D = At(1, 1);
	double Sum0 = pow(A + B, 1)/2 + pow(A + C, 1)/2,
		Sum1 = pow(A + B, 1)/2 + pow(B + D, 1)/2,
		Sum2 = pow(C + D, 1)/2 + pow(C + A, 1)/2,
		Sum3 = pow(C + D, 1)/2 + pow(D + B, 1)/2;
		
	//printf("%f %f %f %f\n", Sum0, Sum1, Sum2, Sum3);
	double MaxSum = 0;
	if (Sum0 > Sum1 && Sum0 > Sum2 && Sum0 > Sum3) {MaxSum = Sum0; }
	else if (Sum1 > Sum0 && Sum1 > Sum2 && Sum1 > Sum3) {MaxSum = Sum1; }
	else if (Sum2 > Sum0 && Sum2 > Sum1 && Sum2 > Sum3) {MaxSum = Sum2; }
	else if (Sum3 > Sum0 && Sum3 > Sum1 && Sum3 > Sum2) {MaxSum = Sum3; }

	return MaxSum;
}

double TKronMtx::GetMaxExpectedDeg(const int& NIter)
{
	int BestRow, BestCol;
	return GetMaxExpectedDeg(At(0,0), At(0,1), At(1,0), At(1,1), NIter, BestRow, BestCol);
}

void TKronMtx::Normalize()
{
	double Sum = GetMtxSum();
	for (int i = 0; i < MtxDim; i++)
	{
		for (int j = 0; j < MtxDim; j++)
		{
			double &Val = At(i,j);
			if (Val > 1)
			{
				int I1 = i, J1 = j ? 0 : 1, I2 = i ? 0 :  1, J2 = j;
				double Diff = Val - 1;
				At(I1, J1) += Diff / 2;
				At(I2, J2) += Diff / 2;
				At(i, j) = 1;
			}
		}
	}
	//if (Sum != GetMtxSum())
	//	printf("Error while normalizing\n");
}


 double TKronMtx::GetMaxExpectedDeg(const double&A, const double&B, const double&C, const double&D, const int& NIter, int& BestRow, int&BestCol){
	double E = pow(A + B + C + D, NIter);
	double NewNIter = log10(E/2) / log10(A + B + C + D);

	double Possible = pow(A + B, NewNIter);
	double Possible1 = pow(A + B, NIter);

	//double Sum0 = pow(A + B, NewNIter)/2 + pow(A + C, NewNIter)/2,
	//	Sum1 = pow(A + B, NewNIter)/2 + pow(B + D, NewNIter)/2,
	//	Sum2 = pow(C + D, NewNIter)/2 + pow(C + A, NewNIter)/2,
	//	Sum3 = pow(C + D, NewNIter)/2 + pow(D + B, NewNIter)/2;

	double Sum0 = pow(A + B, NIter) + pow(A + C, NIter),
	Sum1 = pow(A + B, NIter) + pow(B + D, NIter),
	Sum2 = pow(C + D, NIter) + pow(C + A, NIter),
	Sum3 = pow(C + D, NIter) + pow(D + B, NIter);

	/*double Sum0 = pow(A + B, NIter)/2 + pow(A + C, NIter)/2,
		Sum1 = pow(A + B, NIter)/2 + pow(B + D, NIter)/2,
		Sum2 = pow(C + D, NIter)/2 + pow(C + A, NIter)/2,
		Sum3 = pow(C + D, NIter)/2 + pow(D + B, NIter)/2;*/

	//printf("%f %f %f %f\n", Sum0, Sum1, Sum2, Sum3);
	double MaxSum = 0;
	if (Sum0 > Sum1 && Sum0 > Sum2 && Sum0 > Sum3) {MaxSum = Sum0; BestRow = 0; BestCol = 0;}
	else if (Sum1 > Sum0 && Sum1 > Sum2 && Sum1 > Sum3) {MaxSum = Sum1; BestRow = 0; BestCol = 1;}
	else if (Sum2 > Sum0 && Sum2 > Sum1 && Sum2 > Sum3) {MaxSum = Sum2; BestRow = 1; BestCol = 0;}
	else if (Sum3 > Sum0 && Sum3 > Sum1 && Sum3 > Sum2) {MaxSum = Sum3; BestRow = 1; BestCol = 1;}

	return MaxSum;
 }

// works only for 2x2 size matrix!
void TKronMtx::SetForMaxDeg(const double& MaxDeg, const int& NIter)
{
	if (MaxDeg <= 0)
		Error("SetForMaxDeg", "MaxDeg < 0");
	double A = At(0,0), B = At(0,1), C = At(1, 0), D = At(1, 1), DeltaA = 0, DeltaB = 0;
	if (A > 1) Error("SetForMaxDeg", "A > 1"); if (A < 0) Error("SetForMaxDeg", "A < 0");
	if (B > 1) Error("SetForMaxDeg", "B > 1"); if (B < 0) Error("SetForMaxDeg", "B < 0");
	if (C > 1) Error("SetForMaxDeg", "C > 1"); if (C < 0) Error("SetForMaxDeg", "C < 0");
	if (D > 1) Error("SetForMaxDeg", "D > 1"); if (D < 0) Error("SetForMaxDeg", "D < 0");
	int BestRow = 0, BestCol = 0;
	double MaxExpDeg = GetMaxExpectedDeg(A,B,C,D, NIter, BestRow, BestCol);
	
	if (MaxExpDeg == MaxDeg) return;

	TIntPr Corner, Diag1, Diag2, Least;
	Corner.Val1 = BestRow; Corner.Val2 = BestCol;
	Diag1.Val1 = BestRow; Diag1.Val2 = (BestCol == 1) ? 0 : 1;
	Diag2.Val1 = (BestRow == 1) ? 0 : 1; Diag2.Val2 = BestCol;
	Least.Val1 = (BestRow == 1) ? 0 : 1; Least.Val2 = (BestCol == 1) ? 0 : 1;
	double CornerV = At(Corner.Val1, Corner.Val2),
		Diag1V = At(Diag1.Val1, Diag1.Val2),
		Diag2V = At(Diag2.Val1, Diag2.Val2),
		LeastV = At(Least.Val1, Least.Val2);

	bool DecFound = false;
	double Step = 0.001;
	double SensCoeff = 0.01; 
	if (MaxDeg < MaxExpDeg) 
		Step *= -1;

	// increase Corner, decrease Diag1 and Diag2
	// if MaxExpDeg begins to decrease, stop
	for (CornerV = CornerV + Step; CornerV <= 1; CornerV += Step){
		Diag1V = Diag1V - Step / 2;
		Diag2V = Diag2V - Step / 2;
		double CurrMaxDeg = pow(CornerV + Diag1V, NIter) / 2 + pow(CornerV + Diag2V, NIter) / 2 ;
		//printf("CurrMaxDeg = %d CornerV = %f\n", CurrMaxDeg, CornerV);
		if (CurrMaxDeg < MaxExpDeg){
			// reset values for previous step
			CornerV -= Step;
			Diag1V += Step / 2; Diag2V += Step / 2;
			break;
		}
		else 
			MaxExpDeg = CurrMaxDeg;
		if (abs(CurrMaxDeg - MaxDeg) / MaxDeg <= SensCoeff){
			DecFound = true;
			break;
		}
	}
	if (CornerV > 1) CornerV -= Step;
	printf("MaxDeg after step 1: %f\n", pow(CornerV + Diag1V, NIter) / 2 + pow(CornerV + Diag2V, NIter) / 2);

	// fix Corner, increase Diag1 and Diag2, decrease Least
	if (!DecFound){
		for (LeastV = LeastV - Step; LeastV >= 0; LeastV -= Step){
			Diag1V = Diag1V + Step / 2;
			Diag2V = Diag2V + Step / 2;
			MaxExpDeg = pow(CornerV + Diag1V, NIter) / 2 + pow(CornerV + Diag2V, NIter) / 2 ;
			if (abs(MaxExpDeg - MaxDeg) / MaxDeg <= SensCoeff){
				DecFound = true;
				break;
			}
		}
	}
	if (LeastV < 0) LeastV += Step;

	//printf("MaxDeg after step 2: %f\n", pow(CornerV + Diag1V, NIter) / 2 + pow(CornerV + Diag2V, NIter) / 2);
	//printf("CornerV = %f\n", CornerV);
	// increase max(Diag1, Diag2), decrease min(Diag1,Diag2)
	// if MaxExpDeg begins to decrease, stop
	if (!DecFound){
		double MaxDiag = (Diag1V > Diag2V) ? Diag1V : Diag2V,
			MinDiag = (MaxDiag == Diag1V) ? Diag2V : Diag1V;
		for (MaxDiag = MaxDiag + Step; MaxDiag <= 1; MaxDiag += Step){
			if (MinDiag - Step < 0)
				break;
			MinDiag -= Step;
			double CurrMaxDeg = pow(CornerV + MaxDiag, NIter) / 2 + pow(CornerV + MinDiag, NIter) / 2 ;
			//printf("CurrMaxDeg: %f MaxExpDeg: %f\n", CurrMaxDeg, MaxExpDeg);
			if (CurrMaxDeg < MaxExpDeg){
				// reset values for previous step
				if (Diag1 > Diag2){
					Diag1V = MaxDiag - Step;
					Diag2V = MinDiag + Step;
				}
				else {
					Diag2V = MaxDiag - Step;
					Diag1V = MinDiag + Step;
				}
				break;
				//printf("Err. CurrMaxDeg: %f MaxExpDeg: %f\n", CurrMaxDeg, MaxExpDeg);
			}
			else 
			MaxExpDeg = CurrMaxDeg;
			if (abs(MaxExpDeg - MaxDeg) / MaxDeg <= SensCoeff){
				DecFound = true;
				break;
			}
		}
	}

	//printf("MaxDeg after step 3: %f\n", pow(CornerV + Diag1V, NIter) / 2 + pow(CornerV + Diag2V, NIter) / 2);

	//int BestRow = (At(0,0) + At(0,1) > At(1,0) + At(1,1)) ? 0 : 1;
	//int LeastRow = (BestRow == 0) ? 1 : 0;
	//double A = At(BestRow,0), B = At(BestRow,1), C = At(LeastRow, 0), D = At(LeastRow, 1), DeltaA = 0, DeltaB = 0;
	//bool DecFound = false;
	//double MaxDeltaA = (1 - A < D) ? 1 - A : D, MaxDeltaB = (1 - B < C) ? 1 - B : C;
	//double Step = 0.001;
	//for (double i = Step; i <= MaxDeltaA; i+= Step){
	//	DeltaA = i;
	//	printf("%f\n", pow(A + DeltaA + B + DeltaB, NIter));
	//	if (static_cast<int>(pow(A + DeltaA + B, NIter)) == MaxDeg){
	//		DecFound = true;
	//		break;
	//	}
	//}
	//if (!DecFound){
	//	for (double i = Step; i <= MaxDeltaB; i+=Step){
	//		DeltaB = i;
	//		printf("%f\n", pow(A + DeltaA + B + DeltaB, NIter));
	//		if (static_cast<int>(pow(A + DeltaA + B + DeltaB, NIter)) == MaxDeg){
	//			DecFound = true;
	//			break;
	//		}
	//	}
	//}
	
	if (DecFound == false){
		printf("%f %f %f %f MaxDeg: %f Expected: %f\n", CornerV, Diag1V, Diag2V, Least, pow(CornerV + Diag1V, NIter) / 2 + pow(CornerV + Diag2V, NIter) / 2, MaxDeg);
		Error("SetForMaxDeg", "Cannot find solution");
	}
	printf("CornerV = %f\n", CornerV);
	if (CornerV > 1 || CornerV < 0)
		Error("SetForMaxDeg", "CornerV > 1 || CornerV < 0");
	if (Diag1V > 1 || Diag1V < 0)
		Error("SetForMaxDeg", "Diag1V > 1 || Diag1V < 0");
	if (Diag2V > 1 || Diag2V < 0)
		Error("SetForMaxDeg", "Diag2V > 1 || Diag2V < 0");
	if (LeastV > 1 || LeastV < 0)
		Error("SetForMaxDeg", "LeastV > 1 || LeastV < 0");
	/*if (A + DeltaA > 1)
		Error("SetForMaxDeg", "A + DeltaA > 1");
	if (B + DeltaB > 1)
		Error("SetForMaxDeg", "B + DeltaB > 1");
	if (C - DeltaB < 0)
		Error("SetForMaxDeg", "C - DeltaB < 0");
	if (D - DeltaA < 0)
		Error("SetForMaxDeg", "A + DeltaA > 1");*/
	At(Corner.Val1,Corner.Val2) = CornerV; At(Diag1.Val1, Diag1.Val2) = Diag1V; At(Diag2.Val1,Diag2.Val2) = Diag2V; At(Least.Val1, Least.Val2) = LeastV; 
	//printf("Maximum expected degree: %f\n", GetMaxExpectedDeg(NIter));
}

void TKronMtx::AddRndNoise(const double& SDev) {
  Dump("before");
  double NewVal;
  int c =0;
  for (int i = 0; i < Len(); i++) {
    for(c = 0; ((NewVal = At(i)*Rnd.GetNrmDev(1, SDev, 0.8, 1.2)) < 0.01 || NewVal>0.99) && c <1000; c++) { }
    if (c < 999) { At(i) = NewVal; } else { printf("XXXXX\n"); }
  }
  Dump("after");
}

TStr TKronMtx::GetMtxStr() const {
  TChA ChA("[");
  for (int i = 0; i < Len(); i++) {
    ChA += TStr::Fmt("%g", At(i));
    if ((i+1)%GetDim()==0 && (i+1<Len())) { ChA += "; "; }
    else if (i+1<Len()) { ChA += ", "; }
  }
  ChA += "]";
  return TStr(ChA);
}

void TKronMtx::ToOneMinusMtx() {
  for (int i = 0; i < Len(); i++) {
    IAssert(At(i) >= 0.0 && At(i) <= 1.0);
    At(i) = 1.0 - At(i);
  }
}

void TKronMtx::GetLLMtx(TKronMtx& LLMtx) {
  LLMtx.GenMtx(MtxDim);
  for (int i = 0; i < Len(); i++) {
    if (At(i) != 0.0) { LLMtx.At(i) = log(At(i)); }
    else { LLMtx.At(i) = NInf; }
  }
}

void TKronMtx::GetProbMtx(TKronMtx& ProbMtx) {
  ProbMtx.GenMtx(MtxDim);
  for (int i = 0; i < Len(); i++) {
    if (At(i) != NInf) { ProbMtx.At(i) = exp(At(i)); }
    else { ProbMtx.At(i) = 0.0; }
  }
}

void TKronMtx::Swap(TKronMtx& KronMtx) {
  ::Swap(MtxDim, KronMtx.MtxDim);
  SeedMtx.Swap(KronMtx.SeedMtx);
}

int TKronMtx::GetNodes(const int& NIter) const {
  return (int) pow(double(GetDim()), double(NIter));
}

int TKronMtx::GetEdges(const int& NIter) const {
  return (int) pow(double(GetMtxSum()), double(NIter));
}

int TKronMtx::GetKronIter(const int& Nodes) const {
  return (int) ceil(log(double(Nodes)) / log(double(GetDim()))); // upper bound
  //return (int) TMath::Round(log(double(Nodes)) / log(double(GetDim()))); // round to nearest power
}

int TKronMtx::GetNZeroK(const PNGraph& Graph) const {
 return GetNodes(GetKronIter(Graph->GetNodes()));
}

double TKronMtx::GetEZero(const int& Edges, const int& KronIters) const {
  return pow((double) Edges, 1.0/double(KronIters));
}

double TKronMtx::GetMtxSum() const {
  double Sum = 0;
  for (int i = 0; i < Len(); i++) {
    Sum += At(i); }
  return Sum;
}

double TKronMtx::GetRowSum(const int& RowId) const {
  double Sum = 0;
  for (int c = 0; c < GetDim(); c++) {
    Sum += At(RowId, c); }
  return Sum;
}

double TKronMtx::GetColSum(const int& ColId) const {
  double Sum = 0;
  for (int r = 0; r < GetDim(); r++) {
    Sum += At(r, ColId); }
  return Sum;
}

double TKronMtx::GetEdgeProb(int NId1, int NId2, const int& NKronIters) const {
  double Prob = 1.0;
  for (int level = 0; level < NKronIters; level++) {
    Prob *= At(NId1 % MtxDim, NId2 % MtxDim);
    if (Prob == 0.0) { return 0.0; }
    NId1 /= MtxDim;  NId2 /= MtxDim;
  }
  return Prob;
}

double TKronMtx::GetNoEdgeProb(int NId1, int NId2, const int& NKronIters) const {
  return 1.0 - GetEdgeProb(NId1, NId2, NKronIters);
}

double TKronMtx::GetEdgeLL(int NId1, int NId2, const int& NKronIters) const {
  double LL = 0.0;
  for (int level = 0; level < NKronIters; level++) {
    const double& LLVal = At(NId1 % MtxDim, NId2 % MtxDim);
    if (LLVal == NInf) return NInf;
    LL += LLVal;
    NId1 /= MtxDim;  NId2 /= MtxDim;
  }
  return LL;
}

double TKronMtx::GetNoEdgeLL(int NId1, int NId2, const int& NKronIters) const {
  return log(1.0 - exp(GetEdgeLL(NId1, NId2, NKronIters)));
}

// 2nd order Taylor approximation log(1-x) ~ -x - 0.5x^2
double TKronMtx::GetApxNoEdgeLL(int NId1, int NId2, const int& NKronIters) const {
  const double EdgeLL = GetEdgeLL(NId1, NId2, NKronIters);
  return -exp(EdgeLL) - 0.5*exp(2*EdgeLL);
}

bool TKronMtx::IsEdgePlace(int NId1, int NId2, const int& NKronIters, const double& ProbTresh) const {
  double Prob = 1.0;
  for (int level = 0; level < NKronIters; level++) {
    Prob *= At(NId1 % MtxDim, NId2 % MtxDim);
    if (ProbTresh > Prob) { return false; }
    NId1 /= MtxDim;  NId2 /= MtxDim;
  }
  return true;
}

// deriv a*log(x) = a/x
double TKronMtx::GetEdgeDLL(const int& ParamId, int NId1, int NId2, const int& NKronIters) const {
  const int ThetaX = ParamId % GetDim();
  const int ThetaY = ParamId / GetDim();
  int ThetaCnt = 0;
  for (int level = 0; level < NKronIters; level++) {
    if ((NId1 % MtxDim) == ThetaX && (NId2 % MtxDim) == ThetaY) {
      ThetaCnt++; }
    NId1 /= MtxDim;  NId2 /= MtxDim;
  }
  return double(ThetaCnt) / exp(At(ParamId));
}

// deriv log(1-x^a*y^b..) = -x'/(1-x) = (-a*x^(a-1)*y^b..) / (1-x^a*y^b..)
double TKronMtx::GetNoEdgeDLL(const int& ParamId, int NId1, int NId2, const int& NKronIters) const {
  const int& ThetaX = ParamId % GetDim();
  const int& ThetaY = ParamId / GetDim();
  int ThetaCnt = 0;
  double DLL = 0, LL = 0;
  for (int level = 0; level < NKronIters; level++) {
    const int X = NId1 % MtxDim;
    const int Y = NId2 % MtxDim;
    const double LVal = At(X, Y);
    if (X == ThetaX && Y == ThetaY) {
      if (ThetaCnt != 0) { DLL += LVal; }
      ThetaCnt++;
    } else { DLL += LVal; }
    LL += LVal;
    NId1 /= MtxDim;  NId2 /= MtxDim;
  }
  return -ThetaCnt*exp(DLL) / (1.0 - exp(LL));
}

// 2nd order Taylor approximation log(1-x) ~ -x - 0.5x^2
double TKronMtx::GetApxNoEdgeDLL(const int& ParamId, int NId1, int NId2, const int& NKronIters) const {
  const int& ThetaX = ParamId % GetDim();
  const int& ThetaY = ParamId / GetDim();
  int ThetaCnt = 0;
  double DLL = 0;//, LL = 0;
  for (int level = 0; level < NKronIters; level++) {
    const int X = NId1 % MtxDim;
    const int Y = NId2 % MtxDim;
    const double LVal = At(X, Y); IAssert(LVal > NInf);
    if (X == ThetaX && Y == ThetaY) {
      if (ThetaCnt != 0) { DLL += LVal; }
      ThetaCnt++;
    } else { DLL += LVal; }
    //LL += LVal;
    NId1 /= MtxDim;  NId2 /= MtxDim;
  }
  //return -ThetaCnt*exp(DLL)*(1.0 + exp(LL)); // -x'/(1+x) WRONG!
  // deriv = -(ax^(a-1)*y^b..) - a*x^(2a-1)*y^2b..
  //       = - (ax^(a-1)*y^b..) - a*x*(x^(a-1)*y^b..)^2
  return -ThetaCnt*exp(DLL) - ThetaCnt*exp(At(ThetaX, ThetaY)+2*DLL);
}

uint TKronMtx::GetNodeSig(const double& OneProb) {
  uint Sig = 0;
  for (int i = 0; i < (int)(8*sizeof(uint)); i++) {
    if (TKronMtx::Rnd.GetUniDev() < OneProb) {
      Sig |= (1u<<i); }
  }
  return Sig;
}

double TKronMtx::GetEdgeProb(const uint& NId1Sig, const uint& NId2Sig, const int& NIter) const {
  Assert(GetDim() == 2);
  double Prob = 1.0;
  for (int i = 0; i < NIter; i++) {
    const uint Mask = (1u<<i);
    const uint Bit1 = NId1Sig & Mask;
    const uint Bit2 = NId2Sig & Mask;
    Prob *= At(int(Bit1!=0), int(Bit2!=0));
  }
  return Prob;
}

PNGraph TKronMtx::GenThreshGraph(const double& Thresh) const {
  PNGraph Graph = TNGraph::New();
  for (int i = 0; i < GetDim(); i++) {
    Graph->AddNode(i); }
  for (int r = 0; r < GetDim(); r++) {
    for (int c = 0; c < GetDim(); c++) {
      if (At(r, c) >= Thresh) { Graph->AddEdge(r, c); }
    }
  }
  return Graph;
}

PNGraph TKronMtx::GenRndGraph(const double& RndFact) const {
  PNGraph Graph = TNGraph::New();
  for (int i = 0; i < GetDim(); i++) {
    Graph->AddNode(i); }
  for (int r = 0; r < GetDim(); r++) {
    for (int c = 0; c < GetDim(); c++) {
      if (RndFact * At(r, c) >= TKronMtx::Rnd.GetUniDev()) { Graph->AddEdge(r, c); }
    }
  }
  return Graph;
}

int TKronMtx::GetKronIter(const int& GNodes, const int& SeedMtxSz) {
  return (int) ceil(log(double(GNodes)) / log(double(SeedMtxSz)));
}

// slow but exaxt procedure (we flip all O(N^2) edges)
PNGraph TKronMtx::GenKronecker(const TKronMtx& SeedMtx, const int& NIter, const bool& IsDir, const int& Seed) {
  const TKronMtx& SeedGraph = SeedMtx;
  const int NNodes = SeedGraph.GetNodes(NIter);
  printf("  Kronecker: %d nodes, %s...\n", NNodes, IsDir ? "Directed":"UnDirected");
  PNGraph Graph = TNGraph::New(NNodes, -1);
  TExeTm ExeTm;
  TRnd Rnd(Seed);
  int edges = 0;
  for (int node1 = 0; node1 < NNodes; node1++) {
    Graph->AddNode(node1); }
  if (IsDir) {
    for (int node1 = 0; node1 < NNodes; node1++) {
      for (int node2 = 0; node2 < NNodes; node2++) {
        if (SeedGraph.IsEdgePlace(node1, node2, NIter, Rnd.GetUniDev())) {
          Graph->AddEdge(node1, node2);
          edges++;
        }
      }
      if (node1 % 1000 == 0) printf("\r...%dk, %dk", node1/1000, edges/1000);
    }
  } else {
    for (int node1 = 0; node1 < NNodes; node1++) {
      for (int node2 = node1; node2 < NNodes; node2++) {
        if (SeedGraph.IsEdgePlace(node1, node2, NIter, Rnd.GetUniDev())) {
          Graph->AddEdge(node1, node2);
          Graph->AddEdge(node2, node1);
          edges++;
        }
      }
      if (node1 % 1000 == 0) printf("\r...%dk, %dk", node1/1000, edges/1000);
    }
  }
  printf("\r             %d edges [%s]\n", Graph->GetEdges(), ExeTm.GetTmStr());
  return Graph;
}

double TKronMtx::GetEigMax() const{
	double A = At(0,0), B = At(0,1), C = At(1,0), D = At(1,1);
	double X = (A + D) / 2, Y =  sqrt(pow(A+D,2) - 4 * (A * D - B*C)) / 2;
	double EigMax = X + Y;
	/*double EigMin = X - Y;
	double Member = sqrt(NIter * pow(A + B, NIter));
	printf("EigMax + Member: %f\n", EigMax + Member);
	printf("EigMin: %f\n", EigMin);*/
	return EigMax;
}

double TKronMtx::GetEigMin() const {
	double A = At(0,0), B = At(0,1), C = At(1,0), D = At(1,1);
	double X = (A + D) / 2, Y =  sqrt(pow(A+D,2) - 4 * (A * D - B*C)) / 2;
	double EigMin = X - Y;
	return EigMin;
}

int TKronMtx::AddEdges(const TKronMtx& SeedMtx, const int&NIter, const bool& IsDir, TRnd& Rnd, PNGraph& G, const int& NEdges, const int&InDegMax, const int& OutDegMax, double ModelClustCf){
	const int MtxDim = SeedMtx.GetDim();
	const double MtxSum = SeedMtx.GetMtxSum();
	const int NNodes = SeedMtx.GetNodes(NIter);
	printf("  AddEdges: %d nodes, %d edges, %s...\n", NNodes, NEdges, IsDir ? "Directed":"UnDirected");
	// prepare cell probability vector
	TVec<TFltIntIntTr> ProbToRCPosV; // row, col position
	double CumProb = 0.0;
	for (int r = 0; r < MtxDim; r++) {
		for (int c = 0; c < MtxDim; c++) {
			const double Prob = SeedMtx.At(r, c);
			if (Prob > 0.0) {
				CumProb += Prob;
				ProbToRCPosV.Add(TFltIntIntTr(CumProb/MtxSum, r, c));
				printf("%d%d %f\n", r, c, CumProb/MtxSum);
			}
		}
	}
	
	// add edges
	int Rng, Row, Col, Collision=0, n = 0, edges = 0, ClosedTriads = 0, ClustCollision = 0;
	for (edges = 0; edges < NEdges; ) {
		Rng=NNodes;  Row=0;  Col=0;
		for (int iter = 0; iter < NIter; iter++) {
			const double& Prob = Rnd.GetUniDev();
			n = 0; while(Prob > ProbToRCPosV[n].Val1) { n++; }
			const int MtxRow = ProbToRCPosV[n].Val2;
			const int MtxCol = ProbToRCPosV[n].Val3;
			Rng /= MtxDim;
			Row += MtxRow * Rng;
			Col += MtxCol * Rng;
		}
		if (! G->IsEdge(Row, Col) && Row != Col) { // allow self-loops
			// check for consistency for maximum degrees
			if (G->GetNI(Row).GetOutDeg() < OutDegMax && G->GetNI(Col).GetInDeg() < InDegMax){
				G->AddEdge(Row, Col);  edges++;
				if (! IsDir) {
					G->AddEdge(Col, Row);
					edges++;
				}
				if (CheckClustCf(G, Row, Col, ModelClustCf, Rnd, false, ClustCollision)){
				  edges++;
				  edges++;
				  ClosedTriads++;
			    }
			}
			else { Collision++; }
		} else { Collision++; }
		//if (edges % 1000 == 0) printf("\r...%dk", edges/1000);
	}
	std::string s = "Edges added=" + std::to_string((long long)edges) +", edges to add=" + std::to_string((long long)NEdges) + "\n";
	printf("%s",s.c_str());
	printf("AddEdges: ClosedTriads %d ClustCollision %d\n", ClosedTriads, ClustCollision);
	//PrintDeg(G, "AddEdges");
	//printf("             %d edges [%s]\n", Graph->GetEdges(), ExeTm.GetTmStr());
	return Collision;
}

void TKronMtx::SetForMaxEigen(const double K, const int& NIter){
	double S = GetMtxSum();
	double A = At(0,0), B = At(0,1), C = At(1,0), D = At(1,1);
	double Sum1 = pow(K, 1.00 / NIter) * (A + D),
		Sum2 = pow(K, 2.00 / NIter) * (A + C),
		Sum3 = pow(K, 2.00 / NIter) * B * C,
		Sum4 = S;
	printf("Sum1: %f, Sum2: %f, Sum3: %f, Sum4: %f\n", Sum1, Sum2, Sum3, Sum4);
	double Disc = pow (2 * Sum2 - Sum1 - Sum4, 2) * D * D - 4 * D * D * (Sum3 - (Sum2 - Sum1) * (Sum4 - Sum2));
	double DeltaFM = -1 * D * (2 * Sum2 - Sum1 - Sum4) / (2 * D * D),
		DeltaSM = sqrt(Disc) / (2 * D * D);
	double Delta1 = DeltaFM - DeltaSM, Delta2 = DeltaFM + DeltaSM;
	double Delta = Delta1 * D < 1 ? Delta1 : Delta2;
	double Alpha = (Sum1 - D * Delta) / A;
	double Beta = (Sum2 - A * Alpha) / B;
	double Gamma = Sum3 / (Beta * B) / C;
	double NewA = A * Alpha, NewB = B * Beta, NewC = C * Gamma, NewD = D * Delta;
	if (NewA > 1 && NewD > 1)
		printf("Error: NewA > 1 && NewD > 1\n");
	if (NewB > 1 && NewC > 1)
		printf("Error: NewB > 1 && NewC > 1\n");
	/*if (NewA > 1) {NewD += NewA - 1; if (NewD > 1) printf("Error. NewD > 1\n"); NewA = 1;}
	if (NewB > 1) {NewC += NewB - 1; if (NewC > 1) printf("Error. NewC > 1\n"); NewB = 1;}
	if (NewC > 1) {NewB += NewC - 1; if (NewB > 1) printf("Error. NewB > 1\n"); NewC = 1;}
	if (NewD > 1) {NewA += NewD - 1; if (NewA > 1) printf("Error. NewA > 1\n"); NewD = 1;}*/
	double NewSum = NewA + NewB + NewC + NewD;
	At(0,0) = NewA; At(0,1) = NewB; At(1,0) = NewC; At(1,1) = NewD;
	return;
}

void TKronMtx::Transpose(TKronMtx& SeedMtx){
	TKronMtx Mtx(SeedMtx);
	const int Dim = SeedMtx.GetDim();
	for (int r = 0; r < Dim; r++){
		for (int c = 0; c < Dim; c++)
			SeedMtx.At(r,c) = Mtx.At(c,r);
	}
}

void TKronMtx::GetRowProbCumV(const TKronMtx& Mtx, TVec<TVec<TFltIntIntTr>>& RowProbCumV){
	const int MtxDim = Mtx.GetDim();
	for (int r = 0; r < MtxDim; r++) {
		double CumProb = 0.0, RowSum = 0;
		for (int c = 0; c < MtxDim; c++) RowSum += Mtx.At(r,c);
		TVec<TFltIntIntTr> ProbCumV;
		for (int c = 0; c < MtxDim; c++) {
			const double Prob = Mtx.At(r, c);
			if (Prob > 0.0) {
				CumProb += Prob;
				ProbCumV.Add(TFltIntIntTr(CumProb/RowSum, r, c));
			}
		}
		RowProbCumV.Add(ProbCumV);
	}
}

int GetCol(const TVec<TVec<TFltIntIntTr>>& RowProbCumV, const int &Row, const int& NIter, TRnd& Rnd){
	return 0;
}

int GetCol(const TVec<TVec<TVec<TFltIntIntTr>>>& RowProbCumV, const int &Row, const int& NIter, TRnd& Rnd){
	if (RowProbCumV.Len() == 0)
		Error("GetCol", "Empty probability vector");
	int n = 0, MtxDim = RowProbCumV[0].Len(); 
	int RngBegin = 0, RngEnd = pow(double(MtxDim), NIter)-1, PartRngSize = ( RngEnd - RngBegin + 1 ) / MtxDim, RowC = Row;
	for (int iter = 0; iter < NIter; iter++) {
		const double& Prob = Rnd.GetUniDev();
		int RowProb = RowC  / PartRngSize;
		RowC -= RowProb * PartRngSize;
		n = 0; while(Prob > RowProbCumV[iter][RowProb][n].Val1) { n++; }
		const int MtxCol = RowProbCumV[iter][RowProb][n].Val3;
		RngBegin += MtxCol * PartRngSize;
		RngEnd = RngBegin + PartRngSize - 1;
		//printf("ProbRow=%d PartRngSize=%d ", RowProb, PartRngSize);
		PartRngSize /= MtxDim;
		//printf("MtxCol=%d RngBegin=%d RngEnd=%d\n", MtxCol, RngBegin, RngEnd);
	}
	try {
		if (RngBegin != RngEnd){
			throw "GetCol() error. Failed to find column index, RngBegin=" + std::to_string((long long)RngBegin) + ", RngEnd=" + std::to_string((long long)RngEnd) + "\n";
		}
		//printf("\n");
		return RngBegin;
	}
	catch(const TStr& ex){
		printf("%s\n", ex.CStr());
		system("pause");
		exit(0);
	}
}

bool CheckEdges(TIntV& DegCount, TInt& S, const int&DegAdd, const int&DegToCheck, const int& LeastEdges){
	try {
		if (DegAdd + 1 <= DegToCheck) {
			DegCount[DegAdd] -= 1; 
			S -= 1;
			if (DegAdd + 1 != DegToCheck) 
				DegCount[DegAdd+1] += 1;
			return true;
		}
		else {
			if (LeastEdges < S){
				std::string s = "CheckEdges() error. Least edges < S: Least=" + std::to_string((long long)LeastEdges) + ", S=" + 
					std::to_string((long long)S) + "\n";
				throw s;
			}
			if (LeastEdges == S)
				return false;
			return true;
		}
	}
	catch (const std::string& ex){
		printf("%s\n", ex.c_str());
		system("pause");
		exit(0);
	}
}

int TKronMtx::AddSecondDir(bool IsOut, const TIntPr& InDegR, const TIntPr& OutDegR, PNGraph& G, const TKronMtx& SeedMtx, const int&NIter, TRnd&Rnd){
	int Collision = 0;
	const int NNodes = SeedMtx.GetNodes(NIter);
	const int NEdges = SeedMtx.GetEdges(NIter);
	const TInt& InMin = InDegR.Val1, & InMax = InDegR.Val2, & OutMin = OutDegR.Val1, & OutMax = OutDegR.Val2;
	TKronMtx Mtx(SeedMtx);
	// for in-degree transpose the initiator matrix
	if (!IsOut) Transpose(Mtx);
	// get row prob accum vectors
	TVec<TVec<TFltIntIntTr>> RowProbCumV;
	GetRowProbCumV(Mtx, RowProbCumV);
	
	TInt EdgesToAdd = IsOut ? OutMin : InMin;
	
	int EdgesAdded = 0; //int Rows = 0;

	for (int i = 0; i < NNodes; i++){
		int Row = i;
		int Deg = IsOut? G->GetNI(Row).GetOutDeg(): G->GetNI(Row).GetInDeg();
		//if (Deg > 2) printf("Deg > 2: %d\n", Deg);
		if (Deg >= EdgesToAdd) continue;
		//Rows++;
		//printf("\nRow %d Deg %d", Row, Deg);
		for (int j = 0; j < EdgesToAdd - Deg; j++){
			Row = i;
			int Col = GetCol(RowProbCumV, Row, NIter, Rnd);
			if (!IsOut) {TInt Add = Row; Row = Col; Col = Add;}
			if (Row != Col && !G->IsEdge(Row, Col)){
				int InDeg = G->GetNI(Col).GetInDeg(), OutDeg = G->GetNI(Row).GetOutDeg();
				if (InDeg + 1 > InMax || OutDeg + 1 > OutMax) {
					Collision++; j--; continue;
				}
				G->AddEdge(Row,Col);
				//printf("\nEdge %d,%d was added", Row, Col);
				EdgesAdded++;
			}
			else {Collision++; j--;}
		}

	}
	//printf("Rows: %d\n", Rows);
	std::string s = "Edges added=" + std::to_string((long long)EdgesAdded) +", edges to add=" + std::to_string((long long)EdgesToAdd * NNodes) + "\n";
	printf("%s",s.c_str());
	printf("Collisions (AddSecondDir): %d\n", Collision);
	//PrintDeg(G, "AddSecond");
	return Collision;
}


int TKronMtx::AddFirstDir(bool IsOut, const TIntPr& InDegR, const TIntPr& OutDegR, PNGraph& G, const TKronMtx& SeedMtx, const int&NIter, TRnd&Rnd){
	int Collision = 0;
	const int NNodes = SeedMtx.GetNodes(NIter);
	const int NEdges = SeedMtx.GetEdges(NIter);
	const TInt& InMin = InDegR.Val1, & InMax = InDegR.Val2, & OutMin = OutDegR.Val1, & OutMax = OutDegR.Val2;
	if (InMin == 0) return 0;
	TKronMtx Mtx(SeedMtx);
	// for in-degree transpose the initiator matrix
	if (!IsOut) Transpose(Mtx);
	// get row prob accum vectors
	TVec<TVec<TFltIntIntTr>> RowProbCumV;
	GetRowProbCumV(Mtx, RowProbCumV);
	TIntV DegCount;
	TInt S, EdgesToAdd, DegToCheck;
	if (IsOut){
		S = NNodes * InMin; EdgesToAdd = OutMin; DegToCheck = InMin; 
	}
	else {
		S = NNodes * OutMin; EdgesToAdd = InMin; DegToCheck = OutMin;
	}
	for (int i = 0; i < DegToCheck; i++)
		DegCount.Add(0);
	DegCount[0] = DegToCheck * NNodes;

	int EdgesAdded = 0;

	TIntV vec;
	for (int i = 0; i < G->GetNodes(); i++)
		vec.Add(G->GetNodes() - i - 1);

	for (int i = 0; i < EdgesToAdd; i++){
		for (int j = 0; j < NNodes; j++){
			int Row = vec[j];
			int Col = GetCol(RowProbCumV, Row, NIter, Rnd);
			if (!IsOut) {TInt Add = Row; Row = Col; Col = Add;}
			if (Row != Col && !G->IsEdge(Row, Col)){
				int InDeg = G->GetNI(Col).GetInDeg(), OutDeg = G->GetNI(Row).GetOutDeg();
				if (InDeg + 1 > InMax || OutDeg + 1 > OutMax) {Collision++; j--; continue;}
				// check for non-violation of EdgesToCheck
				int DegAdded = IsOut ? InDeg : OutDeg;
				
				/*printf("Least edges %d, S %d\n", NEdges - EdgesAdded, S.Val);
				if (NEdges - EdgesAdded < S.Val)
					printf("error\n");*/
				//printf("DegAdded: %d, DegToCheck: %d\n", DegAdded, DegToCheck.Val);
				if (CheckEdges(DegCount, S, DegAdded, DegToCheck, NEdges - EdgesAdded)){
					G->AddEdge(Row,Col);
					if (InDeg > 2) printf("InDeg > 2: %d\n", InDeg);
					EdgesAdded++;
				}
				else {Collision++; j--;}
			}
			else {Collision++; j--;}
		}
	}
	printf("Collisions (AddFirstDir): %d\n", Collision);
	std::string s = "Edges added=" + std::to_string((long long)EdgesAdded) +", edges to add=" + std::to_string((long long)EdgesToAdd * NNodes) + "\n";
	printf("%s",s.c_str());
	TSnap::PlotHops(G, "AddFirstHops.tab");
	/*for (int i = 0; i < G->GetNodes(); i++)
		printf("Node %d: InDeg %d, OutDeg %d\n", i, G->GetNI(i).GetInDeg(),G->GetNI(i).GetOutDeg());*/
	//PrintDeg(G, "AddFirst");

	try {
		if (EdgesAdded != EdgesToAdd * NNodes){
			throw "AddFirst() error. Edges added=" + std::to_string((long long)EdgesAdded) +", edges to add=" + std::to_string((long long)EdgesToAdd * NNodes) + "\n";
		}
	}
	catch (const TStr& ex){
		printf("%s\n", ex.CStr());
		system("pause");
		exit(0);
	}

	return Collision;
}

int TKronMtx::AddUnDir(const TIntPr& DegR, PNGraph& G, const TKronMtx& SeedMtx, const int& NIter, TRnd& Rnd, double NoiseCoeff){
	int Collision = 0;
	const int NNodes = SeedMtx.GetNodes(NIter);
	const int NEdges = SeedMtx.GetEdges(NIter);
	const int Dim = SeedMtx.GetDim();
	const TInt& DegReq = DegR.Val1;
	TKronMtx Mtx(SeedMtx);

	// get NIter vectors of cumulative probabilities
	TVec<TVec<TFltIntIntTr>> ProbToRCPosV;
    double AvgExpectedDeg = GetNoisedProbV(ProbToRCPosV, NoiseCoeff, Rnd, NIter, SeedMtx);
	// get row prob accum vectors [NIter x MtxDim * MtxDim]
	TVec<TVec<TVec<TFltIntIntTr>>> RowProbCumV;
	for (int i = 0; i < NIter; i++){
		TVec<TVec<TFltIntIntTr>> IterRowProbCumV;
		TVec<TFltIntIntTr> Instance;
		// to get row prob accum, we should subtract sum from previous row
		double RowSum = 0.0;
		for (int j = 0; j < Dim * Dim; j++){
			TFltIntIntTr T(ProbToRCPosV[i][j]); 
			if (j != 0) T.Val1 -= ProbToRCPosV[i][j-1].Val1;
			RowSum += T.Val1;
			Instance.Add(T);
			//printf("T.Val1 = %f\n", T.Val1);
			// if it is the last element of the row
			if ( (j + 1) % Dim == 0) {
				for (int k = 0; k < Dim; k++) {
					Instance[k].Val1 /= RowSum;
					if (k != 0) Instance[k].Val1 += Instance[k-1].Val1;
					//printf("Instance[%d] = %f\n", k, Instance[k].Val1);
				}
				IterRowProbCumV.Add(Instance);
				Instance.Clr();
				RowSum = 0;
			}
		}
		RowProbCumV.Add(IterRowProbCumV);
	}
	//GetRowProbCumV(Mtx, RowProbCumV); Transpose(Mtx); GetRowProbCumV(Mtx, ColProbCumV);
	int EdgesAdded = 0, ClustCollision = 0, ClosedTriads = 0;

	/*vector<int> NodeIDs;
	for (int i = 0; i < NNodes; i++) NodeIDs.push_back(i);
	for (int i = 0; i < NNodes; i++) {
		int Node1 = Rnd.GetUniDev() * (NNodes-1),
			Node2 = Rnd.GetUniDev() * (NNodes-1);
		int Add = NodeIDs[Node1];
		NodeIDs[Node1] = NodeIDs[Node2];
		NodeIDs[Node2] = Add;
	}*/
	
	int MissedRows = 0;

	for (int i = 0; i < NNodes; i++){
		int Row = i;
		int InDeg = G->GetNI(Row).GetInDeg(), OutDeg = G->GetNI(Row).GetOutDeg();
		if (InDeg != OutDeg) printf("InDeg != OutDeg. Error");
		//printf("InDeg %d OutDeg %d\n", InDeg, OutDeg);
		if (InDeg >= DegReq || OutDeg >= DegReq) {
			MissedRows++; 
			continue;
		}
		for (int j = 0; j < DegReq-OutDeg; j++){	
			//printf("Collision %d\ edges added %d \n", Collision, EdgesAdded);
			int Col = GetCol(RowProbCumV, Row, NIter, Rnd);
			if (Row != Col && !G->IsEdge(Row, Col)){
				//int InDegCol = G->GetNI(Col).GetInDeg(), OutDegCol = G->GetNI(Col).GetOutDeg();
				//if (InDegCol + 1 > DegReq || OutDegCol + 1 > DegReq) {Collision++;  j--; continue;}
				G->AddEdge(Row,Col);
				//printf("(%d %d)\t", Row, Col);
				EdgesAdded++;
				G->AddEdge(Col, Row);
				//printf("(%d %d)\t", Col, Row);
				EdgesAdded++;
				
			}
			else {Collision++; j--;}//printf("Collision1\n"); }	
		}
	}
	std::string s = "Edges added=" + std::to_string((long long)EdgesAdded) +", edges to add=" + std::to_string((long long) DegReq * NNodes) + "\n";
	printf("%s",s.c_str());
	//printf("Missed rows: %d\n", MissedRows); system("pause");
	//printf("AddUnDir: ClosedTriads %d ClustCollision %d\n", ClosedTriads, ClustCollision);
	//PrintDeg(G, "AddFirst");

	/*try {
		if (EdgesAdded != 2 * DegMin * NNodes){
			throw "AddUnDir() error. Edges added=" + std::to_string((long long)EdgesAdded) +", edges to add=" + std::to_string((long long)2 * DegMin * NNodes) + "\n";
		}
	}
	catch (const TStr& ex){
		printf("%s\n", ex);
		system("pause");
		exit(0);
	}*/

	return Collision;
}

//void TKronMtx::DelEdges(PNGraph& Graph, const TKronMtx& Mtx, const int& NIter, const int& NEdges){
//  TRnd Rnd(0);
//  TExeTm ExeTm;
//  // prepare cell probability vector
//  TVec<TFltIntIntTr> ProbToRCPosV; // row, col position
//  double CumProb = 0.0;
//  for (int r = 0; r < MtxDim; r++) {
//    for (int c = 0; c < MtxDim; c++) {
//      const double Prob = Mtx.At(r, c);
//      if (Prob > 0.0) {
//        CumProb += Prob;
//        ProbToRCPosV.Add(TFltIntIntTr(CumProb/Mtx.GetMtxSum(), r, c));
//      }
//    }
//  }
// 
//  int NNodes = pow(static_cast<double>(2), NIter);
//
//  int Rng, Row, Col, Collision=0, n = 0;
//  for (int edges = 0; edges < NEdges; ) {
//    Rng=NNodes;  Row=0;  Col=0;
//    for (int iter = 0; iter < NIter; iter++) {
//      const double& Prob = Rnd.GetUniDev();
//      n = 0; while(Prob > ProbToRCPosV[n].Val1) { n++; }
//      const int MtxRow = ProbToRCPosV[n].Val2;
//      const int MtxCol = ProbToRCPosV[n].Val3;
//      Rng /= MtxDim;
//      Row += MtxRow * Rng;
//      Col += MtxCol * Rng;
//    }
//    if (Graph->IsEdge(Row, Col)) { // allow self-loops
//      Graph->DelEdge(Row, Col,false);  edges++;
//    } else { Collision++; }
//    //if (edges % 1000 == 0) printf("\r...%dk", edges/1000);
//  }
//}

void TKronMtx::RemoveZeroDegreeNodes(PNGraph& out, const TKronMtx& Mtx, const int& NIter, const int& MinDeg, const int&MaxDeg){
	TRnd rnd;
	rnd.GetUniDev();
	int nodesCount = out->GetNodes();
	int edgesCount = out->GetEdges();
	TVec<TVec<TFltIntIntTr>> RowProbCumV, ColProbCumV;
	//TIntV DegBefore, DegAfter; TFltV Prob;
	int MaxExpectedDeg = 0;
	
	for (int i = 0; i < nodesCount; i++){
		if (out->GetNI(i).GetInDeg() > MaxExpectedDeg)
			MaxExpectedDeg = out->GetNI(i).GetInDeg();
	}

	GetRowProbCumV(Mtx, RowProbCumV);

	
	int MaxAddDeg = 0;
	int EdgesToDel = 0;
	for (int i = 0; i < nodesCount; i++){
		int InDeg = out->GetNI(i).GetInDeg();
		
		if (InDeg < MinDeg){
			int EdgesToAdd = MinDeg - InDeg + rnd.GetUniDev() * (MinDeg - InDeg + 2 * MinDeg);
			if (InDeg + EdgesToAdd > MaxDeg) EdgesToAdd = InDeg - MaxDeg;
			EdgesToDel += EdgesToAdd;
			//printf("Edges to add: %d\n", EdgesToAdd);
			for (int j = 0; j < EdgesToAdd; j++){
				while (1){
					double val = rnd.GetUniDev();
					// get neighbour node using probability matrix
					int nodeId = GetCol(RowProbCumV, i, NIter, rnd);
					//int nodeId = rnd.GetUniDev() * (nodesCount-1);
					if (nodeId == i || out->IsEdge(nodeId,i)) continue;
					auto NI = out->GetNI(nodeId);
					int OutDeg = NI.GetOutDeg();
					if (OutDeg == 0 || OutDeg >= MaxDeg) continue;
					int NbId = static_cast<int>(rnd.GetUniDev() * (OutDeg - 1));
					int NbOutId = NI.GetNbrNId(NbId);
					if (out->GetNI(NbOutId).GetInDeg() <= MinDeg || (MaxExpectedDeg - out->GetNI(NbOutId).GetInDeg()) / static_cast<double>(MaxExpectedDeg) < 0.01  ) 
						continue;

					out->DelEdge(nodeId, NbOutId, false);

					/*printf("Edge (%d,%d) was deleted\n", nodeId, NbOutId);
					printf("Edge (%d,%d) was added\n", nodeId, i);*/

					out->AddEdge(nodeId, i);
					out->AddEdge(i,nodeId);
					//printf("Degrees after: i: %d node: %d nb: %d\n", out->GetNI(i).GetInDeg(), out->GetNI(nodeId).GetInDeg(), out->GetNI(NbOutId).GetInDeg());
					//system("pause");
					break;
				}
			}

		}
	}



	for (int i = 0; i < nodesCount; i++){
		int InDeg = out->GetNI(i).GetInDeg();
		if (InDeg > MaxDeg){
			//int EToChange = rnd.GetUniDev() * (MaxDeg - MinDeg) + InDeg - MaxDeg;
			int EToChange = InDeg - MaxDeg;
			//if (EToChange > InDeg - MinDeg) EToChange = InDeg - MinDeg;
			//printf("EToChange = %d InDeg - MaxDeg = %d\n", EToChange, InDeg - MaxDeg);
			
			auto CurrNode = out->GetNI(i);
			//printf("Current node degree %d\n", CurrNode.GetInDeg());
			// rewire EToChange edges
			for (int j = InDeg-1; j >=0 ; j--){
				// rewire nodes from the end of list of edges (check!)
				int NodeToRewire = CurrNode.GetNbrNId(j);
				while (1){
					// get neighbour node using probability matrix
					//int NbNode = GetCol(RowProbCumV, NodeToRewire, NIter, rnd);
					// get random neighbour node
					int NbNode = rnd.GetUniDev() * (nodesCount-1);
					if (NbNode == i || out->IsEdge(NodeToRewire,NbNode)) continue;
					// if neighbour node degree is less than MaxDeg, rewire NodeToRewire
					if (out->GetNI(NbNode).GetInDeg() < MaxDeg){
						out->AddEdge(NodeToRewire, NbNode);
						out->AddEdge(NbNode,NodeToRewire);
						break;
					}
				}
				out->DelEdge(i, NodeToRewire, false);
				//printf("%d Edge was deleted\n", EToChange);
				if (out->GetNI(i).GetInDeg() < MinDeg){
					printf("Current node degree is less than min deg: %d < %d InDeg = %d EToChange=%d\n", CurrNode.GetInDeg(), MinDeg, InDeg, EToChange);
					system("pause");
				}
				if (out->GetNI(NodeToRewire).GetInDeg() < MinDeg){
					printf("Node to rewire degree is less than min deg:%d < %d\n", out->GetNI(NodeToRewire).GetInDeg(), MinDeg);
					system("pause");
				}
				EToChange--;
				if (EToChange == 0) break;
			}
		}
	}
	/*bool IsMaxDeg = false;
	int MaxDegNode = 0;
	
	for (int i = 0; i < EdgesToDel; i++){
		int nodeId = rnd.GetUniDev() * (nodesCount - 1);
		auto NI = out->GetNI(nodeId);
		int OutDeg = NI.GetOutDeg();
		double Val1 = rnd.GetUniDev(), Val2 = (OutDeg - MinDeg) / 99.0;
		if (Val2 < Val1 || OutDeg == 99 ) {i--; continue;}
		int NbId = static_cast<int>(rnd.GetUniDev() * (OutDeg - 1));
		int NbOutId = NI.GetNbrNId(NbId);
		if (rnd.GetUniDev() > (out->GetNI(NbOutId).GetInDeg() - MinDeg) / 99.0 || out->GetNI(NbOutId).GetInDeg() == 99 ) 
		{i--; continue;}

		out->DelEdge(nodeId, NbOutId, false);
	}
	*/
	/*for (int i = 0; i < nodesCount; i++){
		DegAfter.Add(out->GetNI(i).GetInDeg());
	}
	FILE *f = fopen("DegTest.dat", "w");
	fprintf(f, "Node\tDegree before\tDegree after\n");
	for (int i = 0; i < nodesCount; i++){
		fprintf(f, "%d\t%d\t%d\n", i, DegBefore[i], DegAfter[i]);
	}
	fclose(f);*/
	
}

PNGraph TKronMtx::GenFastKronecker(const TKronMtx& SeedMtx, const int& NIter, const bool& IsDir, const int& Seed, const TIntPr& InDegR, const TIntPr& OutDegR, PNGraph& Graph, double ModelClustCf){
	const int NNodes = SeedMtx.GetNodes(NIter);
	const int NEdges = SeedMtx.GetEdges(NIter);
	printf("GenFastKronecker() with restrictions. Nodes %d, edges %d\n", NNodes, NEdges);
	Graph = TNGraph::New(NNodes, -1);
	// add nodes
	for (int i = 0; i < NNodes; i++) {
		Graph->AddNode(i); }
	const TInt& InMin = InDegR.Val1, & OutMin = OutDegR.Val1;
	bool IsOut = InMin < OutMin ? true : false;
	IsOut = false;
	TRnd Rnd(Seed);
	int Collisions = 0;
	if (IsDir){
		Collisions += AddFirstDir(IsOut, InDegR, OutDegR, Graph, SeedMtx, NIter, Rnd);
		Collisions += AddSecondDir(!IsOut, InDegR, OutDegR, Graph, SeedMtx, NIter, Rnd);
	}
	else {
		Collisions += AddUnDir(InDegR, Graph, SeedMtx, NIter, Rnd, ModelClustCf);
	}
	const int Least = NEdges - Graph->GetEdges();
	Collisions += AddEdges(SeedMtx, NIter, IsDir, Rnd, Graph, Least, InDegR.Val2, OutDegR.Val2, ModelClustCf);
	printf("             collisions: %d (%.4f)\n", Collisions, Collisions/(double)Graph->GetEdges());
	//RemoveZeroDegreeNodes(Graph, SeedMtx, NIter, InDegR.Val1, InDegR.Val2);
	return Graph;
}

int TKronMtx::CheckClustCf(const PNGraph& Graph, int Row, int Col, double ModelClustCf, TRnd& Rnd, const bool& IsDir, int& Collision)
{
	int NeighboursCount = Graph->GetNI(Row).GetOutDeg();
	int Triangles = NeighboursCount * (NeighboursCount - 1);
	if (NeighboursCount  > 1){
		const double& ProbCF = Rnd.GetUniDev();
		int NewRow, NewCol;
		if (ProbCF <= ModelClustCf){
			/*int NeighID = floor(Rnd.GetUniDev() * NeighboursCount);
			NewRow = Graph->GetNI(Row).GetNbrNId(NeighID);*/
			NewRow = Col;
			int NeighID = floor(Rnd.GetUniDev() * NeighboursCount);
			NewCol = Graph->GetNI(Row).GetNbrNId(NeighID); 
			//printf("%d %d\n", NewRow, NewCol);
			//ClosedTriads++;
			if (! Graph->IsEdge(NewRow, NewCol) && NewRow != NewCol ){
				//printf("Row %d col %d newrow %d newcol %d\n", Row, Col, NewRow, NewCol);
				Graph->AddEdge(NewRow, NewCol);  
				if (! IsDir) 
					Graph->AddEdge(NewCol, NewRow);
				return 1;
			}
			else Collision++;
			return 0;
		}
		return 0;
	}
	return 0;
}

bool CheckMu(double T1, double T2, double T3, double T4, double Mu){
	if (Mu > (T1 + T4) / 2) 
		return false;
	if (Mu < (T4 - 1) * (T1 + T4) / (2 * T4)) 
		return false;
	if (Mu < (T1 - 1) * (T1 + T4) / (2 * T1)) 
		return false;
	if (Mu > 1 - T2) 
		return false;
	if (Mu > 1 - T3) 
		return false;
	return true;
}

  double TKronMtx::GetMinPossibleDeg(){
	  double Val = GetMtxSum() / 4; 
	  return (2 * Val);
  }
  double TKronMtx::GetMaxPossibleDeg(){
	  double Val = (GetMtxSum() - 1) / 2;
	  return (1 + Val);
  }

// works for 2x2 matrix
double TKronMtx::GetNoisedProbV(TVec<TVec<TFltIntIntTr>>&ProbToRCPosV, const TFlt& NoiseCoeff, TRnd& Rnd, const int& NIter, const TKronMtx& SeedMtx){
	double T1 = SeedMtx.At(0,0), T2 = SeedMtx.At(0,1), T3 = SeedMtx.At(1,0), T4 = SeedMtx.At(1,1);
	TKronMtx BaseMtx(SeedMtx);
	double BaseExpectedDeg = BaseMtx.GetMaxExpectedDeg(),
		MinPossibleDeg = BaseMtx.GetMinPossibleDeg(), MaxPossibleDeg = BaseMtx.GetMaxPossibleDeg(),
		UpperLimit = MaxPossibleDeg - BaseExpectedDeg, LowerLimit =  BaseExpectedDeg - MinPossibleDeg;
	double Limit = (UpperLimit < LowerLimit) ? UpperLimit : LowerLimit;
	//double Limit = MaxPossibleDeg / 2;
	double MtxSum = T1 + T2 + T3 + T4;
	double B = (T1 + T4) / 2 < T2 ? (T1 + T4) / 2 : T2;
	B *= NoiseCoeff.Val;
	Limit *= NoiseCoeff.Val;
	double AvgExpectedDeg = 1;
	double PrevDeg = 0, CurrentDeg = 0;
	double AccumExpected = 1, Step = BaseExpectedDeg, AccumReal = 1;

	bool CanBeNegative = true;

	for (int i = 0; i < NIter; i++){
		TKronMtx NewMtx(BaseMtx);
		// if i == NIter-1 and i is even, leave the matrix without change
		//if (i % 2 == 0 && i != NIter-1){
			//[Base - Limit; Base + Limit]
			//CurrentDeg = BaseExpectedDeg - Limit + Rnd.GetUniDev() * 2 * Limit;
			//if (CurrentDeg < BaseExpectedDeg) CurrentDeg *= 1.05;
			//NewMtx.SetForMaxDeg(CurrentDeg, 1);
			double Mu = 0; 
			do{
				Mu = B * (-1) + Rnd.GetUniDev() * (2 * B);
				if (AccumExpected < AccumReal) Mu *= -1;
				//printf("Mu = %f\n", Mu);
			}
			while (!CheckMu(T1, T2, T3, T4, Mu));
		
			if (AccumExpected > AccumReal) {
				double M = MaxPossibleDeg;
				NewMtx.SetForMaxDeg(MaxPossibleDeg , 1); 
				CurrentDeg = MaxPossibleDeg;
			}
			else {
				NewMtx.At(0,0) -= (2 * Mu * T1) / (T1 + T4);
				NewMtx.At(0,1) += Mu;
				NewMtx.At(1,0) += Mu;
				NewMtx.At(1,1) -= (2 * Mu * T4) / (T1 + T4);
				CurrentDeg = NewMtx.GetMaxExpectedDeg();
			}
			if (i == NIter - 1){
				double RequiredDeg = AccumExpected * Step / AccumReal;
				if (RequiredDeg > MaxPossibleDeg){
					RequiredDeg = MaxPossibleDeg;
				}
				else if (RequiredDeg < MinPossibleDeg)
					RequiredDeg = MinPossibleDeg;
				NewMtx.SetForMaxDeg(RequiredDeg , 1); 
				CurrentDeg = RequiredDeg;

			}
			AccumReal *= CurrentDeg;
			AccumExpected *= Step;
			printf("CurrentDeg = %f, AccumExpected = %f, AccumReal = %f\n", CurrentDeg, AccumExpected, AccumReal);
			if (AccumReal > AccumExpected) CanBeNegative = true;
			else CanBeNegative = false;

		//}
		//else {
		//	CurrentDeg = BaseExpectedDeg / (CurrentDeg / BaseExpectedDeg );
			//if (CurrentDeg < BaseExpectedDeg) CurrentDeg *= 1.05;
		//	NewMtx.SetForMaxDeg(CurrentDeg, 1);
		//}
		
		/*FILE *f = fopen("NoiseMtx.tab", "a");
		NewMtx.Dump(f);
		AvgExpectedDeg *= CurrentDeg;
		fprintf(f, "Maximum expected degree: %f\n", CurrentDeg);
		fclose(f);*/

		TVec<TFltIntIntTr> MtxVec;
		double CumProb = 0.0;
		for (int r = 0; r < NewMtx.GetDim(); r++) {
			for (int c = 0; c < NewMtx.GetDim(); c++) {
				double Prob = NewMtx.At(r, c);
				if (Prob > 0.0) {
					CumProb += Prob;
					MtxVec.Add(TFltIntIntTr(CumProb/MtxSum, r, c));
				}
			}
		}
		
		ProbToRCPosV.Add(MtxVec);
	}
	printf("AccumExpected = %f, AccumReal = %f\n", AccumExpected, AccumReal);
	/*FILE *f = fopen("NoiseMtx.tab", "a");
	fprintf(f, "------------------------------------------------");
	fprintf(f, "Average expected deg: %f\n", AvgExpectedDeg);
	fclose(f);*/
	return AvgExpectedDeg;
}

void TKronMtx::GetLemma3Estimates(const TKronMtx& SeedMtx, const int& NIter, const int& MaxDeg){
	vector<vector<double>> Data;
	vector<double> Deg; vector<double> NC;
	for (int i = 0; i < static_cast<int>(MaxDeg + 0.5); i++){
		int NodesCount = TKronMtx::GetExpectedNodesCount(SeedMtx, NIter, i);
		Deg.push_back(i);
		NC.push_back(NodesCount);
	}
	Data.push_back(Deg); Data.push_back(NC);
	TStrV RowNames; RowNames.Add("Degree"); RowNames.Add("Nodes count");
	MakeDatFile("Lemma3", "Expected degrees", RowNames, Data, SeedMtx.GetNodes(NIter), SeedMtx.GetEdges(NIter));
}

// get E(X_d) for degree Deg
int TKronMtx::GetExpectedNodesCount(const TKronMtx& Mtx, const int& NIter, const int& Deg){
	if (Deg == 0) 
		return 0;
	double ExpDeg = 0;
	double Dim = static_cast<double>(Mtx.GetDim());
	if (Dim != 2)
		Error("GetExpectedNodesCount", "Function works only with 2x2 matrix");
	int NodesCount = pow (Dim, NIter);
	double MtxSum = Mtx.GetMtxSum();
	double EdgesCount = pow (MtxSum, NIter);
	double DeltaBig = EdgesCount / NodesCount * 2;
	// in the paper T1 is the largest value
	int Row = 0, Col = 0;
	double T1 = 0.0;
	for (size_t i = 0; i < Dim; i++)
		for (size_t j = 0; j < Dim; j++){
			if (Mtx.At(i,j) > T1){
				T1 = Mtx.At(i,j);
				Row = i;
				Col = j;
			}
		}
	Col = (Col == 0) ? 1 : 0;
	double T2 = Mtx.At(Row, Col);
	Row = (Row == 0) ? 1 : 0;
	double T4 = Mtx.At(Row, Col);
	Col = (Col == 0) ? 1 : 0;
	double T3 = Mtx.At(Row, Col);
	// in the paper Delta = (T1 + T2) - 0.5, but MtxSum is equal to 1 in that case (use 1 instead of MtxSum)
	double Delta = (T1 + T2) - MtxSum / 2;
	double Tau = (MtxSum + 2 * Delta) / (MtxSum - 2 * Delta);
	// floor() is my assumption
	double Lambda = DeltaBig * pow(1 - 4 * Delta * Delta, NIter / 2.00);
	double ThetaD = log(Deg / Lambda) / log(Tau);
	int RD = floor(ThetaD);
	double DeltaD = ThetaD - RD;
	if (Deg < exp(1.0) * log(2.0) * NIter || Deg > pow (NodesCount, 0.5))
		printf("Warning. Deg < e ln(2) * l or Deg > sqrt(N)\n");
	// floor() is my assumption
 	if (RD >= floor(NIter / 2.00)) return 0;
	else {
		// floor() is my assumption
		ExpDeg = 1 / pow(2 * M_PI * Deg, 0.5) * exp(-1 * Deg * DeltaD * DeltaD * log(Tau) * log(Tau) / 2) * GetBinomCoeff(NIter, floor(NIter / 2.00) + RD);
		// floor() is my assumption
		ExpDeg += 1 / pow(2 * M_PI * Deg, 0.5) * exp(-1 * Deg * (1 - DeltaD) * (1 - DeltaD) * log(Tau) * log(Tau) / 2) * GetBinomCoeff(NIter, floor(NIter / 2.00) + RD + 1);
	}
	return static_cast<int>(ExpDeg + 0.5);
}

void ProbToRCPosVDump(const TVec<TVec<TFltIntIntTr>>& Prob, const int& MtxDim, const int& NIter){
	for (int i = 0; i < NIter; i++){
		for (int j = 0; j < MtxDim * MtxDim; j++){
			printf("%f ", Prob[i][j].Val1);
		}
		printf("\n");
	}
}

// use RMat like recursive descent to quickly generate a Kronecker graph
PNGraph TKronMtx::GenFastKronecker(const TKronMtx& SeedMtx, const int& NIter, const bool& IsDir, double &AvgExpectedDeg, const int& Seed, double NoiseCoeff) {
  const TKronMtx& SeedGraph = SeedMtx;
  const int MtxDim = SeedGraph.GetDim();
  const double MtxSum = SeedGraph.GetMtxSum();
  const int NNodes = SeedGraph.GetNodes(NIter);
  const int NEdges = SeedGraph.GetEdges(NIter);
  //const double DiagEdges = NNodes * pow(SeedGraph.At(0,0), double(NIter));
  //const int NEdges = (int) TMath::Round(((pow(double(SeedGraph.GetMtxSum()), double(NIter)) - DiagEdges) /2.0));
  printf("  FastKronecker: %d nodes, %d edges, %s...\n", NNodes, NEdges, IsDir ? "Directed":"UnDirected");
  PNGraph Graph = TNGraph::New(NNodes, -1);
  TRnd Rnd(Seed);
  TExeTm ExeTm;
  TVec<TVec<TFltIntIntTr>> ProbToRCPosV;
  AvgExpectedDeg = GetNoisedProbV(ProbToRCPosV, NoiseCoeff, Rnd, NIter, SeedMtx);
  
  //ProbToRCPosVDump(ProbToRCPosV, MtxDim, NIter);

  // add nodes
  for (int i = 0; i < NNodes; i++) {
    Graph->AddNode(i); }
  int ClosedTriads = 0;
  int ClustCollision = 0;
  // add edges
  int Rng, Row, Col, Collision=0, n = 0;
  
  for (int edges = 0; edges < NEdges ; ) {
    Rng=NNodes;  Row=0;  Col=0;
    for (int iter = 0; iter < NIter; iter++) {
	  const double& Prob = Rnd.GetUniDev();
      n = 0; while(Prob > ProbToRCPosV[iter][n].Val1) { n++; }
      const int MtxRow = ProbToRCPosV[iter][n].Val2;
      const int MtxCol = ProbToRCPosV[iter][n].Val3;
      Rng /= MtxDim;
      Row += MtxRow * Rng;
      Col += MtxCol * Rng;
    }
    if (! Graph->IsEdge(Row, Col) && Row != Col ) { // allow self-loops
      Graph->AddEdge(Row, Col);  edges++;
	  if (! IsDir) {
        Graph->AddEdge(Col, Row);
        edges++;
      }
    } else { Collision++; }
    //if (edges % 1000 == 0) printf("\r...%dk", edges/1000);
  }

  printf("             %d edges [%s]\n", Graph->GetEdges(), ExeTm.GetTmStr());
  printf("             collisions: %d (%.4f)\n", Collision, Collision/(double)Graph->GetEdges());
  printf("ClosedTriads %d ClustCollision %d\n", ClosedTriads, ClustCollision);
  return Graph;
}

// use RMat like recursive descent to quickly generate a Kronecker graph
PNGraph TKronMtx::GenFastKronecker(const TKronMtx& SeedMtx, const int& NIter, const int& Edges, const bool& IsDir, const int& Seed) {
  const TKronMtx& SeedGraph = SeedMtx;
  const int MtxDim = SeedGraph.GetDim();
  const double MtxSum = SeedGraph.GetMtxSum();
  const int NNodes = SeedGraph.GetNodes(NIter);
  const int NEdges = Edges;
  //const double DiagEdges = NNodes * pow(SeedGraph.At(0,0), double(NIter));
  //const int NEdges = (int) TMath::Round(((pow(double(SeedGraph.GetMtxSum()), double(NIter)) - DiagEdges) /2.0));
  printf("  RMat Kronecker: %d nodes, %d edges, %s...\n", NNodes, NEdges, IsDir ? "Directed":"UnDirected");
  PNGraph Graph = TNGraph::New(NNodes, -1);
  TRnd Rnd(Seed);
  TExeTm ExeTm;
  // prepare cell probability vector
  TVec<TFltIntIntTr> ProbToRCPosV; // row, col position
  double CumProb = 0.0;
  for (int r = 0; r < MtxDim; r++) {
    for (int c = 0; c < MtxDim; c++) {
      const double Prob = SeedGraph.At(r, c);
      if (Prob > 0.0) {
        CumProb += Prob;
        ProbToRCPosV.Add(TFltIntIntTr(CumProb/MtxSum, r, c));
      }
    }
  }
  // add nodes
  for (int i = 0; i < NNodes; i++) {
    Graph->AddNode(i); }
  // add edges
  int Rng, Row, Col, Collision=0, n = 0;
  for (int edges = 0; edges < NEdges; ) {
    Rng=NNodes;  Row=0;  Col=0;
    for (int iter = 0; iter < NIter; iter++) {
      const double& Prob = Rnd.GetUniDev();
      n = 0; while(Prob > ProbToRCPosV[n].Val1) { n++; }
      const int MtxRow = ProbToRCPosV[n].Val2;
      const int MtxCol = ProbToRCPosV[n].Val3;
      Rng /= MtxDim;
      Row += MtxRow * Rng;
      Col += MtxCol * Rng;
    }
    if (! Graph->IsEdge(Row, Col)) { // allow self-loops
      Graph->AddEdge(Row, Col);  edges++;
      if (! IsDir) {
        if (Row != Col) Graph->AddEdge(Col, Row);
        edges++;
      }
    } else { Collision++; }
    //if (edges % 1000 == 0) printf("\r...%dk", edges/1000);
  }
  //printf("             %d edges [%s]\n", Graph->GetEdges(), ExeTm.GetTmStr());
  printf("             collisions: %d (%.4f)\n", Collision, Collision/(double)Graph->GetEdges());
  return Graph;
}

PNGraph TKronMtx::GenDetKronecker(const TKronMtx& SeedMtx, const int& NIter, const bool& IsDir) {
  const TKronMtx& SeedGraph = SeedMtx;
  const int NNodes = SeedGraph.GetNodes(NIter);
  printf("  Deterministic Kronecker: %d nodes, %s...\n", NNodes, IsDir ? "Directed":"UnDirected");
  PNGraph Graph = TNGraph::New(NNodes, -1);
  TExeTm ExeTm;
  int edges = 0;
  for (int node1 = 0; node1 < NNodes; node1++) { Graph->AddNode(node1); }

  for (int node1 = 0; node1 < NNodes; node1++) {
    for (int node2 = 0; node2 < NNodes; node2++) {
      if (SeedGraph.IsEdgePlace(node1, node2, NIter, Rnd.GetUniDev())) {
        Graph->AddEdge(node1, node2);
        edges++;
      }
    }
    if (node1 % 1000 == 0) printf("\r...%dk, %dk", node1/1000, edges/1000);
  }
  return Graph;
}

void TKronMtx::PlotCmpGraphs(const TKronMtx& SeedMtx, const PNGraph& Graph, const TStr& FNmPref, const TStr& Desc) {
  const int KronIters = SeedMtx.GetKronIter(Graph->GetNodes());
  PNGraph KronG, WccG;
  const bool FastGen = true;
  if (FastGen) { KronG = TKronMtx::GenFastKronecker(SeedMtx, KronIters, true, 0); }
  else { KronG = TKronMtx::GenKronecker(SeedMtx, KronIters, true, 0); }
  TSnap::DelZeroDegNodes(KronG);
  WccG = TSnap::GetMxWcc(KronG);
  const TStr Desc1 = TStr::Fmt("%s", Desc.CStr());
  TGStatVec GS(tmuNodes, TFSet() | gsdInDeg | gsdOutDeg | gsdWcc | gsdHops | gsdScc | gsdClustCf | gsdSngVec | gsdSngVal);
  //gsdHops
  //gsWccHops, gsdSngVal, gsdSngVec, gsdClustCf
  GS.Add(Graph, TSecTm(1), TStr::Fmt("GRAPH  G(%d, %d)", Graph->GetNodes(), Graph->GetEdges()));
  GS.Add(KronG, TSecTm(2), TStr::Fmt("KRONECKER  K(%d, %d)", KronG->GetNodes(), KronG->GetEdges()));
  GS.Add(WccG, TSecTm(3),  TStr::Fmt("KRONECKER  wccK(%d, %d)", WccG->GetNodes(), WccG->GetEdges()));
  const TStr Style = "linewidth 1 pointtype 6 pointsize 1";
  GS.ImposeDistr(gsdInDeg, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdInDeg, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdOutDeg, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdOutDeg, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdHops, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdClustCf, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdClustCf, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVal, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVal, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVec, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVec, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdWcc, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdWcc, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdScc, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdScc, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
//    typedef enum { distrUndef, distrInDeg, distrOutDeg, distrWcc, distrScc,
//  distrHops, distrWccHops, distrSngVal, distrSngVec, distrClustCf, distrMx } TGraphDistr;*/
}

void TKronMtx::PlotCmpGraphs(const TKronMtx& SeedMtx1, const TKronMtx& SeedMtx2, const PNGraph& Graph, const TStr& FNmPref, const TStr& Desc) {
  const int KronIters1 = SeedMtx1.GetKronIter(Graph->GetNodes());
  const int KronIters2 = SeedMtx2.GetKronIter(Graph->GetNodes());
  PNGraph KronG1, KronG2;
  const bool FastGen = true;
  if (FastGen) {
    KronG1 = TKronMtx::GenFastKronecker(SeedMtx1, KronIters1, true, 0);
    KronG2 = TKronMtx::GenFastKronecker(SeedMtx2, KronIters2, false, 0); } //!!!
  else {
    KronG1 = TKronMtx::GenKronecker(SeedMtx1, KronIters1, true, 0);
    KronG2 = TKronMtx::GenKronecker(SeedMtx2, KronIters2, true, 0);  }
  TSnap::DelZeroDegNodes(KronG1);
  TSnap::DelZeroDegNodes(KronG2);
  const TStr Desc1 = TStr::Fmt("%s", Desc.CStr());
  TGStatVec GS(tmuNodes, TFSet() | gsdInDeg | gsdOutDeg | gsdWcc | gsdScc | gsdHops | gsdClustCf | gsdSngVec | gsdSngVal | gsdTriadPart);
  //gsdHops
  //gsWccHops, gsdSngVal, gsdSngVec, gsdClustCf
  GS.Add(Graph, TSecTm(1), TStr::Fmt("GRAPH  G(%d, %d)", Graph->GetNodes(), Graph->GetEdges()));
  GS.Add(KronG1, TSecTm(2), TStr::Fmt("KRONECKER1  K(%d, %d) %s", KronG1->GetNodes(), KronG1->GetEdges(), SeedMtx1.GetMtxStr().CStr()));
  GS.Add(KronG2, TSecTm(3),  TStr::Fmt("KRONECKER2  K(%d, %d) %s", KronG2->GetNodes(), KronG2->GetEdges(), SeedMtx2.GetMtxStr().CStr()));
  const TStr Style = "linewidth 1 pointtype 6 pointsize 1";
  // raw data
  GS.ImposeDistr(gsdInDeg, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdOutDeg, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdHops, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdClustCf, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVal, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVec, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdWcc, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdScc, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdTriadPart, FNmPref, Desc1, false, false, gpwLinesPoints, Style);
  // smooth
  GS.ImposeDistr(gsdInDeg, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdOutDeg, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdClustCf, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdScc, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdWcc, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVec, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdSngVal, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
  GS.ImposeDistr(gsdTriadPart, FNmPref+"-B", Desc1, true, false, gpwLinesPoints, Style);
}

void TKronMtx::PlotCmpGraphs(const TVec<TKronMtx>& SeedMtxV, const PNGraph& Graph, const TStr& FNmPref, const TStr& Desc) {
  const TStr Desc1 = TStr::Fmt("%s", Desc.CStr());
  TGStatVec GS(tmuNodes, TFSet() | gsdInDeg | gsdOutDeg | gsdWcc | gsdScc | gsdHops | gsdClustCf | gsdSngVec | gsdSngVal);
  GS.Add(Graph, TSecTm(1), TStr::Fmt("GRAPH  G(%d, %d)", Graph->GetNodes(), Graph->GetEdges()));
  //gsdHops
  //gsWccHops, gsdSngVal, gsdSngVec, gsdClustCf
  for (int m = 0; m < SeedMtxV.Len(); m++) {
    const int KronIters = SeedMtxV[m].GetKronIter(Graph->GetNodes());
    PNGraph KronG1 = TKronMtx::GenFastKronecker(SeedMtxV[m], KronIters, true, 0);
    printf("*** K(%d, %d) n0=%d\n", KronG1->GetNodes(), KronG1->GetEdges(), SeedMtxV[m].GetDim());
    TSnap::DelZeroDegNodes(KronG1);
    printf(" del zero deg K(%d, %d) n0=%d\n", KronG1->GetNodes(), KronG1->GetEdges(), m);
    GS.Add(KronG1, TSecTm(m+2), TStr::Fmt("K(%d, %d) n0^k=%d n0=%d", KronG1->GetNodes(), KronG1->GetEdges(), SeedMtxV[m].GetNZeroK(Graph), SeedMtxV[m].GetDim()));
    // plot after each Kronecker is done
    const TStr Style = "linewidth 1 pointtype 6 pointsize 1";
    GS.ImposeDistr(gsdInDeg, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdInDeg, FNmPref+"-B", Desc1, true, false, gpwLines, Style);
    GS.ImposeDistr(gsdOutDeg, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdOutDeg, FNmPref+"-B", Desc1, true, false, gpwLines, Style);
    GS.ImposeDistr(gsdHops, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdClustCf, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdClustCf, FNmPref+"-B", Desc1, true, false, gpwLines, Style);
    GS.ImposeDistr(gsdSngVal, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdSngVal, FNmPref+"-B", Desc1, true, false, gpwLines, Style);
    GS.ImposeDistr(gsdSngVec, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdSngVec, FNmPref+"-B", Desc1, true, false, gpwLines, Style);
    GS.ImposeDistr(gsdWcc, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdWcc, FNmPref+"-B", Desc1, true, false, gpwLines, Style);
    GS.ImposeDistr(gsdScc, FNmPref, Desc1, false, false, gpwLines, Style);
    GS.ImposeDistr(gsdScc, FNmPref+"-B", Desc1, true, false, gpwLines, Style);
  }
  //    typedef enum { distrUndef, distrInDeg, distrOutDeg, distrWcc, distrScc,
  //  distrHops, distrWccHops, distrSngVal, distrSngVec, distrClustCf, distrMx } TGraphDistr;*/
}

void TKronMtx::KronMul(const TKronMtx& Left, const TKronMtx& Right, TKronMtx& Result) {
  const int LDim = Left.GetDim();
  const int RDim = Right.GetDim();
  Result.GenMtx(LDim * RDim);
  for (int r1 = 0; r1 < LDim; r1++) {
    for (int c1 = 0; c1 < LDim; c1++) {
      const double& Val = Left.At(r1, c1);
      for (int r2 = 0; r2 < RDim; r2++) {
        for (int c2 = 0; c2 < RDim; c2++) {
          Result.At(r1*RDim+r2, c1*RDim+c2) = Val * Right.At(r2, c2);
        }
      }
    }
  }
}

void TKronMtx::KronSum(const TKronMtx& Left, const TKronMtx& Right, TKronMtx& Result) {
  const int LDim = Left.GetDim();
  const int RDim = Right.GetDim();
  Result.GenMtx(LDim * RDim);
  for (int r1 = 0; r1 < LDim; r1++) {
    for (int c1 = 0; c1 < LDim; c1++) {
      const double& Val = Left.At(r1, c1);
      for (int r2 = 0; r2 < RDim; r2++) {
        for (int c2 = 0; c2 < RDim; c2++) {
          if (Val == NInf || Right.At(r2, c2) == NInf) {
            Result.At(r1*RDim+r2, c1*RDim+c2) = NInf; }
          else {
            Result.At(r1*RDim+r2, c1*RDim+c2) = Val + Right.At(r2, c2); }
        }
      }
    }
  }
}

void TKronMtx::KronPwr(const TKronMtx& KronMtx, const int& NIter, TKronMtx& OutMtx) {
  OutMtx = KronMtx;
  TKronMtx NewOutMtx;
  for (int iter = 0; iter < NIter; iter++) {
    KronMul(OutMtx, KronMtx, NewOutMtx);
    NewOutMtx.Swap(OutMtx);
  }

}

void TKronMtx::Dump(ofstream& TFile) const{
	double Sum=0.0;
	TFltV ValV = SeedMtx;
	for (int i = 0; i < ValV.Len(); i++) {
		TFile << ValV[i]() << " ";
		Sum += ValV[i];
		if ((i+1) % GetDim() == 0) TFile << endl; 
	}
	TFile << " (sum: " << Sum << ")" << endl;
}

void TKronMtx::Dump(FILE* TFile) const{
	double Sum=0.0;
	TFltV ValV = SeedMtx;
	for (int i = 0; i < ValV.Len(); i++) {
		fprintf(TFile, "%f ",ValV[i]());
		Sum += ValV[i];
		if ((i+1) % GetDim() == 0) fprintf(TFile, "\n"); 
	}
	fprintf(TFile, " (sum: %f)\n",Sum);
}

void TKronMtx::Dump(const TStr& MtxNm, const bool& Sort) const {
  /*printf("%s: %d x %d\n", MtxNm.Empty()?"Mtx":MtxNm.CStr(), GetDim(), GetDim());
  for (int r = 0; r < GetDim(); r++) {
    for (int c = 0; c < GetDim(); c++) { printf("  %8.2g", At(r, c)); }
    printf("\n");
  }*/
  if (! MtxNm.Empty()) printf("%s\n", MtxNm.CStr());
  double Sum=0.0;
  TFltV ValV = SeedMtx;
  if (Sort) { ValV.Sort(false); }
  for (int i = 0; i < ValV.Len(); i++) {
    printf("  %10.4g", ValV[i]());
    Sum += ValV[i];
    if ((i+1) % GetDim() == 0) { printf("\n"); }
  }
  printf(" (sum:%.4f)\n", Sum);
}

// average difference in the parameters
double TKronMtx::GetAvgAbsErr(const TKronMtx& Kron1, const TKronMtx& Kron2) {
  TFltV P1 = Kron1.GetMtx();
  TFltV P2 = Kron2.GetMtx();
  IAssert(P1.Len() == P2.Len());
  P1.Sort();  P2.Sort();
  double delta = 0.0;
  for (int i = 0; i < P1.Len(); i++) {
    delta += fabs(P1[i] - P2[i]);
  }
  return delta/P1.Len();
}

// average L2 difference in the parameters
double TKronMtx::GetAvgFroErr(const TKronMtx& Kron1, const TKronMtx& Kron2) {
  TFltV P1 = Kron1.GetMtx();
  TFltV P2 = Kron2.GetMtx();
  IAssert(P1.Len() == P2.Len());
  P1.Sort();  P2.Sort();
  double delta = 0.0;
  for (int i = 0; i < P1.Len(); i++) {
    delta += pow(P1[i] - P2[i], 2);
  }
  return sqrt(delta/P1.Len());
}

// get matrix from matlab matrix notation
TKronMtx TKronMtx::GetMtx(TStr MatlabMtxStr) {
  TStrV RowStrV, ColStrV;
  MatlabMtxStr.ChangeChAll(',', ' ');
  MatlabMtxStr.SplitOnAllCh(';', RowStrV);  IAssert(! RowStrV.Empty());
  RowStrV[0].SplitOnWs(ColStrV);    IAssert(! ColStrV.Empty());
  const int Rows = RowStrV.Len();
  const int Cols = ColStrV.Len();
  IAssert(Rows == Cols);
  TKronMtx Mtx(Rows);
  for (int r = 0; r < Rows; r++) {
    RowStrV[r].SplitOnWs(ColStrV);
    IAssert(ColStrV.Len() == Cols);
    for (int c = 0; c < Cols; c++) {
      Mtx.At(r, c) = (double) ColStrV[c].GetFlt(); }
  }
  return Mtx;
}

TKronMtx TKronMtx::GetRndMtx(const int& Dim, const double& MinProb) {
  TKronMtx Mtx;
  Mtx.SetRndMtx(Dim, MinProb);
  return Mtx;
}

TKronMtx TKronMtx::GetInitMtx(const int& Dim, const int& Nodes, const int& Edges) {
  const double MxParam = 0.8+TKronMtx::Rnd.GetUniDev()/5.0;
  const double MnParam = 0.2-TKronMtx::Rnd.GetUniDev()/5.0;
  const double Step = (MxParam-MnParam) / (Dim*Dim-1);
  TFltV ParamV(Dim*Dim);
  if (Dim == 1) { ParamV.PutAll(0.5); } // random graph
  else {
    for (int p = 0; p < ParamV.Len(); p++) {
      ParamV[p] = MxParam - p*Step; }
  }
  //IAssert(ParamV[0]==MxParam && ParamV.Last()==MnParam);
  TKronMtx Mtx(ParamV);
  Mtx.SetForEdges(Nodes, Edges);
  return Mtx;
}

TKronMtx TKronMtx::GetInitMtx(const TStr& MtxStr, const int& Dim, const int& Nodes, const int& Edges) {
  TKronMtx Mtx(Dim);
  if (TCh::IsNum(MtxStr[0])) { Mtx = TKronMtx::GetMtx(MtxStr); }
  else if (MtxStr[0] == 'r') { Mtx = TKronMtx::GetRndMtx(Dim, 0.1); }
  else if (MtxStr[0] == 'a') {
    const double Prob = TKronMtx::Rnd.GetUniDev();
    if (Prob < 0.4) {
      Mtx = TKronMtx::GetInitMtx(Dim, Nodes, Edges); }
    else { // interpolate so that there are in the corners 0.9, 0.5, 0.1, 0.5
      const double Max = 0.9+TKronMtx::Rnd.GetUniDev()/10.0;
      const double Min = 0.1-TKronMtx::Rnd.GetUniDev()/10.0;
      const double Med = (Max-Min)/2.0;
      Mtx.At(0,0)      = Max;       Mtx.At(0,Dim-1) = Med;
      Mtx.At(Dim-1, 0) = Med;  Mtx.At(Dim-1, Dim-1) = Min;
      for (int i = 1; i < Dim-1; i++) {
        Mtx.At(i,i) = Max - double(i)*(Max-Min)/double(Dim-1);
        Mtx.At(i, 0) = Mtx.At(0, i) = Max - double(i)*(Max-Med)/double(Dim-1);
        Mtx.At(i, Dim-1) = Mtx.At(Dim-1, i) = Med - double(i)*(Med-Min)/double(Dim-1);
      }
      for (int i = 1; i < Dim-1; i++) {
        for (int j = 1; j < Dim-1; j++) {
          if (i >= j) { continue; }
          Mtx.At(i,j) = Mtx.At(j,i) = Mtx.At(i,i) - (j-i)*(Mtx.At(i,i)-Mtx.At(i,Dim-1))/(Dim-i-1);
        }
      }
      Mtx.AddRndNoise(0.1);
    }
  } else { FailR("Wrong mtx: matlab str, or random (r), or all (a)"); }
  Mtx.SetForEdges(Nodes, Edges);
  return Mtx;
}

TKronMtx TKronMtx::GetMtxFromNm(const TStr& MtxNm) {
  if (MtxNm == "3chain") return TKronMtx::GetMtx("1 1 0; 1 1 1; 0 1 1");
  else if (MtxNm == "4star") return TKronMtx::GetMtx("1 1 1 1; 1 1 0 0 ; 1 0 1 0; 1 0 0 1");
  else if (MtxNm == "4chain") return TKronMtx::GetMtx("1 1 0 0; 1 1 1 0 ; 0 1 1 1; 0 0 1 1");
  else if (MtxNm == "4square") return TKronMtx::GetMtx("1 1 0 1; 1 1 1 0 ; 0 1 1 1; 1 0 1 1");
  else if (MtxNm == "5star") return TKronMtx::GetMtx("1 1 1 1 1; 1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1");
  else if (MtxNm == "6star") return TKronMtx::GetMtx("1 1 1 1 1 1; 1 1 0 0 0 0; 1 0 1 0 0 0; 1 0 0 1 0 0; 1 0 0 0 1 0; 1 0 0 0 0 1");
  else if (MtxNm == "7star") return TKronMtx::GetMtx("1 1 1 1 1 1 1; 1 1 0 0 0 0 0; 1 0 1 0 0 0 0; 1 0 0 1 0 0 0; 1 0 0 0 1 0 0; 1 0 0 0 0 1 0; 1 0 0 0 0 0 1");
  else if (MtxNm == "5burst") return TKronMtx::GetMtx("1 1 1 1 0; 1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 1; 0 0 0 1 1");
  else if (MtxNm == "7burst") return TKronMtx::GetMtx("1 0 0 1 0 0 0; 0 1 0 1 0 0 0; 0 0 1 1 0 0 0; 1 1 1 1 1 0 0; 0 0 0 1 1 1 1; 0 0 0 0 1 1 0; 0 0 0 0 1 0 1");
  else if (MtxNm == "7cross") return TKronMtx::GetMtx("1 0 0 1 0 0 0; 0 1 0 1 0 0 0; 0 0 1 1 0 0 0; 1 1 1 1 1 0 0; 0 0 0 1 1 1 0; 0 0 0 0 1 1 1; 0 0 0 0 0 1 1");
  FailR(TStr::Fmt("Unknow matrix: '%s'", MtxNm.CStr()).CStr());
  return TKronMtx();
}

TKronMtx TKronMtx::LoadTxt(const TStr& MtxFNm) {
  PSs Ss = TSs::LoadTxt(ssfTabSep, MtxFNm);
  IAssertR(Ss->GetXLen() == Ss->GetYLen(), "Not a square matrix");
  IAssert(Ss->GetYLen() == Ss->GetXLen());
  TKronMtx Mtx(Ss->GetYLen());
  for (int r = 0; r < Ss->GetYLen(); r++) {
    for (int c = 0; c < Ss->GetXLen(); c++) {
      Mtx.At(r, c) = (double) Ss->At(c, r).GetFlt(); }
  }
  return Mtx;
}


/////////////////////////////////////////////////
// Kronecker Log Likelihood
TKroneckerLL::TKroneckerLL(const PNGraph& GraphPt, const TFltV& ParamV, const double& PermPSwapNd): PermSwapNodeProb(PermPSwapNd) {
  InitLL(GraphPt, TKronMtx(ParamV));
}

TKroneckerLL::TKroneckerLL(const PNGraph& GraphPt, const TKronMtx& ParamMtx, const double& PermPSwapNd) : PermSwapNodeProb(PermPSwapNd) {
  InitLL(GraphPt, ParamMtx);
}

TKroneckerLL::TKroneckerLL(const PNGraph& GraphPt, const TKronMtx& ParamMtx, const TIntV& NodeIdPermV, const double& PermPSwapNd) : PermSwapNodeProb(PermPSwapNd) {
  InitLL(GraphPt, ParamMtx);
  NodePerm = NodeIdPermV;
  SetIPerm(NodePerm);
}

PKroneckerLL TKroneckerLL::New(const PNGraph& GraphPt, const TKronMtx& ParamMtx, const double& PermPSwapNd) {
  return new TKroneckerLL(GraphPt, ParamMtx, PermPSwapNd);
}

PKroneckerLL TKroneckerLL::New(const PNGraph& GraphPt, const TKronMtx& ParamMtx, const TIntV& NodeIdPermV, const double& PermPSwapNd) {
  return new TKroneckerLL(GraphPt, ParamMtx, NodeIdPermV, PermPSwapNd);
}

void TKroneckerLL::SetPerm(const char& PermId) {
  if (PermId == 'o') { SetOrderPerm(); }
  else if (PermId == 'd') { SetDegPerm(); }
  else if (PermId == 'r') { SetRndPerm(); }
  else if (PermId == 'b') { SetBestDegPerm(); }
  else FailR("Unknown permutation type (o,d,r)");
}

void TKroneckerLL::SetOrderPerm() {
  NodePerm.Gen(Nodes, 0);
  for (int i = 0; i < Graph->GetNodes(); i++) {
    NodePerm.Add(i); }
  SetIPerm(NodePerm);
}

void TKroneckerLL::SetRndPerm() {
  NodePerm.Gen(Nodes, 0);
  for (int i = 0; i < Graph->GetNodes(); i++) {
    NodePerm.Add(i); }
  NodePerm.Shuffle(TKronMtx::Rnd);
  SetIPerm(NodePerm);
}

void TKroneckerLL::SetDegPerm() {
  TIntPrV DegNIdV;
  for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    DegNIdV.Add(TIntPr(NI.GetDeg(), NI.GetId()));
  }
  DegNIdV.Sort(false);
  NodePerm.Gen(DegNIdV.Len(), 0);
  for (int i = 0; i < DegNIdV.Len(); i++) {
    NodePerm.Add(DegNIdV[i].Val2);
  }
  SetIPerm(NodePerm);
}

/// !!!!! MYUNGHWAN, CHECK!
void TKroneckerLL::SetBestDegPerm() {
  NodePerm.Gen(Nodes);
  const int NZero = ProbMtx.GetDim();
  TFltIntPrV DegV(Nodes), CDegV(Nodes);
  TFltV Row(NZero);
  TFltV Col(NZero);
  for(int i = 0; i < NZero; i++) {
	  for(int j = 0; j < NZero; j++) {
		  Row[i] += ProbMtx.At(i, j);
		  Col[i] += ProbMtx.At(j, i);
	  }
  }

  for(int i = 0; i < Nodes; i++) {
	  TNGraph::TNodeI NodeI = Graph->GetNI(i);
	  int NId = i;
	  double RowP = 1.0, ColP = 1.0;
	  for(int j = 0; j < KronIters; j++) {
		  int Bit = NId % NZero;
		  RowP *= Row[Bit];		ColP *= Col[Bit];
		  NId /= NZero;
	  }
	  CDegV[i] = TFltIntPr(RowP + ColP, i);
	  DegV[i] = TFltIntPr(NodeI.GetDeg(), i);
  }
  DegV.Sort(false);		CDegV.Sort(false);
  for(int i = 0; i < Nodes; i++) {
	  NodePerm[DegV[i].Val2] = CDegV[i].Val2;
  }
  SetIPerm(NodePerm);
}

/// !!!!! MYUNGHWAN, CHECK!
void TKroneckerLL::SetIPerm(const TIntV& Perm) {
	InvertPerm.Gen(Perm.Len());
	for (int i = 0; i < Perm.Len(); i++) {
		InvertPerm[Perm[i]] = i;
	}
}

void TKroneckerLL::SetGraph(const PNGraph& GraphPt) {
  Graph = GraphPt;
  bool NodesOk = true;
  // check that nodes IDs are {0,1,..,Nodes-1}
  for (int nid = 0; nid < Graph->GetNodes(); nid++) {
    if (! Graph->IsNode(nid)) { NodesOk=false; break; } }
  if (! NodesOk) {
    TIntV NIdV;  GraphPt->GetNIdV(NIdV);
    Graph = TSnap::GetSubGraph(GraphPt, NIdV, true);
    for (int nid = 0; nid < Graph->GetNodes(); nid++) {
      IAssert(Graph->IsNode(nid)); }
  }
  Nodes = Graph->GetNodes();
  IAssert(LLMtx.GetDim() > 1 && LLMtx.Len() == ProbMtx.Len());
  KronIters = (int) ceil(log(double(Nodes)) / log(double(ProbMtx.GetDim())));
  // edge vector (for swap-edge permutation proposal)
//  if (PermSwapNodeProb < 1.0) { /// !!!!! MYUNGHWAN, CHECK! WHY IS THIS COMMENTED OUT
    GEdgeV.Gen(Graph->GetEdges(), 0);
    for (TNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
      if (EI.GetSrcNId() != EI.GetDstNId()) {
        GEdgeV.Add(TIntTr(EI.GetSrcNId(), EI.GetDstNId(), -1));
      }
    }
//  }

  RealNodes = Nodes;
  RealEdges = Graph->GetEdges();
  LEdgeV = TIntTrV();
  LSelfEdge = 0;
}


void TKroneckerLL::AppendIsoNodes() {
  Nodes = (int) pow((double)ProbMtx.GetDim(), KronIters);
  // add nodes until filling the Kronecker graph model
  for (int nid = Graph->GetNodes(); nid < Nodes; nid++) {
	  Graph->AddNode(nid);
  }
}

/// !!!!! MYUNGHWAN, CHECK!
void TKroneckerLL::RestoreGraph(const bool RestoreNodes) {
	//	remove from Graph
	int NId1, NId2;
	for (int e = 0; e < LEdgeV.Len(); e++) {
    	NId1 = LEdgeV[e].Val1;  NId2 = LEdgeV[e].Val2;
		Graph->DelEdge(NId1, NId2);
//		GEdgeV.DelIfIn(LEdgeV[e]);
	}
	if(LEdgeV.Len() - LSelfEdge)
		GEdgeV.Del(GEdgeV.Len() - LEdgeV.Len() + LSelfEdge, GEdgeV.Len() - 1);
	LEdgeV.Clr();
	LSelfEdge = 0;

	if(RestoreNodes) {
		for(int i = Graph->GetNodes()-1; i >= RealNodes; i--) {
			Graph->DelNode(i);
		}
	}
}

double TKroneckerLL::GetFullGraphLL() const {
  // the number of times a seed matrix element appears in
  // the full kronecker adjacency matrix after KronIter
  // kronecker multiplications
  double ElemCnt = 1;
  const double dim = LLMtx.GetDim();
  // count number of times x appears in the full kronecker matrix
  for (int i = 1; i < KronIters; i++) {
    ElemCnt = dim*dim*ElemCnt + TMath::Power(dim, 2*i);
  }
  return ElemCnt * LLMtx.GetMtxSum();
}

double TKroneckerLL::GetFullRowLL(int RowId) const {
  double RowLL = 0.0;
  const int MtxDim = LLMtx.GetDim();
  for (int level = 0; level < KronIters; level++) {
    RowLL += LLMtx.GetRowSum(RowId % MtxDim);
    RowId /= MtxDim;
  }
  return RowLL;
}

double TKroneckerLL::GetFullColLL(int ColId) const {
  double ColLL = 0.0;
  const int MtxDim = LLMtx.GetDim();
  for (int level = 0; level < KronIters; level++) {
    ColLL += LLMtx.GetColSum(ColId % MtxDim);
    ColId /= MtxDim;
  }
  return ColLL;
}

double TKroneckerLL::GetEmptyGraphLL() const {
  double LL = 0;
  for (int NId1 = 0; NId1 < LLMtx.GetNodes(KronIters); NId1++) {
    for (int NId2 = 0; NId2 < LLMtx.GetNodes(KronIters); NId2++) {
      LL = LL + LLMtx.GetNoEdgeLL(NId1, NId2, KronIters);
    }
  }
  return LL;
}

// 2nd prder Taylor approximation, log(1-x) ~ -x - 0.5x^2
double TKroneckerLL::GetApxEmptyGraphLL() const {
  double Sum=0.0, SumSq=0.0;
  for (int i = 0; i < ProbMtx.Len(); i++) {
    Sum += ProbMtx.At(i);
    SumSq += TMath::Sqr(ProbMtx.At(i));
  }
  return -pow(Sum, KronIters) - 0.5*pow(SumSq, KronIters);
}

void TKroneckerLL::InitLL(const TFltV& ParamV) {
  InitLL(TKronMtx(ParamV));
}

void TKroneckerLL::InitLL(const TKronMtx& ParamMtx) {
  IAssert(ParamMtx.IsProbMtx());
  ProbMtx = ParamMtx;
  ProbMtx.GetLLMtx(LLMtx);
  LogLike = TKronMtx::NInf;
  if (GradV.Len() != ProbMtx.Len()) {
    GradV.Gen(ProbMtx.Len()); }
  GradV.PutAll(0.0);
}

void TKroneckerLL::InitLL(const PNGraph& GraphPt, const TKronMtx& ParamMtx) {
  IAssert(ParamMtx.IsProbMtx());
  ProbMtx = ParamMtx;
  ProbMtx.GetLLMtx(LLMtx);
  SetGraph(GraphPt);
  LogLike = TKronMtx::NInf;
  if (GradV.Len() != ProbMtx.Len()) {
    GradV.Gen(ProbMtx.Len()); }
  GradV.PutAll(0.0);
}

// exact graph log-likelihood, takes O(N^2 + E)
double TKroneckerLL::CalcGraphLL() {
  LogLike = GetEmptyGraphLL(); // takes O(N^2)
  for (int nid = 0; nid < Nodes; nid++) {
    const TNGraph::TNodeI Node = Graph->GetNI(nid);
    const int SrcNId = NodePerm[nid];
    for (int e = 0; e < Node.GetOutDeg(); e++) {
      const int DstNId = NodePerm[Node.GetOutNId(e)];
      LogLike = LogLike - LLMtx.GetNoEdgeLL(SrcNId, DstNId, KronIters)
        + LLMtx.GetEdgeLL(SrcNId, DstNId, KronIters);
    }
  }
  return LogLike;
}

// approximate graph log-likelihood, takes O(E + N_0)
double TKroneckerLL::CalcApxGraphLL() {
  LogLike = GetApxEmptyGraphLL(); // O(N_0)
  for (int nid = 0; nid < Nodes; nid++) {
    const TNGraph::TNodeI Node = Graph->GetNI(nid);
    const int SrcNId = NodePerm[nid];
    for (int e = 0; e < Node.GetOutDeg(); e++) {
      const int DstNId = NodePerm[Node.GetOutNId(e)];
      LogLike = LogLike - LLMtx.GetApxNoEdgeLL(SrcNId, DstNId, KronIters)
        + LLMtx.GetEdgeLL(SrcNId, DstNId, KronIters);
    }
  }
  return LogLike;
}

// Used in TKroneckerLL::SwapNodesLL: DeltaLL if we
// add the node to the matrix (node gets/creates all
// of its in- and out-edges).
// Zero is for the empty row/column (isolated node)
double TKroneckerLL::NodeLLDelta(const int& NId) const {
  if (! Graph->IsNode(NId)) { return 0.0; } // zero degree node
  double Delta = 0.0;
  const TNGraph::TNodeI Node = Graph->GetNI(NId);
  // out-edges
  const int SrcRow = NodePerm[NId];
  for (int e = 0; e < Node.GetOutDeg(); e++) {
    const int DstCol = NodePerm[Node.GetOutNId(e)];
    Delta += - LLMtx.GetApxNoEdgeLL(SrcRow, DstCol, KronIters)
      + LLMtx.GetEdgeLL(SrcRow, DstCol, KronIters);
  }
  //in-edges
  const int SrcCol = NodePerm[NId];
  for (int e = 0; e < Node.GetInDeg(); e++) {
    const int DstRow = NodePerm[Node.GetInNId(e)];
    Delta += - LLMtx.GetApxNoEdgeLL(DstRow, SrcCol, KronIters)
      + LLMtx.GetEdgeLL(DstRow, SrcCol, KronIters);
  }
  // double counted self-edge
  if (Graph->IsEdge(NId, NId)) {
    Delta += + LLMtx.GetApxNoEdgeLL(SrcRow, SrcCol, KronIters)
      - LLMtx.GetEdgeLL(SrcRow, SrcCol, KronIters);
    IAssert(SrcRow == SrcCol);
  }
  return Delta;
}

// swapping two nodes, only need to go over two rows and columns
double TKroneckerLL::SwapNodesLL(const int& NId1, const int& NId2) {
  // subtract old LL (remove nodes)
  LogLike = LogLike - NodeLLDelta(NId1) - NodeLLDelta(NId2);
  const int PrevId1 = NodePerm[NId1], PrevId2 = NodePerm[NId2];
  // double-counted edges
  if (Graph->IsEdge(NId1, NId2)) {
    LogLike += - LLMtx.GetApxNoEdgeLL(PrevId1, PrevId2, KronIters)
      + LLMtx.GetEdgeLL(PrevId1, PrevId2, KronIters); }
  if (Graph->IsEdge(NId2, NId1)) {
    LogLike += - LLMtx.GetApxNoEdgeLL(PrevId2, PrevId1, KronIters)
      + LLMtx.GetEdgeLL(PrevId2, PrevId1, KronIters); }
  // swap
  NodePerm.Swap(NId1, NId2);
  InvertPerm.Swap(NodePerm[NId1], NodePerm[NId2]);
  // add new LL (add nodes)
  LogLike = LogLike + NodeLLDelta(NId1) + NodeLLDelta(NId2);
  const int NewId1 = NodePerm[NId1], NewId2 = NodePerm[NId2];
  // correct for double-counted edges
  if (Graph->IsEdge(NId1, NId2)) {
    LogLike += + LLMtx.GetApxNoEdgeLL(NewId1, NewId2, KronIters)
      - LLMtx.GetEdgeLL(NewId1, NewId2, KronIters); }
  if (Graph->IsEdge(NId2, NId1)) {
    LogLike += + LLMtx.GetApxNoEdgeLL(NewId2, NewId1, KronIters)
      - LLMtx.GetEdgeLL(NewId2, NewId1, KronIters); }
  return LogLike;
}

// metropolis sampling from P(permutation|graph)
bool TKroneckerLL::SampleNextPerm(int& NId1, int& NId2) {
  // pick 2 uniform nodes and swap
  if (TKronMtx::Rnd.GetUniDev() < PermSwapNodeProb) {
    NId1 = TKronMtx::Rnd.GetUniDevInt(Nodes);
    NId2 = TKronMtx::Rnd.GetUniDevInt(Nodes);
    while (NId2 == NId1) { NId2 = TKronMtx::Rnd.GetUniDevInt(Nodes); }
  } else {
    // pick uniform edge and swap endpoints (slow as it moves around high degree nodes)
    const int e = TKronMtx::Rnd.GetUniDevInt(GEdgeV.Len());
    NId1 = GEdgeV[e].Val1;  NId2 = GEdgeV[e].Val2;
  }
  const double U = TKronMtx::Rnd.GetUniDev();
  const double OldLL = LogLike;
  const double NewLL = SwapNodesLL(NId1, NId2);
  const double LogU = log(U);
  if (LogU > NewLL - OldLL) { // reject
    LogLike = OldLL;
    NodePerm.Swap(NId2, NId1); //swap back
	InvertPerm.Swap(NodePerm[NId2], NodePerm[NId1]); // swap back
    return false;
  }
  return true; // accept new sample
}

// exact gradient of an empty graph, O(N^2)
double TKroneckerLL::GetEmptyGraphDLL(const int& ParamId) const {
  double DLL = 0.0;
  for (int NId1 = 0; NId1 < Nodes; NId1++) {
    for (int NId2 = 0; NId2 < Nodes; NId2++) {
      DLL += LLMtx.GetNoEdgeDLL(ParamId, NodePerm[NId1], NodePerm[NId2], KronIters);
    }
  }
  return DLL;
}

// approx gradient, using 2nd order Taylor approximation, O(N_0^2)
double TKroneckerLL::GetApxEmptyGraphDLL(const int& ParamId) const {
  double Sum=0.0, SumSq=0.0;
  for (int i = 0; i < ProbMtx.Len(); i++) {
    Sum += ProbMtx.At(i);
    SumSq += TMath::Sqr(ProbMtx.At(i));
  }
  // d/dx -sum(x_i) - 0.5sum(x_i^2) = d/dx sum(theta)^k - 0.5 sum(theta^2)^k
  return -KronIters*pow(Sum, KronIters-1) - KronIters*pow(SumSq, KronIters-1)*ProbMtx.At(ParamId);
}

// exact graph gradient, runs O(N^2)
const TFltV& TKroneckerLL::CalcGraphDLL() {
  for (int ParamId = 0; ParamId < LLMtx.Len(); ParamId++) {
    double DLL = 0.0;
    for (int NId1 = 0; NId1 < Nodes; NId1++) {
      for (int NId2 = 0; NId2 < Nodes; NId2++) {
        if (Graph->IsEdge(NId1, NId2)) {
          DLL += LLMtx.GetEdgeDLL(ParamId, NodePerm[NId1], NodePerm[NId2], KronIters);
        } else {
          DLL += LLMtx.GetNoEdgeDLL(ParamId, NodePerm[NId1], NodePerm[NId2], KronIters);
        }
      }
    }
    GradV[ParamId] = DLL;
  }
  return GradV;
}

// slow (but correct) approximate gradient, runs O(N^2)
const TFltV& TKroneckerLL::CalcFullApxGraphDLL() {
  for (int ParamId = 0; ParamId < LLMtx.Len(); ParamId++) {
    double DLL = 0.0;
    for (int NId1 = 0; NId1 < Nodes; NId1++) {
      for (int NId2 = 0; NId2 < Nodes; NId2++) {
        if (Graph->IsEdge(NId1, NId2)) {
          DLL += LLMtx.GetEdgeDLL(ParamId, NodePerm[NId1], NodePerm[NId2], KronIters);
        } else {
          DLL += LLMtx.GetApxNoEdgeDLL(ParamId, NodePerm[NId1], NodePerm[NId2], KronIters);
        }
      }
    }
    GradV[ParamId] = DLL;
  }
  return GradV;
}

// fast approximate gradient, runs O(E)
const TFltV& TKroneckerLL::CalcApxGraphDLL() {
  for (int ParamId = 0; ParamId < LLMtx.Len(); ParamId++) {
    double DLL = GetApxEmptyGraphDLL(ParamId);
    for (int nid = 0; nid < Nodes; nid++) {
      const TNGraph::TNodeI Node = Graph->GetNI(nid);
      const int SrcNId = NodePerm[nid];
      for (int e = 0; e < Node.GetOutDeg(); e++) {
        const int DstNId = NodePerm[Node.GetOutNId(e)];
        DLL = DLL - LLMtx.GetApxNoEdgeDLL(ParamId, SrcNId, DstNId, KronIters)
          + LLMtx.GetEdgeDLL(ParamId, SrcNId, DstNId, KronIters);
      }
    }
    GradV[ParamId] = DLL;
  }
  return GradV;
}

// Used in TKroneckerLL::UpdateGraphDLL: DeltaDLL if we
// add the node to the empty matrix/graph (node
// gets/creates all of its in- and out-edges).
double TKroneckerLL::NodeDLLDelta(const int ParamId, const int& NId) const {
  if (! Graph->IsNode(NId)) { return 0.0; } // zero degree node
  double Delta = 0.0;
  const TNGraph::TNodeI Node = Graph->GetNI(NId);
  const int SrcRow = NodePerm[NId];
  for (int e = 0; e < Node.GetOutDeg(); e++) {
    const int DstCol = NodePerm[Node.GetOutNId(e)];
    Delta += - LLMtx.GetApxNoEdgeDLL(ParamId, SrcRow, DstCol, KronIters)
      + LLMtx.GetEdgeDLL(ParamId, SrcRow, DstCol, KronIters);
  }
  const int SrcCol = NodePerm[NId];
  for (int e = 0; e < Node.GetInDeg(); e++) {
    const int DstRow = NodePerm[Node.GetInNId(e)];
    Delta += - LLMtx.GetApxNoEdgeDLL(ParamId, DstRow, SrcCol, KronIters)
      + LLMtx.GetEdgeDLL(ParamId, DstRow, SrcCol, KronIters);
  }
  // double counter self-edge
  if (Graph->IsEdge(NId, NId)) {
    Delta += + LLMtx.GetApxNoEdgeDLL(ParamId, SrcRow, SrcCol, KronIters)
      - LLMtx.GetEdgeDLL(ParamId, SrcRow, SrcCol, KronIters);
    IAssert(SrcRow == SrcCol);
  }
  return Delta;
}

// given old DLL and new permutation, efficiently updates the DLL
// permutation is new, but DLL is old
void TKroneckerLL::UpdateGraphDLL(const int& SwapNId1, const int& SwapNId2) {
  for (int ParamId = 0; ParamId < LLMtx.Len(); ParamId++) {
    // permutation before the swap (swap back to previous position)
    NodePerm.Swap(SwapNId1, SwapNId2);
    // subtract old DLL
    TFlt& DLL = GradV[ParamId];
    DLL = DLL - NodeDLLDelta(ParamId, SwapNId1) - NodeDLLDelta(ParamId, SwapNId2);
    // double-counted edges
    const int PrevId1 = NodePerm[SwapNId1], PrevId2 = NodePerm[SwapNId2];
    if (Graph->IsEdge(SwapNId1, SwapNId2)) {
      DLL += - LLMtx.GetApxNoEdgeDLL(ParamId, PrevId1, PrevId2, KronIters)
        + LLMtx.GetEdgeDLL(ParamId, PrevId1, PrevId2, KronIters); }
    if (Graph->IsEdge(SwapNId2, SwapNId1)) {
      DLL += - LLMtx.GetApxNoEdgeDLL(ParamId, PrevId2, PrevId1, KronIters)
        + LLMtx.GetEdgeDLL(ParamId, PrevId2, PrevId1, KronIters); }
    // permutation after the swap (restore the swap)
    NodePerm.Swap(SwapNId1, SwapNId2);
    // add new DLL
    DLL = DLL + NodeDLLDelta(ParamId, SwapNId1) + NodeDLLDelta(ParamId, SwapNId2);
    const int NewId1 = NodePerm[SwapNId1], NewId2 = NodePerm[SwapNId2];
    // double-counted edges
    if (Graph->IsEdge(SwapNId1, SwapNId2)) {
      DLL += + LLMtx.GetApxNoEdgeDLL(ParamId, NewId1, NewId2, KronIters)
        - LLMtx.GetEdgeDLL(ParamId, NewId1, NewId2, KronIters); }
    if (Graph->IsEdge(SwapNId2, SwapNId1)) {
      DLL += + LLMtx.GetApxNoEdgeDLL(ParamId, NewId2, NewId1, KronIters)
        - LLMtx.GetEdgeDLL(ParamId, NewId2, NewId1, KronIters); }
  }
}

void TKroneckerLL::SampleGradient(const int& WarmUp, const int& NSamples, double& AvgLL, TFltV& AvgGradV) {
  printf("SampleGradient: %s (%s warm-up):", TInt::GetMegaStr(NSamples).CStr(), TInt::GetMegaStr(WarmUp).CStr());
  int NId1=0, NId2=0, NAccept=0;
  TExeTm ExeTm1;
  if (WarmUp > 0) {
    CalcApxGraphLL();
    for (int s = 0; s < WarmUp; s++) { SampleNextPerm(NId1, NId2); }
    printf("  warm-up:%s,", ExeTm1.GetTmStr());  ExeTm1.Tick();
  }
  CalcApxGraphLL(); // re-calculate LL (due to numerical errors)
  CalcApxGraphDLL();
  AvgLL = 0;
  AvgGradV.Gen(LLMtx.Len());  AvgGradV.PutAll(0.0);
  printf("  sampl");
  for (int s = 0; s < NSamples; s++) {
    if (SampleNextPerm(NId1, NId2)) { // new permutation
      UpdateGraphDLL(NId1, NId2);  NAccept++; }
    for (int m = 0; m < LLMtx.Len(); m++) { AvgGradV[m] += GradV[m]; }
    AvgLL += GetLL();
  }
  printf("ing");
  AvgLL = AvgLL / double(NSamples);
  for (int m = 0; m < LLMtx.Len(); m++) {
    AvgGradV[m] = AvgGradV[m] / double(NSamples); }
  printf(":%s (%.0f/s), accept %.1f%%\n", ExeTm1.GetTmStr(), double(NSamples)/ExeTm1.GetSecs(),
    double(100*NAccept)/double(NSamples));
}

double TKroneckerLL::GradDescent(const int& NIter, const double& LrnRate, double MnStep, double MxStep, const int& WarmUp, const int& NSamples) {
  printf("\n----------------------------------------------------------------------\n");
  printf("Fitting graph on %d nodes, %d edges\n", Graph->GetNodes(), Graph->GetEdges());
  printf("Kron iters:  %d (== %d nodes)\n\n", KronIters(), ProbMtx.GetNodes(KronIters()));
  TExeTm IterTm, TotalTm;
  double OldLL=-1e10, CurLL=0;
  const double EZero = pow((double) Graph->GetEdges(), 1.0/double(KronIters));
  TFltV CurGradV, LearnRateV(GetParams()), LastStep(GetParams());
  LearnRateV.PutAll(LrnRate);
  TKronMtx NewProbMtx = ProbMtx;

  if(DebugMode) {  /// !!!!! MYUNGHWAN, CHECK!
	  LLV.Gen(NIter, 0);
	  MtxV.Gen(NIter, 0);
  }

  for (int Iter = 0; Iter < NIter; Iter++) {
    printf("%03d] ", Iter);
    SampleGradient(WarmUp, NSamples, CurLL, CurGradV);
    for (int p = 0; p < GetParams(); p++) {
      LearnRateV[p] *= 0.95;
      if (Iter < 1) {
        while (fabs(LearnRateV[p]*CurGradV[p]) > MxStep) { LearnRateV[p] *= 0.95; }
        while (fabs(LearnRateV[p]*CurGradV[p]) < 0.02) { LearnRateV[p] *= (1.0/0.95); } // move more
      } else {
        // set learn rate so that move for each parameter is inside the [MnStep, MxStep]
        while (fabs(LearnRateV[p]*CurGradV[p]) > MxStep) { LearnRateV[p] *= 0.95; printf(".");}
        while (fabs(LearnRateV[p]*CurGradV[p]) < MnStep) { LearnRateV[p] *= (1.0/0.95); printf("*");}
        if (MxStep > 3*MnStep) { MxStep *= 0.95; }
      }
      NewProbMtx.At(p) = ProbMtx.At(p) + LearnRateV[p]*CurGradV[p];
      if (NewProbMtx.At(p) > 0.9999) { NewProbMtx.At(p)=0.9999; }
      if (NewProbMtx.At(p) < 0.0001) { NewProbMtx.At(p)=0.0001; }
    }
    printf("  trueE0: %.2f (%d),  estE0: %.2f (%d),  ERR: %f\n", EZero, Graph->GetEdges(),
      ProbMtx.GetMtxSum(), ProbMtx.GetEdges(KronIters), fabs(EZero-ProbMtx.GetMtxSum()));
    printf("  currLL: %.4f, deltaLL: %.4f\n", CurLL, CurLL-OldLL); // positive is good
    for (int p = 0; p < GetParams(); p++) {
      printf("    %d]  %f  <--  %f + %9f   Grad: %9.1f   Rate: %g\n", p, NewProbMtx.At(p),
        ProbMtx.At(p), (double)(LearnRateV[p]*CurGradV[p]), CurGradV[p](), LearnRateV[p]());
    }
    if (Iter+1 < NIter) { // skip last update
      ProbMtx = NewProbMtx;  ProbMtx.GetLLMtx(LLMtx); }
    OldLL=CurLL;
    printf("\n");  fflush(stdout);

	if(DebugMode) {  /// !!!!! MYUNGHWAN, CHECK!
		LLV.Add(CurLL);
		MtxV.Add(NewProbMtx);
	}
  }
  printf("TotalExeTm: %s %g\n", TotalTm.GetStr(), TotalTm.GetSecs());
  ProbMtx.Dump("FITTED PARAMS", false);
  return CurLL;
}

double TKroneckerLL::GradDescent2(const int& NIter, const double& LrnRate, double MnStep, double MxStep, const int& WarmUp, const int& NSamples) {
  printf("\n----------------------------------------------------------------------\n");
  printf("GradDescent2\n");
  printf("Fitting graph on %d nodes, %d edges\n", Graph->GetNodes(), Graph->GetEdges());
  printf("Skip moves that make likelihood smaller\n");
  printf("Kron iters:  %d (== %d nodes)\n\n", KronIters(), ProbMtx.GetNodes(KronIters()));
  TExeTm IterTm, TotalTm;
  double CurLL=0, NewLL=0;
  const double EZero = pow((double) Graph->GetEdges(), 1.0/double(KronIters));
  TFltV CurGradV, NewGradV, LearnRateV(GetParams()), LastStep(GetParams());
  LearnRateV.PutAll(LrnRate);
  TKronMtx NewProbMtx=ProbMtx, CurProbMtx=ProbMtx;
  bool GoodMove = false;
  // Start
  for (int Iter = 0; Iter < NIter; Iter++) {
    printf("%03d] ", Iter);
    if (! GoodMove) { SampleGradient(WarmUp, NSamples, CurLL, CurGradV); }
    CurProbMtx = ProbMtx;
    // update parameters
    for (int p = 0; p < GetParams(); p++) {
      while (fabs(LearnRateV[p]*CurGradV[p]) > MxStep) { LearnRateV[p] *= 0.95; printf(".");}
      while (fabs(LearnRateV[p]*CurGradV[p]) < MnStep) { LearnRateV[p] *= (1.0/0.95); printf("*");}
      NewProbMtx.At(p) = CurProbMtx.At(p) + LearnRateV[p]*CurGradV[p];
      if (NewProbMtx.At(p) > 0.9999) { NewProbMtx.At(p)=0.9999; }
      if (NewProbMtx.At(p) < 0.0001) { NewProbMtx.At(p)=0.0001; }
      LearnRateV[p] *= 0.95;
    }
    printf("  ");
    ProbMtx=NewProbMtx;  ProbMtx.GetLLMtx(LLMtx);
    SampleGradient(WarmUp, NSamples, NewLL, NewGradV);
    if (NewLL > CurLL) { // accept the move
      printf("== Good move:\n");
      printf("  trueE0: %.2f (%d),  estE0: %.2f (%d),  ERR: %f\n", EZero, Graph->GetEdges(),
        ProbMtx.GetMtxSum(), ProbMtx.GetEdges(KronIters), fabs(EZero-ProbMtx.GetMtxSum()));
      printf("  currLL: %.4f  deltaLL: %.4f\n", CurLL, NewLL-CurLL); // positive is good
      for (int p = 0; p < GetParams(); p++) {
        printf("    %d]  %f  <--  %f + %9f   Grad: %9.1f   Rate: %g\n", p, NewProbMtx.At(p),
          CurProbMtx.At(p), (double)(LearnRateV[p]*CurGradV[p]), CurGradV[p](), LearnRateV[p]()); }
      CurLL = NewLL;
      CurGradV = NewGradV;
      GoodMove = true;
    } else {
      printf("** BAD move:\n");
      printf("  *trueE0: %.2f (%d),  estE0: %.2f (%d),  ERR: %f\n", EZero, Graph->GetEdges(),
        ProbMtx.GetMtxSum(), ProbMtx.GetEdges(KronIters), fabs(EZero-ProbMtx.GetMtxSum()));
      printf("  *curLL:  %.4f  deltaLL: %.4f\n", CurLL, NewLL-CurLL); // positive is good
      for (int p = 0; p < GetParams(); p++) {
        printf("   b%d]  %f  <--  %f + %9f   Grad: %9.1f   Rate: %g\n", p, NewProbMtx.At(p),
          CurProbMtx.At(p), (double)(LearnRateV[p]*CurGradV[p]), CurGradV[p](), LearnRateV[p]()); }
      // move to old position
      ProbMtx = CurProbMtx;  ProbMtx.GetLLMtx(LLMtx);
      GoodMove = false;
    }
    printf("\n");  fflush(stdout);
  }
  printf("TotalExeTm: %s %g\n", TotalTm.GetStr(), TotalTm.GetSecs());
  ProbMtx.Dump("FITTED PARAMS\n", false);
  return CurLL;
}

/// !!!!! MYUNGHWAN, CHECK!
// filling in random edges for KronEM
void TKroneckerLL::SetRandomEdges(const int& NEdges, const bool isDir) {
	int count = 0, added = 0, collision = 0;
	const int MtxDim = ProbMtx.GetDim();
	const double MtxSum = ProbMtx.GetMtxSum();
	TVec<TFltIntIntTr> ProbToRCPosV; // row, col position
	double CumProb = 0.0;

	for(int r = 0; r < MtxDim; r++) {
		for(int c = 0; c < MtxDim; c++) {
			const double Prob = ProbMtx.At(r, c);
			if (Prob > 0.0) {
				CumProb += Prob;
				ProbToRCPosV.Add(TFltIntIntTr(CumProb/MtxSum, r, c));
			}
		}
	}

	int Rng, Row, Col, n, NId1, NId2;
	while(added < NEdges) {
		Rng = Nodes;	Row = 0;	Col = 0;
		for (int iter = 0; iter < KronIters; iter++) {
			const double& Prob = TKronMtx::Rnd.GetUniDev();
			n = 0; while(Prob > ProbToRCPosV[n].Val1) { n++; }
			const int MtxRow = ProbToRCPosV[n].Val2;
			const int MtxCol = ProbToRCPosV[n].Val3;
			Rng /= MtxDim;
			Row += MtxRow * Rng;
			Col += MtxCol * Rng;
		}

		count++;

		NId1 = InvertPerm[Row];	NId2 = InvertPerm[Col];

		//	Check conflicts
		if(EMType != kronEdgeMiss && IsObsEdge(NId1, NId2)) {
			continue;
		}

		if (! Graph->IsEdge(NId1, NId2)) {
			Graph->AddEdge(NId1, NId2);
			if(NId1 != NId2)	{ GEdgeV.Add(TIntTr(NId1, NId2, LEdgeV.Len())); }
			else { LSelfEdge++; }
			LEdgeV.Add(TIntTr(NId1, NId2, GEdgeV.Len()-1));
			added++;
			if (! isDir) {
				if (NId1 != NId2) {
				   Graph->AddEdge(NId2, NId1);
				   GEdgeV.Add(TIntTr(NId2, NId1, LEdgeV.Len()));
				   LEdgeV.Add(TIntTr(NId2, NId1, GEdgeV.Len()-1));
				   added++;
				}
			}
		} else { collision ++; }
	}
//	printf("total = %d / added = %d / collision = %d\n", count, added, collision);
}

// sampling setup for KronEM
void TKroneckerLL::MetroGibbsSampleSetup(const int& WarmUp) {
	double alpha = log(ProbMtx.GetMtxSum()) / log(double(ProbMtx.GetDim()));
	int NId1 = 0, NId2 = 0;
	int NMissing;

	RestoreGraph(false);
	if(EMType == kronEdgeMiss) {
		CalcApxGraphLL();
		for (int s = 0; s < WarmUp; s++)	SampleNextPerm(NId1, NId2);
	}

	if(EMType == kronFutureLink) {
		NMissing = (int) (pow(ProbMtx.GetMtxSum(), KronIters) - pow(double(RealNodes), alpha));
	} else if(EMType == kronEdgeMiss) {
		NMissing = MissEdges;
	} else {
		NMissing = (int) (pow(ProbMtx.GetMtxSum(), KronIters) * (1.0 - pow(double(RealNodes) / double(Nodes), 2)));
	}
	NMissing = (NMissing < 1) ? 1 : NMissing;

	SetRandomEdges(NMissing, true);

	CalcApxGraphLL();
	for (int s = 0; s < WarmUp; s++)	SampleNextPerm(NId1, NId2);
}

// Metropolis-Hastings steps for KronEM
void TKroneckerLL::MetroGibbsSampleNext(const int& WarmUp, const bool DLLUpdate) {
	int NId1 = 0, NId2 = 0, hit = 0, GId = 0;
	TIntTr EdgeToRemove, NewEdge;
	double RndAccept;

	if(LEdgeV.Len()) {
		for(int i = 0; i < WarmUp; i++) {
			hit = TKronMtx::Rnd.GetUniDevInt(LEdgeV.Len());

			NId1 = LEdgeV[hit].Val1;	NId2 = LEdgeV[hit].Val2;
			GId = LEdgeV[hit].Val3;
			SetRandomEdges(1, true);
			NewEdge = LEdgeV.Last();

			RndAccept = (1.0 - exp(LLMtx.GetEdgeLL(NewEdge.Val1, NewEdge.Val2, KronIters))) / (1.0 - exp(LLMtx.GetEdgeLL(NId1, NId2, KronIters)));
			RndAccept = (RndAccept > 1.0) ? 1.0 : RndAccept;

			if(TKronMtx::Rnd.GetUniDev() > RndAccept) { //	reject
				Graph->DelEdge(NewEdge.Val1, NewEdge.Val2);
				if(NewEdge.Val1 != NewEdge.Val2) {	GEdgeV.DelLast();	}
				else {	LSelfEdge--;	}
				LEdgeV.DelLast();
			} else {	//	accept
				Graph->DelEdge(NId1, NId2);
				LEdgeV[hit] = LEdgeV.Last();
				LEdgeV.DelLast();
				if(NId1 == NId2) {
					LSelfEdge--;
					if(NewEdge.Val1 != NewEdge.Val2) {
						GEdgeV[GEdgeV.Len()-1].Val3 = hit;
					}
				} else {
					IAssertR(GEdgeV.Last().Val3 >= 0, "Invalid indexing");

					GEdgeV[GId] = GEdgeV.Last();
					if(NewEdge.Val1 != NewEdge.Val2) {
						GEdgeV[GId].Val3 = hit;
					}
					LEdgeV[GEdgeV[GId].Val3].Val3 = GId;
					GEdgeV.DelLast();
				}

      			LogLike += LLMtx.GetApxNoEdgeLL(EdgeToRemove.Val1, EdgeToRemove.Val2, KronIters) - LLMtx.GetEdgeLL(EdgeToRemove.Val1, EdgeToRemove.Val2, KronIters);
      			LogLike += -LLMtx.GetApxNoEdgeLL(NewEdge.Val1, NewEdge.Val2, KronIters) + LLMtx.GetEdgeLL(NewEdge.Val1, NewEdge.Val2, KronIters);

				if(DLLUpdate) {
  					for (int p = 0; p < LLMtx.Len(); p++) {
						GradV[p] += LLMtx.GetApxNoEdgeDLL(p, EdgeToRemove.Val1, EdgeToRemove.Val2, KronIters) - LLMtx.GetEdgeDLL(p, EdgeToRemove.Val1, EdgeToRemove.Val2, KronIters);
						GradV[p] += -LLMtx.GetApxNoEdgeDLL(p, NewEdge.Val1, NewEdge.Val2, KronIters) + LLMtx.GetEdgeDLL(p, NewEdge.Val1, NewEdge.Val2, KronIters);
					}
				}
			}
		}
	}

//	CalcApxGraphLL();
	for (int s = 0; s < WarmUp; s++) {
		if(SampleNextPerm(NId1, NId2)) {
			if(DLLUpdate)	UpdateGraphDLL(NId1, NId2);
		}
	}
}

// E-step in KronEM
void TKroneckerLL::RunEStep(const int& GibbsWarmUp, const int& WarmUp, const int& NSamples, TFltV& LLV, TVec<TFltV>& DLLV) {
	TExeTm ExeTm, TotalTm;
	LLV.Gen(NSamples, 0);
	DLLV.Gen(NSamples, 0);

	ExeTm.Tick();
	for(int i = 0; i < 2; i++)	MetroGibbsSampleSetup(WarmUp);
	printf("  Warm-Up [%u] : %s\n", WarmUp, ExeTm.GetTmStr());
	CalcApxGraphLL();
	for(int i = 0; i < GibbsWarmUp; i++)	MetroGibbsSampleNext(10, false);
	printf("  Gibbs Warm-Up [%u] : %s\n", GibbsWarmUp, ExeTm.GetTmStr());

	ExeTm.Tick();
	CalcApxGraphLL();
	CalcApxGraphDLL();
	for(int i = 0; i < NSamples; i++) {
		MetroGibbsSampleNext(50, false);

		LLV.Add(LogLike);
		DLLV.Add(GradV);

		int OnePercent = (i+1) % (NSamples / 10);
		if(OnePercent == 0) {
			int TenPercent = ((i+1) / (NSamples / 10)) * 10;
			printf("  %3u%% done : %s\n", TenPercent, ExeTm.GetTmStr());
		}
	}
}

// M-step in KronEM
double TKroneckerLL::RunMStep(const TFltV& LLV, const TVec<TFltV>& DLLV, const int& GradIter, const double& LrnRate, double MnStep, double MxStep) {
	TExeTm IterTm, TotalTm;
	double OldLL=LogLike, CurLL=0;
	const double alpha = log(double(RealEdges)) / log(double(RealNodes));
	const double EZero = pow(double(ProbMtx.GetDim()), alpha);

	TFltV CurGradV(GetParams()), LearnRateV(GetParams()), LastStep(GetParams());
	LearnRateV.PutAll(LrnRate);

	TKronMtx NewProbMtx = ProbMtx;
	const int NSamples = LLV.Len();
	const int ReCalcLen = NSamples / 10;

	for (int s = 0; s < LLV.Len(); s++) {
		CurLL += LLV[s];
		for(int p = 0; p < GetParams(); p++) { CurGradV[p] += DLLV[s][p]; }
	}
	CurLL /= NSamples;
	for(int p = 0; p < GetParams(); p++) { CurGradV[p] /= NSamples; }

	double MaxLL = CurLL;
	TKronMtx MaxProbMtx = ProbMtx;
	TKronMtx OldProbMtx = ProbMtx;

	for (int Iter = 0; Iter < GradIter; Iter++) {
		printf("    %03d] ", Iter+1);
		IterTm.Tick();

		for (int p = 0; p < GetParams(); p++) {
			if (Iter < 1) {
				while (fabs(LearnRateV[p]*CurGradV[p]) > MxStep) { LearnRateV[p] *= 0.95; }
				while (fabs(LearnRateV[p]*CurGradV[p]) < 5 * MnStep) { LearnRateV[p] *= (1.0/0.95); } // move more
			} else {
			// set learn rate so that move for each parameter is inside the [MnStep, MxStep]
				while (fabs(LearnRateV[p]*CurGradV[p]) > MxStep) { LearnRateV[p] *= 0.95; printf(".");}
				while (fabs(LearnRateV[p]*CurGradV[p]) < MnStep) { LearnRateV[p] *= (1.0/0.95); printf("*");}
				if (MxStep > 3*MnStep) { MxStep *= 0.95; }
			}
			NewProbMtx.At(p) = ProbMtx.At(p) + LearnRateV[p]*CurGradV[p];
			if (NewProbMtx.At(p) > 0.9999) { NewProbMtx.At(p)=0.9999; }
			if (NewProbMtx.At(p) < 0.0001) { NewProbMtx.At(p)=0.0001; }
			LearnRateV[p] *= 0.95;
		}
		printf("  trueE0: %.2f (%u from %u),  estE0: %.2f (%u from %u),  ERR: %f\n", EZero, RealEdges(), RealNodes(), ProbMtx.GetMtxSum(), Graph->GetEdges(), Graph->GetNodes(), fabs(EZero-ProbMtx.GetMtxSum()));
		printf("      currLL: %.4f, deltaLL: %.4f\n", CurLL, CurLL-OldLL); // positive is good
		for (int p = 0; p < GetParams(); p++) {
			printf("      %d]  %f  <--  %f + %9f   Grad: %9.1f   Rate: %g\n", p, NewProbMtx.At(p),
			ProbMtx.At(p), (double)(LearnRateV[p]*CurGradV[p]), CurGradV[p](), LearnRateV[p]());
		}

		OldLL=CurLL;
		if(Iter == GradIter - 1) {
			break;
		}

		CurLL = 0;
		CurGradV.PutAll(0.0);
		TFltV OneDLL;

		CalcApxGraphLL();
		CalcApxGraphDLL();

		for(int s = 0; s < NSamples; s++) {
			ProbMtx = OldProbMtx;  ProbMtx.GetLLMtx(LLMtx);
			MetroGibbsSampleNext(10, true);
			ProbMtx = NewProbMtx;  ProbMtx.GetLLMtx(LLMtx);
			if(s % ReCalcLen == ReCalcLen/2) {
				CurLL += CalcApxGraphLL();
				OneDLL = CalcApxGraphDLL();
			} else {
				CurLL += LogLike;
				OneDLL = GradV;
			}
			for(int p = 0; p < GetParams(); p++) {
				CurGradV[p] += OneDLL[p];
			}
		}
		CurLL /= NSamples;

		if(MaxLL < CurLL) {
			MaxLL = CurLL;	MaxProbMtx = ProbMtx;
		}

		printf("    Time: %s\n", IterTm.GetTmStr());
		printf("\n");  fflush(stdout);
	}
	ProbMtx = MaxProbMtx;	ProbMtx.GetLLMtx(LLMtx);

	printf("    FinalLL : %f,   TotalExeTm: %s\n", MaxLL, TotalTm.GetTmStr());
	ProbMtx.Dump("    FITTED PARAMS", false);

	return MaxLL;
}

// KronEM
void TKroneckerLL::RunKronEM(const int& EMIter, const int& GradIter, double LrnRate, double MnStep, double MxStep, const int& GibbsWarmUp, const int& WarmUp, const int& NSamples, const TKronEMType& Type, const int& NMissing) {
	printf("\n----------------------------------------------------------------------\n");
	printf("Fitting graph on %d nodes, %d edges\n", int(RealNodes), int(RealEdges));
	printf("Kron iters:  %d (== %d nodes)\n\n", KronIters(), ProbMtx.GetNodes(KronIters()));

	TFltV LLV(NSamples);
	TVec<TFltV> DLLV(NSamples);
	//int count = 0;

	EMType = Type;
	MissEdges = NMissing;
	AppendIsoNodes();
	SetRndPerm();

	if(DebugMode) {
		LLV.Gen(EMIter, 0);
		MtxV.Gen(EMIter, 0);
	}

	for(int i = 0; i < EMIter; i++) {
		printf("\n----------------------------------------------------------------------\n");
		printf("%03d EM-iter] E-Step\n", i+1);
		RunEStep(GibbsWarmUp, WarmUp, NSamples, LLV, DLLV);
		printf("\n\n");

		printf("%03d EM-iter] M-Step\n", i+1);
		double CurLL = RunMStep(LLV, DLLV, GradIter, LrnRate, MnStep, MxStep);
		printf("\n\n");

		if(DebugMode) {
			LLV.Add(CurLL);
			MtxV.Add(ProbMtx);
		}
	}

	RestoreGraph();
}



void GetMinMax(const TFltPrV& XYValV, double& Min, double& Max, const bool& ResetMinMax) {
  if (ResetMinMax) { Min = TFlt::Mx;  Max = TFlt::Mn; }
  for (int i = 0; i < XYValV.Len(); i++) {
    Min = TMath::Mn(Min, XYValV[i].Val2.Val);
    Max = TMath::Mx(Max, XYValV[i].Val2.Val);
  }
}

void PlotGrad(const TFltPrV& EstLLV, const TFltPrV& TrueLLV, const TVec<TFltPrV>& GradVV, const TFltPrV& AcceptV, const TStr& OutFNm, const TStr& Desc) {
  double Min, Max, Min1, Max1;
  // plot log-likelihood
  { TGnuPlot GP("sLL-"+OutFNm, TStr::Fmt("Log-likelihood (avg 1k samples). %s", Desc.CStr()), true);
  GP.AddPlot(EstLLV, gpwLines, "Esimated LL", "linewidth 1");
  if (! TrueLLV.Empty()) { GP.AddPlot(TrueLLV, gpwLines, "TRUE LL", "linewidth 1"); }
  //GetMinMax(EstLLV, Min, Max, true);  GetMinMax(TrueLLV, Min, Max, false);
  //GP.SetYRange((int)floor(Min-1), (int)ceil(Max+1));
  GP.SetXYLabel("Sample Index (time)", "Log-likelihood");
  GP.SavePng(); }
  // plot accept
  { TGnuPlot GP("sAcc-"+OutFNm, TStr::Fmt("Pct. accepted rnd moves (over 1k samples). %s", Desc.CStr()), true);
  GP.AddPlot(AcceptV, gpwLines, "Pct accepted swaps", "linewidth 1");
  GP.SetXYLabel("Sample Index (time)", "Pct accept permutation swaps");
  GP.SavePng(); }
  // plot grads
  TGnuPlot GPAll("sGradAll-"+OutFNm, TStr::Fmt("Gradient (avg 1k samples). %s", Desc.CStr()), true);
  GetMinMax(GradVV[0], Min1, Max1, true);
  for (int g = 0; g < GradVV.Len(); g++) {
    GPAll.AddPlot(GradVV[g], gpwLines, TStr::Fmt("param %d", g+1), "linewidth 1");
    GetMinMax(GradVV[g], Min1, Max1, false);
    TGnuPlot GP(TStr::Fmt("sGrad%02d-", g+1)+OutFNm, TStr::Fmt("Gradient (avg 1k samples). %s", Desc.CStr()), true);
    GP.AddPlot(GradVV[g], gpwLines, TStr::Fmt("param id %d", g+1), "linewidth 1");
    GetMinMax(GradVV[g], Min, Max, true);
    GP.SetYRange((int)floor(Min-1), (int)ceil(Max+1));
    GP.SetXYLabel("Sample Index (time)", "Gradient");
    GP.SavePng();
  }
  GPAll.SetYRange((int)floor(Min1-1), (int)ceil(Max1+1));
  GPAll.SetXYLabel("Sample Index (time)", "Gradient");
  GPAll.SavePng();
}

void PlotAutoCorrelation(const TFltV& ValV, const int& MaxK, const TStr& OutFNm, const TStr& Desc) {
  double Avg=0.0, Var=0.0;
  for (int i = 0; i < ValV.Len(); i++) { Avg += ValV[i]; }
  Avg /= (double) ValV.Len();
  for (int i = 0; i < ValV.Len(); i++) { Var += TMath::Sqr(ValV[i]-Avg); }
  TFltPrV ACorrV;
  for (int k = 0; k < TMath::Mn(ValV.Len(), MaxK); k++) {
    double corr = 0.0;
    for (int i = 0; i < ValV.Len() - k; i++) {
      corr += (ValV[i]-Avg)*(ValV[i+k]-Avg);
    }
    ACorrV.Add(TFltPr(k, corr/Var));
  }
  // plot grads
  TGnuPlot GP("sAutoCorr-"+OutFNm, TStr::Fmt("AutoCorrelation (%d samples). %s", ValV.Len(), Desc.CStr()), true);
  GP.AddPlot(ACorrV, gpwLines, "", "linewidth 1");
  GP.SetXYLabel("Lag, k", "Autocorrelation, r_k");
  GP.SavePng();
}

// sample permutations and plot the current gradient and log-likelihood as the function
// of the number of samples
TFltV TKroneckerLL::TestSamplePerm(const TStr& OutFNm, const int& WarmUp, const int& NSamples, const TKronMtx& TrueMtx, const bool& DoPlot) {
  printf("Sample permutations: %s (warm-up: %s)\n", TInt::GetMegaStr(NSamples).CStr(), TInt::GetMegaStr(WarmUp).CStr());
  int NId1=0, NId2=0, NAccept=0;
  TExeTm ExeTm;
  const int PlotLen = NSamples/1000+1;
  double TrueLL=-1, AvgLL=0.0;
  TFltV AvgGradV(GetParams());
  TFltPrV TrueLLV(PlotLen, 0); // true log-likelihood (under the correct permutation)
  TFltPrV EstLLV(PlotLen, 0);  // estiamted log-likelihood (averaged over last 1k permutation)
  TFltPrV AcceptV;             // sample acceptance ratio
  TFltV SampleLLV(NSamples, 0);
  TVec<TFltPrV> GradVV(GetParams());
  for (int g = 0; g < GetParams(); g++) {
    GradVV[g].Gen(PlotLen, 0); }
  if (! TrueMtx.Empty()) {
    TIntV PermV=NodePerm;  TKronMtx CurMtx=ProbMtx;  ProbMtx.Dump();
    InitLL(TrueMtx);  SetOrderPerm();  CalcApxGraphLL();  printf("TrueLL: %f\n", LogLike());
    TrueLL=LogLike;  InitLL(CurMtx); NodePerm=PermV;
  }
  CalcApxGraphLL();
  printf("LogLike at start:       %f\n", LogLike());
  if (WarmUp > 0) {
    EstLLV.Add(TFltPr(0, LogLike));
    if (TrueLL != -1) { TrueLLV.Add(TFltPr(0, TrueLL)); }
    for (int s = 0; s < WarmUp; s++) { SampleNextPerm(NId1, NId2); }
    printf("  warm-up:%s,", ExeTm.GetTmStr());  ExeTm.Tick();
  }
  printf("LogLike afterm warm-up: %f\n", LogLike());
  CalcApxGraphLL(); // re-calculate LL (due to numerical errors)
  CalcApxGraphDLL();
  EstLLV.Add(TFltPr(WarmUp, LogLike));
  if (TrueLL != -1) { TrueLLV.Add(TFltPr(WarmUp, TrueLL)); }
  printf("  recalculated:         %f\n", LogLike());
  // start sampling
  printf("  sampling (average per 1000 samples)\n");
  TVec<TFltV> SamplVV(5);
  for (int s = 0; s < NSamples; s++) {
    if (SampleNextPerm(NId1, NId2)) { // new permutation
      UpdateGraphDLL(NId1, NId2);  NAccept++; }
    for (int m = 0; m < AvgGradV.Len(); m++) { AvgGradV[m] += GradV[m]; }
    AvgLL += GetLL();
    SampleLLV.Add(GetLL());
    /*SamplVV[0].Add(GetLL()); // gives worse autocoreelation than the avg below
    SamplVV[1].Add(GradV[0]);
    SamplVV[2].Add(GradV[1]);
    SamplVV[3].Add(GradV[2]);
    SamplVV[4].Add(GradV[3]);*/
    if (s > 0 && s % 1000 == 0) {
      printf(".");
      for (int g = 0; g < AvgGradV.Len(); g++) {
        GradVV[g].Add(TFltPr(WarmUp+s, AvgGradV[g] / 1000.0)); }
      EstLLV.Add(TFltPr(WarmUp+s, AvgLL / 1000.0));
      if (TrueLL != -1) { TrueLLV.Add(TFltPr(WarmUp+s, TrueLL)); }
      AcceptV.Add(TFltPr(WarmUp+s, NAccept/1000.0));
      // better (faster decaying) autocorrelation when one takes avg. of 1000 consecutive samples
      /*SamplVV[0].Add(AvgLL);
      SamplVV[1].Add(AvgGradV[0]);
      SamplVV[2].Add(AvgGradV[1]);
      SamplVV[3].Add(AvgGradV[2]);
      SamplVV[4].Add(AvgGradV[3]); //*/
      if (s % 100000 == 0 && DoPlot) {
        const TStr Desc = TStr::Fmt("P(NodeSwap)=%g. Nodes: %d, Edges: %d, Params: %d, WarmUp: %s, Samples: %s", PermSwapNodeProb(),
          Graph->GetNodes(), Graph->GetEdges(), GetParams(), TInt::GetMegaStr(WarmUp).CStr(), TInt::GetMegaStr(NSamples).CStr());
        PlotGrad(EstLLV, TrueLLV, GradVV, AcceptV, OutFNm, Desc);
        for (int n = 0; n < SamplVV.Len(); n++) {
          PlotAutoCorrelation(SamplVV[n], 1000, TStr::Fmt("%s-n%d", OutFNm.CStr(), n), Desc); }
        printf("  samples: %d, time: %s, samples/s: %.1f\n", s, ExeTm.GetTmStr(), double(s+1)/ExeTm.GetSecs());
      }
      AvgLL = 0;  AvgGradV.PutAll(0);  NAccept=0;
    }
  }
  if (DoPlot) {
    const TStr Desc = TStr::Fmt("P(NodeSwap)=%g. Nodes: %d, Edges: %d, Params: %d, WarmUp: %s, Samples: %s", PermSwapNodeProb(),
      Graph->GetNodes(), Graph->GetEdges(), GetParams(), TInt::GetMegaStr(WarmUp).CStr(), TInt::GetMegaStr(NSamples).CStr());
    PlotGrad(EstLLV, TrueLLV, GradVV, AcceptV, OutFNm, Desc);
    for (int n = 0; n < SamplVV.Len(); n++) {
      PlotAutoCorrelation(SamplVV[n], 1000, TStr::Fmt("%s-n%d", OutFNm.CStr(), n), Desc); }
  }
  return SampleLLV; // seems to work better for potential scale reduction plot
}

void McMcGetAvgAvg(const TFltV& AvgJV, double& AvgAvg) {
  AvgAvg = 0.0;
  for (int j = 0; j < AvgJV.Len(); j++) {
    AvgAvg += AvgJV[j]; }
  AvgAvg /= AvgJV.Len();
}

void McMcGetAvgJ(const TVec<TFltV>& ChainLLV, TFltV& AvgJV) {
  for (int j = 0; j < ChainLLV.Len(); j++) {
    const TFltV& ChainV = ChainLLV[j];
    double Avg = 0;
    for (int i = 0; i < ChainV.Len(); i++) {
      Avg += ChainV[i];
    }
    AvgJV.Add(Avg/ChainV.Len());
  }
}

// calculates the chain potential scale reduction (see Gelman Bayesian Statistics book)
double TKroneckerLL::CalcChainR2(const TVec<TFltV>& ChainLLV) {
  const double J = ChainLLV.Len();
  const double K = ChainLLV[0].Len();
  TFltV AvgJV;    McMcGetAvgJ(ChainLLV, AvgJV);
  double AvgAvg;  McMcGetAvgAvg(AvgJV, AvgAvg);
  IAssert(AvgJV.Len() == ChainLLV.Len());
  double InChainVar=0, OutChainVar=0;
  // between chain var
  for (int j = 0; j < AvgJV.Len(); j++) {
    OutChainVar += TMath::Sqr(AvgJV[j] - AvgAvg); }
  OutChainVar = OutChainVar * (K/double(J-1));
  printf("*** %g chains of len %g\n", J, K);
  printf("  ** between chain var: %f\n", OutChainVar);
  //within chain variance
  for (int j = 0; j < AvgJV.Len(); j++) {
    const TFltV& ChainV = ChainLLV[j];
    for (int k = 0; k < ChainV.Len(); k++) {
      InChainVar += TMath::Sqr(ChainV[k] - AvgJV[j]); }
  }
  InChainVar = InChainVar * 1.0/double(J*(K-1));
  printf("  ** within chain var: %f\n", InChainVar);
  const double PostVar = (K-1)/K * InChainVar + 1.0/K * OutChainVar;
  printf("  ** posterior var: %f\n", PostVar);
  const double ScaleRed = sqrt(PostVar/InChainVar);
  printf("  ** scale reduction (< 1.2): %f\n\n", ScaleRed);
  return ScaleRed;
}

// Gelman-Rubin-Brooks plot: how does potential scale reduction chainge with chain length
void TKroneckerLL::ChainGelmapRubinPlot(const TVec<TFltV>& ChainLLV, const TStr& OutFNm, const TStr& Desc) {
  TFltPrV LenR2V; // how does potential scale reduction chainge with chain length
  TVec<TFltV> SmallLLV(ChainLLV.Len());
  const int K = ChainLLV[0].Len();
  const int Buckets=1000;
  const int BucketSz = K/Buckets;
  for (int b = 1; b < Buckets; b++) {
    const int End = TMath::Mn(BucketSz*b, K-1);
    for (int c = 0; c < ChainLLV.Len(); c++) {
      ChainLLV[c].GetSubValV(0, End, SmallLLV[c]); }
    LenR2V.Add(TFltPr(End, TKroneckerLL::CalcChainR2(SmallLLV)));
  }
  LenR2V.Add(TFltPr(K, TKroneckerLL::CalcChainR2(ChainLLV)));
  TGnuPlot::PlotValV(LenR2V, TStr::Fmt("gelman-%s", OutFNm.CStr()), TStr::Fmt("%s. %d chains of len %d. BucketSz: %d.",
    Desc.CStr(), ChainLLV.Len(), ChainLLV[0].Len(), BucketSz), "Chain length", "Potential scale reduction");
}

// given a Kronecker graph generate from TrueParam, try to recover the parameters
TFltQu TKroneckerLL::TestKronDescent(const bool& DoExact, const bool& TruePerm, double LearnRate, const int& WarmUp, const int& NSamples, const TKronMtx& TrueParam) {
  printf("Test gradient descent on a synthetic kronecker graphs:\n");
  if (DoExact) { printf("  -- Exact gradient calculations\n"); }
  else { printf("  -- Approximate gradient calculations\n"); }
  if (TruePerm) { printf("  -- No permutation sampling (use true permutation)\n"); }
  else { printf("  -- Sample permutations (start with degree permutation)\n"); }
  TExeTm IterTm;
  int Iter;
  double OldLL=0, MyLL=0, AvgAbsErr, AbsSumErr;
  TFltV MyGradV, SDevV;
  TFltV LearnRateV(GetParams());  LearnRateV.PutAll(LearnRate);
  if (TruePerm) {
    SetOrderPerm();
  }
  else {
    /*printf("SET EIGEN VECTOR PERMUTATIONS\n");
    TFltV LeftSV, RightSV;
    TGSvd::GetSngVec(Graph, LeftSV, RightSV);
    TFltIntPrV V;
    for (int v=0; v<LeftSV.Len();v++) { V.Add(TFltIntPr(LeftSV[v], v)); }
    V.Sort(false);
    NodePerm.Gen(Nodes, 0);
    for (int v=0; v < V.Len();v++) { NodePerm.Add(V[v].Val2); } //*/
    //printf("RANDOM PERMUTATION\n");   SetRndPerm();
    printf("DEGREE  PERMUTATION\n");  SetDegPerm();
  }
  for (Iter = 0; Iter < 100; Iter++) {
    if (TruePerm) {
      // don't sample over permutations
      if (DoExact) { CalcGraphDLL();  CalcGraphLL(); } // slow, O(N^2)
      else { CalcApxGraphDLL();  CalcApxGraphLL(); }   // fast
      MyLL = LogLike;  MyGradV = GradV;
    } else {
      printf(".");
      // sample over permutations (approximate calculations)
      SampleGradient(WarmUp, NSamples, MyLL, MyGradV);
    }
    printf("%d] LL: %g, ", Iter, MyLL);
    AvgAbsErr = TKronMtx::GetAvgAbsErr(ProbMtx, TrueParam);
    AbsSumErr = fabs(ProbMtx.GetMtxSum() - TrueParam.GetMtxSum());
    printf("  avgAbsErr: %.4f, absSumErr: %.4f, newLL: %.2f, deltaLL: %.2f\n", AvgAbsErr, AbsSumErr, MyLL, OldLL-MyLL);
    for (int p = 0; p < GetParams(); p++) {
      // set learn rate so that move for each parameter is inside the [0.01, 0.1]
      LearnRateV[p] *= 0.9;
      //printf("%d: rate: %f   delta:%f\n", p, LearnRateV[p], fabs(LearnRateV[p]*MyGradV[p]));
      while (fabs(LearnRateV[p]*MyGradV[p]) > 0.1) { LearnRateV[p] *= 0.9; }
      //printf("   rate: %f   delta:%f\n", LearnRateV[p], fabs(LearnRateV[p]*MyGradV[p]));
      while (fabs(LearnRateV[p]*MyGradV[p]) < 0.001) { LearnRateV[p] *= (1.0/0.9); }
      //printf("   rate: %f   delta:%f\n", LearnRateV[p], fabs(LearnRateV[p]*MyGradV[p]));
      printf("    %d]  %f  <--  %f + %f    lrnRate:%g\n", p, ProbMtx.At(p) + LearnRateV[p]*MyGradV[p],
        ProbMtx.At(p), (double)(LearnRateV[p]*MyGradV[p]), LearnRateV[p]());
      ProbMtx.At(p) = ProbMtx.At(p) + LearnRateV[p]*MyGradV[p];
      // box constraints
      if (ProbMtx.At(p) > 0.99) { ProbMtx.At(p)=0.99; }
      if (ProbMtx.At(p) < 0.01) { ProbMtx.At(p)=0.01; }
    }
    ProbMtx.GetLLMtx(LLMtx);  OldLL = MyLL;
    if (AvgAbsErr < 0.01) { printf("***CONVERGED!\n");  break; }
    printf("\n");  fflush(stdout);
  }
  TrueParam.Dump("True  Thetas", true);
  ProbMtx.Dump("Final Thetas", true);
  printf("  AvgAbsErr: %f\n  AbsSumErr: %f\n  Iterations: %d\n", AvgAbsErr, AbsSumErr, Iter);
  printf("Iteration run time: %s, sec: %g\n\n", IterTm.GetTmStr(), IterTm.GetSecs());
  return TFltQu(AvgAbsErr, AbsSumErr, Iter, IterTm.GetSecs());
}

void PlotTrueAndEst(const TStr& OutFNm, const TStr& Desc, const TStr& YLabel, const TFltPrV& EstV, const TFltPrV& TrueV) {
  TGnuPlot GP(OutFNm, Desc.CStr(), true);
  GP.AddPlot(EstV, gpwLinesPoints, YLabel, "linewidth 1 pointtype 6 pointsize 1");
  if (! TrueV.Empty()) { GP.AddPlot(TrueV, gpwLines, "TRUE"); }
  GP.SetXYLabel("Gradient descent iterations", YLabel);
  GP.SavePng();
}

void TKroneckerLL::GradDescentConvergence(const TStr& OutFNm, const TStr& Desc1, const bool& SamplePerm, const int& NIters,
 double LearnRate, const int& WarmUp, const int& NSamples, const int& AvgKGraphs, const TKronMtx& TrueParam) {
  TExeTm IterTm;
  int Iter;
  double OldLL=0, MyLL=0, AvgAbsErr=0, AbsSumErr=0;
  TFltV MyGradV, SDevV;
  TFltV LearnRateV(GetParams());  LearnRateV.PutAll(LearnRate);
  TFltPrV EZeroV, DiamV, Lambda1V, Lambda2V, AvgAbsErrV, AvgLLV;
  TFltPrV TrueEZeroV, TrueDiamV, TrueLambda1V, TrueLambda2V, TrueLLV;
  TFltV SngValV;  TSnap::GetSngVals(Graph, 2, SngValV);  SngValV.Sort(false);
  const double TrueEZero = pow((double) Graph->GetEdges(), 1.0/double(KronIters));
  const double TrueEffDiam = TSnap::GetAnfEffDiam(Graph, false, 10);
  const double TrueLambda1 = SngValV[0];
  const double TrueLambda2 = SngValV[1];
  if (! TrueParam.Empty()) {
    const TKronMtx CurParam = ProbMtx;  ProbMtx.Dump();
    InitLL(TrueParam);  SetOrderPerm();  CalcApxGraphLL(); printf("TrueLL: %f\n", LogLike());
    OldLL = LogLike;  InitLL(CurParam);
  }
  const double TrueLL = OldLL;
  if (! SamplePerm) { SetOrderPerm(); } else { SetDegPerm(); }
  for (Iter = 0; Iter < NIters; Iter++) {
    if (! SamplePerm) {
      // don't sample over permutations
      CalcApxGraphDLL();  CalcApxGraphLL();   // fast
      MyLL = LogLike;  MyGradV = GradV;
    } else {
      // sample over permutations (approximate calculations)
      SampleGradient(WarmUp, NSamples, MyLL, MyGradV);
    }
    double SumDiam=0, SumSngVal1=0, SumSngVal2=0;
    for (int trial = 0; trial < AvgKGraphs; trial++) {
      // generate kronecker graph
      PNGraph KronGraph = TKronMtx::GenFastKronecker(ProbMtx, KronIters, true, 0); // approx
      //PNGraph KronGraph = TKronMtx::GenKronecker(ProbMtx, KronIters, true, 0); // true
      SngValV.Clr(true);  TSnap::GetSngVals(KronGraph, 2, SngValV);  SngValV.Sort(false);
      SumDiam += TSnap::GetAnfEffDiam(KronGraph, false, 10);
      SumSngVal1 += SngValV[0];  SumSngVal2 += SngValV[1];
    }
    // how good is the current fit
    AvgLLV.Add(TFltPr(Iter, MyLL));
    EZeroV.Add(TFltPr(Iter, ProbMtx.GetMtxSum()));
    DiamV.Add(TFltPr(Iter, SumDiam/double(AvgKGraphs)));
    Lambda1V.Add(TFltPr(Iter, SumSngVal1/double(AvgKGraphs)));
    Lambda2V.Add(TFltPr(Iter, SumSngVal2/double(AvgKGraphs)));
    TrueLLV.Add(TFltPr(Iter, TrueLL));
    TrueEZeroV.Add(TFltPr(Iter, TrueEZero));
    TrueDiamV.Add(TFltPr(Iter, TrueEffDiam));
    TrueLambda1V.Add(TFltPr(Iter, TrueLambda1));
    TrueLambda2V.Add(TFltPr(Iter, TrueLambda2));
    if (Iter % 10 == 0) {
      const TStr Desc = TStr::Fmt("%s. Iter: %d, G(%d, %d)  K(%d, %d)", Desc1.Empty()?OutFNm.CStr():Desc1.CStr(),
        Iter, Graph->GetNodes(), Graph->GetEdges(), ProbMtx.GetNodes(KronIters), ProbMtx.GetEdges(KronIters));
      PlotTrueAndEst("LL."+OutFNm, Desc, "Average LL", AvgLLV, TrueLLV);
      PlotTrueAndEst("E0."+OutFNm, Desc, "E0 (expected number of edges)", EZeroV, TrueEZeroV);
      PlotTrueAndEst("Diam."+OutFNm+"-Diam", Desc, "Effective diameter", DiamV, TrueDiamV);
      PlotTrueAndEst("Lambda1."+OutFNm, Desc, "Lambda 1", Lambda1V, TrueLambda1V);
      PlotTrueAndEst("Lambda2."+OutFNm, Desc, "Lambda 2", Lambda2V, TrueLambda2V);
      if (! TrueParam.Empty()) {
        PlotTrueAndEst("AbsErr."+OutFNm, Desc, "Average Absolute Error", AvgAbsErrV, TFltPrV()); }
    }
    if (! TrueParam.Empty()) {
      AvgAbsErr = TKronMtx::GetAvgAbsErr(ProbMtx, TrueParam);
      AvgAbsErrV.Add(TFltPr(Iter, AvgAbsErr));
    } else { AvgAbsErr = 1.0; }
    // update parameters
    AbsSumErr = fabs(ProbMtx.GetMtxSum() - TrueEZero);
    // update parameters
    for (int p = 0; p < GetParams(); p++) {
      LearnRateV[p] *= 0.99;
      while (fabs(LearnRateV[p]*MyGradV[p]) > 0.1) { LearnRateV[p] *= 0.99; printf(".");}
      while (fabs(LearnRateV[p]*MyGradV[p]) < 0.002) { LearnRateV[p] *= (1.0/0.95); printf("*");}
      printf("    %d]  %f  <--  %f + %9f   Grad: %9.1f,  Rate:%g\n", p, ProbMtx.At(p) + LearnRateV[p]*MyGradV[p],
        ProbMtx.At(p), (double)(LearnRateV[p]*MyGradV[p]), MyGradV[p](), LearnRateV[p]());
      ProbMtx.At(p) = ProbMtx.At(p) + LearnRateV[p]*MyGradV[p];
      // box constraints
      if (ProbMtx.At(p) > 1.0) { ProbMtx.At(p)=1.0; }
      if (ProbMtx.At(p) < 0.001) { ProbMtx.At(p)=0.001; }
    }
    printf("%d] LL: %g, ", Iter, MyLL);
    printf("  avgAbsErr: %.4f, absSumErr: %.4f, newLL: %.2f, deltaLL: %.2f\n", AvgAbsErr, AbsSumErr, MyLL, OldLL-MyLL);
    if (AvgAbsErr < 0.001) { printf("***CONVERGED!\n");  break; }
    printf("\n");  fflush(stdout);
    ProbMtx.GetLLMtx(LLMtx);  OldLL = MyLL;
  }
  TrueParam.Dump("True  Thetas", true);
  ProbMtx.Dump("Final Thetas", true);
  printf("  AvgAbsErr: %f\n  AbsSumErr: %f\n  Iterations: %d\n", AvgAbsErr, AbsSumErr, Iter);
  printf("Iteration run time: %s, sec: %g\n\n", IterTm.GetTmStr(), IterTm.GetSecs());
}

// given true N0, fit the parameters, get likelihood and calculate BIC (MDL), plot n0 vs. BIC
void TKroneckerLL::TestBicCriterion(const TStr& OutFNm, const TStr& Desc1, const PNGraph& G, const int& GradIters,
 double LearnRate, const int& WarmUp, const int& NSamples, const int& TrueN0) {
  TFltPrV BicV, MdlV, LLV;
  const double rndGP = G->GetEdges()/TMath::Sqr(double(G->GetNodes()));
  const double RndGLL = G->GetEdges()*log(rndGP )+ (TMath::Sqr(double(G->GetNodes()))-G->GetEdges())*log(1-rndGP);
  LLV.Add(TFltPr(1, RndGLL));
  BicV.Add(TFltPr(1, -RndGLL + 0.5*TMath::Sqr(1)*log(TMath::Sqr(G->GetNodes()))));
  MdlV.Add(TFltPr(1, -RndGLL + 32*TMath::Sqr(1)+2*(log((double)1)+log((double)G->GetNodes()))));
  for (int NZero = 2; NZero < 10; NZero++) {
    const TKronMtx InitKronMtx = TKronMtx::GetInitMtx(NZero, G->GetNodes(), G->GetEdges());
    InitKronMtx.Dump("INIT PARAM", true);
    TKroneckerLL KronLL(G, InitKronMtx);
    KronLL.SetPerm('d'); // degree perm
    const double LastLL = KronLL.GradDescent(GradIters, LearnRate, 0.001, 0.01, WarmUp, NSamples);
    LLV.Add(TFltPr(NZero, LastLL));
    BicV.Add(TFltPr(NZero, -LastLL + 0.5*TMath::Sqr(NZero)*log(TMath::Sqr(G->GetNodes()))));
    MdlV.Add(TFltPr(NZero, -LastLL + 32*TMath::Sqr(NZero)+2*(log((double)NZero)+log((double)KronLL.GetKronIters()))));
    { TGnuPlot GP("LL-"+OutFNm, Desc1);
    GP.AddPlot(LLV, gpwLinesPoints, "Log-likelihood", "linewidth 1 pointtype 6 pointsize 2");
    GP.SetXYLabel("NZero", "Log-Likelihood");  GP.SavePng(); }
    { TGnuPlot GP("BIC-"+OutFNm, Desc1);
    GP.AddPlot(BicV, gpwLinesPoints, "BIC", "linewidth 1 pointtype 6 pointsize 2");
    GP.SetXYLabel("NZero", "BIC");  GP.SavePng(); }
    { TGnuPlot GP("MDL-"+OutFNm, Desc1);
    GP.AddPlot(MdlV, gpwLinesPoints, "MDL", "linewidth 1 pointtype 6 pointsize 2");
    GP.SetXYLabel("NZero", "MDL");  GP.SavePng(); }
  }
}

void TKroneckerLL::TestGradDescent(const int& KronIters, const int& KiloSamples, const TStr& Permutation) {
  const TStr OutFNm = TStr::Fmt("grad-%s%d-%dk", Permutation.CStr(), KronIters, KiloSamples);
  TKronMtx KronParam = TKronMtx::GetMtx("0.8 0.6; 0.6 0.4");
  PNGraph Graph  = TKronMtx::GenFastKronecker(KronParam, KronIters, true, 0);
  TKroneckerLL KronLL(Graph, KronParam);
  TVec<TFltV> GradVV(4), SDevVV(4);  TFltV XValV;
  int NId1 = 0, NId2 = 0, NAccept = 0;
  TVec<TMom> GradMomV(4);
  TExeTm ExeTm;
  if (Permutation == "r") KronLL.SetRndPerm();
  else if (Permutation == "d") KronLL.SetDegPerm();
  else if (Permutation == "o") KronLL.SetOrderPerm();
  else FailR("Unknown permutation (r,d,o)");
  KronLL.CalcApxGraphLL();
  KronLL.CalcApxGraphDLL();
  for (int s = 0; s < 1000*KiloSamples; s++) {
    if (KronLL.SampleNextPerm(NId1, NId2)) { // new permutation
      KronLL.UpdateGraphDLL(NId1, NId2);  NAccept++; }
    if (s > 50000) { //warm up period
      for (int m = 0; m < 4; m++) { GradVV[m].Add(KronLL.GradV[m]); }
      if ((s+1) % 1000 == 0) {
        printf(".");
        for (int m = 0; m < 4; m++) { GradVV[m].Add(KronLL.GradV[m]); }
        XValV.Add((s+1));
        if ((s+1) % 100000 == 0) {
          TGnuPlot GP(OutFNm, TStr::Fmt("Gradient vs. samples. %d nodes, %d edges", Graph->GetNodes(), Graph->GetEdges()), true);
          for (int g = 0; g < GradVV.Len(); g++) {
            GP.AddPlot(XValV, GradVV[g], gpwLines, TStr::Fmt("grad %d", g)); }
          GP.SetXYLabel("sample index","log Gradient");
          GP.SavePng();
        }
      }
    }
  }
  printf("\n");
  for (int m = 0; m < 4; m++) {
    GradMomV[m].Def();
    printf("grad %d: mean: %12f  sDev: %12f  median: %12f\n", m,
      GradMomV[m].GetMean(), GradMomV[m].GetSDev(), GradMomV[m].GetMedian());
  }
}
/*
// sample over permutations
void TKroneckerLL::GradDescent(const double& LearnRate, const int& WarmUp, const int& NSamples, const int& NIter) {
  TFltV GradV, SDevV;
  double AvgLL;
  for (int Iter = 0; Iter < 100; Iter++) {
    //SampleGradient(WarmUp, NSamples, AvgLL, GradV, SDevV, true);
    SampleGradient(WarmUp, NSamples, AvgLL, GradV);
    for (int theta = 0; theta < GetParams(); theta++) {
      printf("%d]  %f  <--  %f + %f\n", theta, ProbMtx.At(theta) + LearnRate*GradV[theta], ProbMtx.At(theta), LearnRate*GradV[theta]);
      ProbMtx.At(theta) = ProbMtx.At(theta) + LearnRate*GradV[theta];
      // box constraints
      if (ProbMtx.At(theta) > 0.99) ProbMtx.At(theta)=0.99;
      if (ProbMtx.At(theta) < 0.01) ProbMtx.At(theta)=0.01;
    }
    ProbMtx.GetLLMtx(LLMtx);
  }
  ProbMtx.Dump("Final Thetas");
}
*/


/////////////////////////////////////////////////
// Add Noise to Graph
/// !!!!! MYUNGHWAN, CHECK!
int TKronNoise::RemoveNodeNoise(PNGraph& Graph, const int& NNodes, const bool Random) {
	IAssert(NNodes > 0 && NNodes < (Graph->GetNodes() / 2));

	int i = 0;
	TIntV ShufflePerm;
	Graph->GetNIdV(ShufflePerm);
	if(Random) {
		ShufflePerm.Shuffle(TKronMtx::Rnd);
		for(i = 0; i < NNodes; i++) {
			Graph->DelNode(int(ShufflePerm[i]));
		}
	} else {
		for(i = 0; i < NNodes; i++) {
			Graph->DelNode(int(ShufflePerm[ShufflePerm.Len() - 1 - i]));
		}
	}

	return Graph->GetNodes();
}

int TKronNoise::RemoveNodeNoise(PNGraph& Graph, const double& Rate, const bool Random) {
	IAssert(Rate > 0 && Rate < 0.5);
	return TKronNoise::RemoveNodeNoise(Graph, (int) floor(Rate * double(Graph->GetNodes())), Random);
}

int TKronNoise::FlipEdgeNoise(PNGraph& Graph, const int& NEdges, const bool Random) {
	IAssert(NEdges > 0 && NEdges < Graph->GetEdges());

	const int Nodes = Graph->GetNodes();
	const int Edges = Graph->GetEdges();
	int Src, Dst;

	TIntV NIdV, TempV;
	TIntPrV ToAdd, ToDel;
	Graph->GetNIdV(NIdV);

	ToAdd.Gen(NEdges / 2, 0);
	for(int i = 0; i < NEdges / 2; i++) {
		Src = NIdV[TKronMtx::Rnd.GetUniDevInt(Nodes)];
		Dst = NIdV[TKronMtx::Rnd.GetUniDevInt(Nodes)];
		if(Graph->IsEdge(Src, Dst)) {	i--;	continue;	}

		ToAdd.Add(TIntPr(Src, Dst));
	}

	ToDel.Gen(Edges, 0);
	for(TNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
		ToDel.Add(TIntPr(EI.GetSrcNId(), EI.GetDstNId()));
	}
	ToDel.Shuffle(TKronMtx::Rnd);

	for(int i = 0; i < NEdges / 2; i++) {
		Graph->DelEdge(ToDel[i].Val1, ToDel[i].Val2);
		Graph->AddEdge(ToAdd[i].Val1, ToAdd[i].Val2);
	}

	return Graph->GetEdges();
}

int TKronNoise::FlipEdgeNoise(PNGraph& Graph, const double& Rate, const bool Random) {
	IAssert(Rate > 0 && Rate < 0.5);
	return TKronNoise::FlipEdgeNoise(Graph, (int) floor(Rate * double(Graph->GetEdges())), Random);
}

int TKronNoise::RemoveEdgeNoise(PNGraph& Graph, const int& NEdges) {
	IAssert(NEdges > 0 && NEdges < Graph->GetEdges());

	TIntPrV ToDel;

	ToDel.Gen(Graph->GetEdges(), 0);
	for(TNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
		if(EI.GetSrcNId() != EI.GetDstNId()) {
			ToDel.Add(TIntPr(EI.GetSrcNId(), EI.GetDstNId()));
		}
	}
	ToDel.Shuffle(TKronMtx::Rnd);

	for(int i = 0; i < NEdges; i++) {
		Graph->DelEdge(ToDel[i].Val1, ToDel[i].Val2);
	}

	return Graph->GetEdges();
}

int TKronNoise::RemoveEdgeNoise(PNGraph& Graph, const double& Rate) {
	IAssert(Rate > 0 && Rate < 0.5);
	return TKronNoise::RemoveEdgeNoise(Graph, (int) floor(Rate * double(Graph->GetEdges())));
}



/////////////////////////////////////////////////
// Kronecker Log Likelihood Maximization
void TKronMaxLL::SetPerm(const char& PermId) {
  if (PermId == 'o') KronLL.SetOrderPerm();
  else if (PermId == 'd') KronLL.SetDegPerm();
  else if (PermId == 'r') KronLL.SetRndPerm();
  else FailR("Unknown permutation type (o,d,r)");
}

/* void TKronMaxLL::EvalNewEdgeP(const TKronMtx& ProbMtx) {
  ProbMtx.Dump("\n***EVAL:");
  for (int i = 0; i < ProbMtx.Len(); i++) {
    // parameters are out of bounds
    if (ProbMtx.At(i) < 0.05 || ProbMtx.At(i) > 0.95) {
      TFltV MxDerivV(ProbMtx.Len()); MxDerivV.PutAll(1e5);
      FEvalH.AddDat(ProbMtx, TFEval(-1e10, MxDerivV));
      return;
    }
  }
  double AvgLL;
  TFltV GradV, SDevV;
  KronLL.InitLL(ProbMtx); // set current point
  //KronLL.SampleGradient(WarmUp, NSamples, AvgLL, GradV, SDevV, true); //sample the gradient
  KronLL.SampleGradient(WarmUp, NSamples, AvgLL, GradV);
  FEvalH.AddDat(ProbMtx, TFEval(AvgLL, GradV));
}

double TKronMaxLL::GetLL(const TFltV& ThetaV) {
  TKronMtx ProbMtx;  RoundTheta(ThetaV, ProbMtx);
  if (! FEvalH.IsKey(ProbMtx)) {
    EvalNewEdgeP(ProbMtx);
  }
  return FEvalH.GetDat(ProbMtx).LogLike;
}

void TKronMaxLL::GetDLL(const TFltV& ThetaV, TFltV& GradV) {
  TKronMtx ProbMtx;  RoundTheta(ThetaV, ProbMtx);
  if (! FEvalH.IsKey(ProbMtx)) {
    EvalNewEdgeP(ProbMtx);
  }
  GradV = FEvalH.GetDat(ProbMtx).GradV;
}

double TKronMaxLL::GetDLL(const TFltV& ThetaV, const int& ParamId) {
  TKronMtx ProbMtx;  RoundTheta(ThetaV, ProbMtx);
  if (! FEvalH.IsKey(ProbMtx)) {
    EvalNewEdgeP(ProbMtx);
  }
  return FEvalH.GetDat(ProbMtx).GradV[ParamId];
}
void TKronMaxLL::MaximizeLL(const int& NWarmUp, const int& Samples) {
  WarmUp = NWarmUp;
  NSamples = Samples;
  TConjGrad<TFunc> ConjGrad(KronLL.GetProbMtx().GetMtx(), TFunc(this));
  //TConjGrad<TLogBarFunc> ConjGrad(KronLL.GetEdgeP().GetV(), TLogBarFunc(this, 0.1));
  ConjGrad.ConjGradMin(0.1);
}*/

// round to 3 decimal places
void TKronMaxLL::RoundTheta(const TFltV& ThetaV, TFltV& NewThetaV) {
  NewThetaV.Gen(ThetaV.Len());
  for (int i = 0; i < ThetaV.Len(); i++) {
    NewThetaV[i] = TMath::Round(ThetaV[i], 3); }
}

// round to 3 decimal places
void TKronMaxLL::RoundTheta(const TFltV& ThetaV, TKronMtx& Kronecker) {
  Kronecker.GenMtx((int)sqrt((double)ThetaV.Len()));
  for (int i = 0; i < ThetaV.Len(); i++) {
    Kronecker.At(i) = TMath::Round(ThetaV[i], 3); }
}

void TKronMaxLL::Test() {
  TKronMtx::PutRndSeed(1);
  TKronMtx KronParam = TKronMtx::GetMtx("0.8 0.7; 0.6 0.5");
  PNGraph Graph  = TKronMtx::GenFastKronecker(KronParam, 8, true, 1);

  TKronMaxLL KronMaxLL(Graph, TFltV::GetV(0.9, 0.7, 0.5, 0.3));
  KronMaxLL.SetPerm('d');
  //KronMaxLL.MaximizeLL(10000, 50000);
  /*TKroneckerLL KronLL(Graph, *TKronMtx::GetMtx("0.9 0.7; 0.5 0.3"));
  KronLL.SetDegPerm();
  KronLL.GradDescent(0.005/double(Graph->GetNodes()));*/
}

/////////////////////////////////////////////////
// Kronecker Phase Plot
/*
void TKronPhasePlot::SaveMatlab(const TStr& OutFNm) const {
  FILE *F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "#Take last graph stats\n");
  fprintf(F, "#i\tAlpha\tBeta\tNodes\tNonZNodes\tEdges\tWccNodes\tWccEdges\tDiam\tEffDiam\tWccDiam\tWccEffDiam\t1StEigVal\tMxEigVal\n");
  for (int i = 0 ; i < PhaseV.Len(); i++) {
    const TPhasePoint& Point = PhaseV[i];
    const TGrowthStat& GrowthStat = Point.GrowthStat;
    const PGraphStat& FirstGrowth = GrowthStat[0];
    const PGraphStat& LastGrowth = GrowthStat.Last();
    fprintf(F, "%d\t%g\t%g\t", i, Point.Alpha, Point.Beta);
    fprintf(F, "%d\t%d\t%d\t", LastGrowth->Nodes, LastGrowth->Edges, LastGrowth->NonZNodes);
    fprintf(F, "%d\t%d\t", LastGrowth->WccNodes, LastGrowth->WccEdges);
    fprintf(F, "%f\t%f\t%f\t%f\t", LastGrowth->FullDiam, LastGrowth->EffDiam, LastGrowth->FullWccDiam, LastGrowth->EffWccDiam);
    //fprintf(F, "%f\t%f", FirstGrowth.MxEigVal, LastGrowth.MxEigVal);
    fprintf(F, "\n");
  }
  fclose(F);
}

void TKronPhasePlot::KroneckerPhase(const TStr& MtxId, const int& MxIter,
 const double& MnAlpha, const double& MxAlpha, const double& AlphaStep,
 const double& MnBeta, const double& MxBeta, const double& BetaStep,
 const TStr& FNmPref) {
  TKronPhasePlot PhasePlot;
  TExeTm KronTm;
  int AlphaCnt=0, BetaCnt=0;
  for (double Alpha = MnAlpha; (Alpha-1e-6) <= MxAlpha; Alpha += AlphaStep) {
    AlphaCnt++;  BetaCnt = 0;
    printf("\n\n****A:%g***********************************************************************", Alpha);
    for (double Beta = MnBeta; (Beta-1e-6) <= MxBeta; Beta += BetaStep) {
      printf("\n\n==A[%d]:%g====B[%d]:%g=====================================================\n", AlphaCnt, Alpha, BetaCnt, Beta);
      BetaCnt++;
      TGrowthStat GrowthStat;
      PNGraph Graph;
      // run Kronecker
      TFullRectMtx SeedMtx = TFullRectMtx::GetMtxFromNm(MtxId);
      SeedMtx.Epsilonize(Alpha, Beta);
      for (int iter = 1; iter < MxIter + 1; iter++) {
        printf("%2d] at %s\n", iter, TExeTm::GetCurTm().CStr());
        Graph = PNGraph();  KronTm.Tick();
        Graph = SeedMtx.GenRMatKronecker(iter, false, 0);
        GrowthStat.Add(Graph, TNodeTm(iter));
        if (KronTm.GetSecs() > 30 * 60) {
          printf("*** TIME LIMIT [%s]\n", KronTm.GetTmStr().CStr());  break; }
      }
      const TStr Desc = TStr::Fmt("%s. Alpha:%g. Beta:%g", MtxId.CStr(), Alpha, Beta);
      const TStr FNmPref1 = TStr::Fmt("%s.a%02d.b%02d", FNmPref.CStr(), AlphaCnt, BetaCnt);
      TGPlot::PlotDegDist(Graph, FNmPref1, Desc, false, true, true);
      TGPlot::PlotWccDist(Graph, FNmPref1, Desc);
      TGPlot::PlotSccDist(Graph, FNmPref1, Desc);
      GrowthStat.PlotAll(FNmPref1, Desc);
      GrowthStat.SaveTxt(FNmPref1, Desc);
      PhasePlot.PhaseV.Add(TPhasePoint(Alpha, Beta, GrowthStat));
    }
    {TFOut FOut(TStr::Fmt("phase.%s.bin", FNmPref.CStr()));
    PhasePlot.Save(FOut); }
  }
}
*/
/*void TKroneckerLL::SetRndThetas() {
  ProbMtx.Dump("TRUE parameters");
  TFltV RndV;
  double SumRnd = 0.0;
  for (int i = 0; i < ProbMtx.Len(); i++) {
    RndV.Add(0.1+TKronMtx::Rnd.GetUniDev());
    SumRnd += RndV.Last();
  }
  RndV.Sort(false);
  for (int i = 0; i < ProbMtx.Len(); i++) { ProbMtx.At(i) = RndV[i]; }
  ProbMtx.Dump("Random parameters");
  const double EdgePSum = pow(Graph->GetEdges(), 1.0/KronIters);
  bool Repeat = true;
  while (Repeat) {
    const double Scale = EdgePSum / SumRnd;
    Repeat=false;  SumRnd = 0.0;
    for (int i = 0; i < ProbMtx.Len(); i++) {
      ProbMtx.At(i) = ProbMtx.At(i)*Scale;
      if (ProbMtx.At(i) > 0.95) { ProbMtx.At(i)=0.95; Repeat=true; }
      SumRnd += ProbMtx.At(i);
    }
  }
  ProbMtx.Dump("INIT parameters");
  ProbMtx.GetLLMtx(LLMtx);
}*/

/*
void TKroneckerLL::TestLL() {
  TExeTm ExeTm;
  // approximation to empty graph log-likelihood
  */
  /*{ PNGraph Graph = TNGraph::New();
  for (uint i = 0; i < TMath::Pow2(4); i++) { Graph->AddNode(i); } //8k nodes
  PKronecker KronParam = TKronMtx::GetMtx("0.8 0.6; 0.7 0.3");
  TKroneckerLL KronLL(Graph, KronParam);
  printf("\nNodes: %d\n", KronLL.GetNodes());
  printf("Full Graph LL:               %f\n", KronLL.GetFullGraphLL());
  printf("Empty Graph Exact LL:        %f\n", KronLL.GetEmptyGraphLL());
  printf("Empty Approx x=log(1-x) LL:  %f\n", KronLL.GetApxEmptyGraphLL());
  printf("Empty Sample LL (100/node):  %f\n", KronLL.GetSampleEmptyGraphLL(Graph->GetNodes() * 100));
  KronLL.SetOrderPerm();
  printf("\nEdge    prob: %f, LL: %f\n", KronParam->GetEdgeProb(0,0,8), log(KronParam->GetEdgeProb(0,0,8)));
  printf("No Edge prob: %f, LL: %f\n", KronParam->GetNoEdgeProb(0,0,8), log(KronParam->GetNoEdgeProb(0,0,8)));
  printf("Empty Graph  LL:     %f\n", KronLL.CalcGraphLL());
  printf("Apx Empty Graph  LL: %f\n", KronLL.CalcApxGraphLL());
  Graph->AddEdge(0, 0);
  printf("add 1 edge.  LL:     %f\n", KronLL.CalcGraphLL());
  printf("Apx add 1 edge.  LL: %f\n", KronLL.CalcApxGraphLL()); }
  */

  // log-likelihood versus different Kronecker parameters
  /*{ PKronecker KronParam = TKronMtx::GetMtx("0.9 0.6; 0.6 0.2");
  PNGraph Graph = TKronMtx::GenKronecker(KronParam, 10, true, 10);
  TVec<PKronecker> ParamV;
  ParamV.Add(KronParam);
  ParamV.Add(TKronMtx::GetMtx("0.6 0.6; 0.6 0.5")); // sum = 2.3
  //ParamV.Add(TKronMtx::GetMtx("0.9 0.9; 0.4 0.1")); // sum = 2.3
  //ParamV.Add(TKronMtx::GetMtx("0.8 0.7; 0.6 0.2")); // sum = 2.3
  ParamV.Add(TKronMtx::GetMtx("0.9 0.9; 0.6 0.2")); // sum = 2.6
  for (int i = 0; i < ParamV.Len(); i++) {
    ParamV[i]->Dump();
    TKroneckerLL KronLL(Graph, ParamV[i]);
    for (int k = 0; k < 3; k++) {
      if (k==0) { KronLL.SetOrderPerm();  printf("Order permutation:\n"); }
      if (k==1) { KronLL.SetDegPerm();  printf("Degree permutation:\n"); }
      if (k==2) { KronLL.SetRndPerm();  printf("Random permutation:\n"); }
      const double LL = KronLL.CalcGraphLL(), aLL = KronLL.CalcApxGraphLL();
      printf("  Exact  Graph LL:  %f\n", LL);
      printf("  Approx Graph LL:  %f\n", aLL);
      printf("  error          :  %.12f\n", -fabs(LL-aLL)/LL);
    }
  } }
  */
  // exact vs. approximate log-likelihood
  /*{ PKronecker KronParam = TKronMtx::GetMtx("0.9 0.6; 0.6 0.2");
  PNGraph Graph = TKronMtx::GenFastKronecker(KronParam, 16, true, 0);
  TKroneckerLL KronLL(Graph, KronParam);
  TMom ExactLL, ApxLL;
  printf("Random permutation:\n");
  for (int i = 0; i < 100; i++) {
    KronLL.SetRndPerm();
    //ExactLL.Add(KronLL.CalcGraphLL());
    ApxLL.Add(KronLL.CalcApxGraphLL());
    //printf("  Exact  Graph LL:  %f\n", ExactLL.GetVal(ExactLL.GetVals()-1));
    printf("  Approx Graph LL:  %f\n", ApxLL.GetVal(ApxLL.GetVals()-1));
  }
  ExactLL.Def();  ApxLL.Def();
  //printf("EXACT:   %f  (%f)\n", ExactLL.GetMean(), ExactLL.GetSDev());
  printf("APPROX:  %f  (%f)\n", ApxLL.GetMean(), ApxLL.GetSDev());
  KronLL.SetOrderPerm();
  printf("Order permutation:\n");
  printf("  Exact  Graph LL:  %f\n", KronLL.CalcGraphLL());
  printf("  Approx Graph LL:  %f\n", KronLL.CalcApxGraphLL());
  }
  */

  // start from random permultation and sort it using bubble sort
  // compare the end result with ordered permutation
  //PKronecker KronParam = TKronMtx::GetMtx("0.9 0.6; 0.6 0.2");
  //PNGraph Graph = TKronMtx::GenFastKronecker(KronParam, 10, true, 1);
  /*PKronecker KronParam = TKronMtx::GetMtx("0.9 0.7; 0.9 0.5");
  PNGraph Graph = TKronMtx::GenFastKronecker(KronParam, 6, true, 2);
  TGAlg::SaveFullMtx(Graph, "kron32.tab");

  TKroneckerLL KronLL(Graph, KronParam);
  KronLL.SetOrderPerm();
  KronLL.LogLike = KronLL.CalcApxGraphLL();
  printf("  Approx Graph LL:   %f\n", KronLL.CalcApxGraphLL());
  printf("  swap 1-20:         %f\n", KronLL.SwapNodesLL(1, 20));
  printf("  swap 20-30:        %f\n", KronLL.SwapNodesLL(20, 30));
  printf("  swap 30-1:         %f\n", KronLL.SwapNodesLL(1, 30));
  printf("  swap 20-30:        %f\n", KronLL.SwapNodesLL(30, 20));
  IAssert(KronLL.GetPerm().IsSorted());
  KronLL.SetRndPerm();
  KronLL.LogLike = KronLL.CalcApxGraphLL();
  for (int i = 0; i < 1000000; i++) {
    const int nid1 = TInt::Rnd.GetUniDevInt(KronLL.Nodes);
    const int nid2 = TInt::Rnd.GetUniDevInt(KronLL.Nodes);
    printf("%3d]  swap LL:     %f\n", i, KronLL.SwapNodesLL(nid1, nid2));
  }
  printf("***   approx LL:   %f\n", KronLL.CalcApxGraphLL());
  printf("***    exact LL:   %f\n", KronLL.CalcGraphLL());
  */
  /*ExeTm.Tick();
  // bubble sort
  for (int i = 0; i < Graph->GetNodes()-1; i++) {
    for (int j = 1; j < Graph->GetNodes(); j++) {
      if (KronLL.GetPerm()[j-1] > KronLL.GetPerm()[j]) {
        const double oldLL = KronLL.GetLL();
        const double newLL = KronLL.SwapNodesLL(j-1, j);
        //const double trueLL = KronLL.CalcApxGraphLL();
        //printf("swap %3d - %3d:   old: %f   new: %f  true:%f\n",
        //    KronLL.GetPerm()[j-1], KronLL.GetPerm()[j], oldLL, newLL, trueLL);
      }
    }
  }
  //for (int i = 0; i < 100000; i++) {
  //  KronLL.SwapNodesLL(TInt::Rnd.GetUniDevInt(TMath::Pow2(16)), TInt::Rnd.GetUniDevInt(TMath::Pow2(16))); }
  printf("\nPermutation is %s\n", KronLL.GetPerm().IsSorted()? "SORTED" : "NOT SORTED");
  printf("  Swap   Graph LL:   %f\n", KronLL.GetLL());
  printf("  Approx Graph LL:   %f\n", KronLL.CalcApxGraphLL());
  KronLL.SetOrderPerm();
  printf("  Order  Graph LL:   %f\n\n", KronLL.CalcApxGraphLL());
  printf("Permutation is %s\n", KronLL.GetPerm().IsSorted()? "SORTED" : "NOT SORTED");
  printf("time: %f\n", ExeTm.GetSecs());
  */
  // evaluate the derivatives
  /*{ PNGraph Graph = TNGraph::New();
  TKronMtx KronParam = TKronMtx::GetMtx("0.8 0.4; 0.4 0.2");
  //for (uint i = 0; i < TMath::Pow2(4); i++) { Graph->AddNode(i); } //8k nodes
  Graph = TKronMtx::GenFastKronecker(KronParam, 8, true, 2); //TGAlg::SaveFullMtx(Graph, "kron16.txt");
  TKroneckerLL KronLL(Graph, KronParam);
  KronLL.SetOrderPerm();
  printf("\nNodes: %d\n", KronLL.GetNodes());
  printf("Full  Graph Exact LL:        %f\n", KronLL.GetFullGraphLL());
  printf("Empty Graph Exact LL:        %f\n", KronLL.GetEmptyGraphLL());
  printf("Empty Approx LL:             %f\n", KronLL.GetApxEmptyGraphLL());
  printf("Exact Graph LL:   %f\n", KronLL.CalcGraphLL());
  printf("Apx   Graph LL:   %f\n\n", KronLL.CalcApxGraphLL());
  // derivatives
  printf("Empty graph Exact DLL:           %f\n", KronLL.GetEmptyGraphDLL(0));
  printf("Empty graph Apx   DLL:           %f\n", KronLL.GetApxEmptyGraphDLL(0));
  printf("Theta0    edge(0,1) DLL:      %f\n", KronLL.LLMtx.GetEdgeDLL(0, 0, 1, 4));
  printf("Theta0 NO edge(0,1) DLL:      %f\n", KronLL.LLMtx.GetNoEdgeDLL(0, 0, 1, 4));
  printf("Theta0 NO edge(0,1) DLL:      %f\n", KronLL.LLMtx.GetApxNoEdgeDLL(0, 0, 1, 4));
  KronLL.CalcGraphDLL();  printf("Exact Theta0 DLL:");  KronLL.GetDLL(0);
  KronLL.CalcApxGraphDLL();  printf("Apx   Theta0 DLL:");  KronLL.GetDLL(0);
  KronLL.CalcFullApxGraphDLL();  printf("Apx   Theta0 DLL:");  KronLL.GetDLL(0);
  // swap
  */
  /*for (int i = 0; i < 100; i++) {
    const int A = TInt::Rnd.GetUniDevInt(KronLL.Nodes);
    KronLL.NodePerm.Swap(i, A);
    //KronLL.UpdateGraphDLL(i, A); printf("Fast  Theta0 DLL:");  KronLL.GetDLL(0);
    KronLL.CalcApxGraphDLL();  printf("Apx   Theta0 DLL:");  KronLL.GetDLL(0);
    //KronLL.CalcFullApxGraphDLL();  printf("Apx   Theta0 DLL:");  KronLL.GetDLL(0);
    //KronLL.CalcGraphDLL();  printf("Exact Theta0 DLL:");  KronLL.GetDLL(0);
    printf("\n");
  } */
  //}
//}

/*void TKroneckerLL::SampleGradient(const int& WarmUp, const int& NSamples, double& AvgLL, TFltV& AvgGradV, TFltV& SDevV, const bool& Plot) {
  printf("Samples: %s (warm-up: %s)\n", TInt::GetMegaStr(NSamples).CStr(), TInt::GetMegaStr(WarmUp).CStr());
  int NId1 = 0, NId2 = 0;
  TExeTm ExeTm;
  CalcApxGraphLL();
  for (int s = 0; s < WarmUp; s++) {
    SampleNextPerm(NId1, NId2); }
  printf("  warm-up:%s", ExeTm.GetTmStr());  ExeTm.Tick();
  CalcApxGraphLL();
  CalcApxGraphDLL();
  AvgLL = 0;
  TVec<TMom> DLLMomV(LLMtx.Len());
  for (int s = 0; s < NSamples; s++) {
    if (SampleNextPerm(NId1, NId2)) { // new permutation
      UpdateGraphDLL(NId1, NId2);
    }
    AvgLL += GetLL();
    for (int m = 0; m < LLMtx.Len(); m++) { DLLMomV[m].Add(GradV[m]); }
  }
  AvgLL = AvgLL / (NSamples*Nodes);
  // plot gradients over sampling time
  if (Plot) {
    TVec<TFltV> FltVV(LLMtx.Len()+1);
    for (int s = 0; s < DLLMomV[0].GetVals(); s += 1000) {
      for (int m = 0; m < LLMtx.Len(); m++) { FltVV[m].Add(DLLMomV[m].GetVal(s)); }
      FltVV.Last().Add(s); }
    const TStr FNm = TFile::GetUniqueFNm(TStr::Fmt("grad%dW%sS%s-#.png", KronIters, TInt::GetMegaStr(WarmUp).CStr(), TInt::GetMegaStr(NSamples).CStr()));
    TGnuPlot GP(FNm.GetFMid(), TStr::Fmt("Gradient vs. Sample Index. Nodes: %d, WarmUp: %s, Samples: %s, Avg LL: %f", Nodes,
      TInt::GetMegaStr(WarmUp).CStr(), TInt::GetMegaStr(NSamples).CStr(), AvgLL), true);
    for (int m = 0; m < LLMtx.Len(); m++) {
      GP.AddPlot(FltVV.Last(), FltVV[m], gpwLines, TStr::Fmt("Grad %d", m+1), "linewidth 5"); }
    GP.SetXYLabel("Sample Index (time)", "Log-likelihood gradient");
    GP.SavePng();
  }
  // average gradients
  printf("  sampling:%s\n", ExeTm.GetTmStr());
  printf("  AverageLL: %f\n", AvgLL);
  printf("Gradients:\n");
  AvgGradV.Gen(LLMtx.Len());
  SDevV.Gen(LLMtx.Len());
  for (int m = 0; m < LLMtx.Len(); m++) {
    DLLMomV[m].Def();
    AvgGradV[m] = DLLMomV[m].GetMean() / (Nodes*Nodes);
    SDevV[m] = DLLMomV[m].GetSDev() / (Nodes*Nodes);
    printf("  %d]  mean: %16f    sDev: %16f\n", m, AvgGradV[m], SDevV[m]);
  }
}

void TKronMaxLL::TFunc::FDeriv(const TFltV& Point, TFltV& GradV) {
  CallBack->GetDLL(Point, GradV);
  for (int i = 0; i < GradV.Len(); i++) { GradV[i] = -GradV[i]; }
}

double TKronMaxLL::TLogBarFunc::FVal(const TFltV& Point) {
  // log-likelihood
  const double LogLL = CallBack->GetLL(Point);
  // log-barrier
  const double MinBarrier = 0.05;
  const double MaxBarrier = 0.95;
  const double T1 = 1.0/T;
  double Barrier = 0.0;
  for (int i = 0; i < Point.Len(); i++) {
    if(Point[i].Val > MinBarrier && Point[i].Val < MaxBarrier) {
      Barrier += - T1 * (log(Point[i]-MinBarrier) + log(MaxBarrier-Point[i])); //log-barrier
    } else { Barrier = 1e5; }
  }
  IAssert(Barrier > 0.0);
  printf("barrrier: %f\n", Barrier);
  return -LogLL + Barrier; // minus LL since we want to maximize it
}

void TKronMaxLL::TLogBarFunc::FDeriv(const TFltV& Point, TFltV& DerivV) {
  // derivative of log-likelihood
  CallBack->GetDLL(Point, DerivV);
  // derivative of log barrier
  const double MinBarrier = 0.05;
  const double MaxBarrier = 0.95;
  const double T1 = 1.0/T;
  for (int i = 0; i < Point.Len(); i++) {
    DerivV[i] = - DerivV[i] + (- T1*(1.0/(Point[i]-MinBarrier) - 1.0/(MaxBarrier-Point[i])));
  }
}
*/
