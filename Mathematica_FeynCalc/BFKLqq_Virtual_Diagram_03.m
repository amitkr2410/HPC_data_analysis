<< /home/amit/.Mathematica/Applications/FeynCalc/FeynCalc.m

LOneS = DiracSlash[p1+q2];
LTwoS = DiracSlash[p2-q];
LThreeS = DiracSlash[p2-q-q2];
LFourS = DiracSlash[p2-q2];
Q2S = DiracSlash[q2];
QS = DiracSlash[q];
POneS = DiracSlash[p1];
PTwoS = DiracSlash[p2];

GRho = DiracMatrix[\[Rho]];
GSigma = DiracMatrix[ \[Sigma]]
GNu = DiracMatrix[\[Nu]];
GMu = DiracMatrix[\[Mu]];
GAlpha = DiracMatrix[\[Alpha]];
GBeta = DiracMatrix[\[Beta]];
GPhi = DiracMatrix[\[Phi]];
GEpsilon = DiracMatrix[\[Epsilon]];
GChi = DiracMatrix[\[Chi]];

ScalarProduct[n, n] = 0;
ScalarProduct[p1, p1] = ScalarProduct[p2, p2] = ScalarProduct[p1+q2,p1+q2]=0;
ScalarProduct[p2-q2,p2-q2]  = 0;
ScalarProduct[q2,q2]  = -qtwo2perp;
ScalarProduct[q,q] = -q2perp;

ScalarProduct[p1, p2] = Pplus*Pminus; 
ScalarProduct[p1, q] = Pplus*Qminus; 
ScalarProduct[p1, q2] = qtwo2perp/2.0;
ScalarProduct[p1, n] = Pplus;


ScalarProduct[p2,q] = Pminus*Qplus;
ScalarProduct[p2, q2] = - qtwo2perp/2.0; 
ScalarProduct[p2, n] = 0;


ScalarProduct[q, q2] = -QperpDotQtwoPerp; 
ScalarProduct[q, n]  = Qplus; 

ScalarProduct[q2,n] = - qtwo2perp/(2.0*Pminus);

PzSum1 = PolarizationSum[\[Mu], \[Nu], q2, n];
PzSum2 = PolarizationSum[\[Beta], \[Chi], q2, n];
PzSum3 = PolarizationSum[\[Alpha], \[Epsilon], q, n];



Simplify[Contract[ Tr[POneS. GMu. LOneS. GBeta]. PzSum1. PzSum2 . Tr[PTwoS. GNu. LFourS . GAlpha . LThreeS . GChi . LTwoS . GEpsilon ]. PzSum3 ]]

