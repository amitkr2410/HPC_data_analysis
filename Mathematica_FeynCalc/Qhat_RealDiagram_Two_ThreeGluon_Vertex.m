<< /home/amit/.Mathematica/Applications/FeynCalc/FeynCalc.m

LThreeS = DiracSlash[p-p2];
PTwoS = DiracSlash[p2];
Q2S = DiracSlash[q-k];
QS = DiracSlash[q];
GFakeS = DiracSlash[p];

GRho = DiracMatrix[\[Rho]];
GSigma = DiracMatrix[ \[Sigma]]
GNu = DiracMatrix[\[Nu]];
GMu = DiracMatrix[\[Mu]];
GAlpha = DiracMatrix[\[Alpha]];
GBeta = DiracMatrix[\[Beta]];
GPhi = DiracMatrix[\[Phi]];


ScalarProduct[n, n] = 0;
ScalarProduct[l, l] = ScalarProduct[p2, p2] = ScalarProduct[q-k,q-k]=0;
ScalarProduct[p,p]  = 0;
ScalarProduct[k,k]  = -k2perp;
ScalarProduct[q,q] = 0;

ScalarProduct[q, k] = -k2perp/2;
ScalarProduct[q, n] = qminus;
ScalarProduct[q, p] = qminus*x*Pplus;
ScalarProduct[q, p2] = qminus*(1-y)*x*Pplus;
ScalarProduct[q,l] = qminus*y*x*Pplus;

ScalarProduct[k,n] = k2perp/(2*y*x*Pplus);
ScalarProduct[k, p] = k2perp/(2*y);
ScalarProduct[k, p2] = (1-y)*k2perp/(2*y);
ScalarProduct[k, l] = -k2perp/2;

ScalarProduct[n, p] = 0;
ScalarProduct[n, p2] = p2perp/(2*(1-y)*x*Pplus);
ScalarProduct[n, l] = k2perp/(2*y*x*Pplus);

ScalarProduct[p, p2] = p2perp/(2*(1-y));
ScalarProduct[p, l] = k2perp/(2*y);
ScalarProduct[p2,l] = (1-y)*k2perp/(2*y);

PzSum1 = PolarizationSum[\[Mu], \[Sigma], k, n];
PzSum2 = PolarizationSum[\[Nu], \[Rho], k, n];

GLVertex1 = GluonVertex[ {k,\[Mu] },  {-l,\[Epsilon]},   {p-p2,\[Chi]}  ] // Explicit 

PzSum3 = PolarizationSum[\[Alpha], \[Chi], p-p2, n];
PzSum4 = PolarizationSum[\[Beta], \[Kappa], p-p2, n];
PzSum5 = PolarizationSum[\[Epsilon], \[Phi], l, n];

GLVertex2 = GluonVertex[ {-k,\[Nu]},   {-p+p2,\[Kappa]}, {l,\[Phi] } ] // Explicit
     

Simplify[Contract[ Tr[QS. GRho. Q2S. GSigma]. PzSum1. PzSum2 . GLVertex1 . GLVertex2 . Tr[GFakeS . GBeta . PTwoS . GAlpha]. PzSum3. PzSum4. PzSum5 ]/(x*Pplus)]

