<< /home/amit/.Mathematica/Applications/FeynCalc/FeynCalc.m

L1 = p - l;
L3 = p - l;
L2 = k + p - l;
PS = DiracMatrix[\[Rho]];
LThreeS = DiracSlash[L3];
LTwoS = DiracSlash[L2];
LOneS = DiracSlash[L1];
GBeta = DiracMatrix[\[Beta]];
GNu = DiracMatrix[\[Nu]];
GMu = DiracMatrix[\[Mu]];
GAlpha = DiracMatrix[\[Alpha]];
ScalarProduct[n, n] = 0;
ScalarProduct[l, l] = ScalarProduct[p, p] = ScalarProduct[L2, L2] 0;
PL = PolarizationSum[\[Alpha], \[Beta], l, n];
ScalarProduct[l, n] = l2p/(2*y*x*Pplus);
ScalarProduct[l, p] = l2p/(2*y);
ScalarProduct[n, p] = 0;
ScalarProduct[k, n] = k2p/(2*(1 - y)*x*Pplus);
ScalarProduct[k, l] = k2p*y/(2*(1 - y));
ScalarProduct[k, p] = k2p/(2*(1 - y));
//Contract[ Tr[   PS. GBeta . LThreeS . GNu . LTwoS . GMu . LOneS .   GAlpha] PolarizationSum[ \[Alpha], \[Beta], l, n]];
HH = GluonVertex[ {k,\[Mu] },  {-l,\[Alpha]},   {p,\[Nu]}] // Explicit

Simplify[Contract[ HH  Tr[GMu. GNu. GAlpha. LOneS] ]]
