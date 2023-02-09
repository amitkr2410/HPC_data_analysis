<< /home/amit/.Mathematica/Applications/FeynCalc/FeynCalc.m

L1 = p - l;
L3 = p - l;
L2 = k + p - l;
PS = DiracSlash[p];
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
ScalarProduct[l, n] = l2perp/(2*y*x*Pplus);
ScalarProduct[l, p] = l2perp/(2*y);
ScalarProduct[n, p] = 0;
ScalarProduct[k, k] = -k2perp;
ScalarProduct[k, n] = k2perp/(2*(1 - y)*x*Pplus);
ScalarProduct[k, l] = k2perp*y/(2*(1 - y));
ScalarProduct[k, p] = k2perp/(2*(1 - y));
Contract[Tr[ PS. GBeta . LThreeS . GNu . LTwoS . GMu . LOneS . GAlpha] PolarizationSum[ \[Alpha], \[Beta], l, n] ]

