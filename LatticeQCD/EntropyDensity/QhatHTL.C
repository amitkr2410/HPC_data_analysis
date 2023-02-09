//Fmunu in data files is (Umunu - Udaggermunu)/(2i)
//Here, We divide Fmunu by barecoupling constant g0=sqrt(6/Beta), and lattice spacing square a^2; See function Subtract()
//Thus, F^2_{munu} should be divided by g0^2 and a^4;
//We use perturbative two-loop function result to set the scale.
//Plot is for Correlator w.r.t Temperature.   
#include"iostream"
#include"Custom.h"

TGraphErrors * Scale(TGraphErrors *g1,double K);

void QhatHTL()
{
  char FileName[3000];
  double YMin=0.0, YMax=20.0;
  double XMin=200, XMax=500;
  TGraphErrors *G1 = new TGraphErrors("EntropyDensity_RawData_QCD.txt","%lf %lf");
  double qhat0=2.0*25, S0=96.0;  //s0=96 fm^-3 == 373MeV
  double T3Qhat0=2.0*0.2, T0=0.362; //T0=362 MeV matches s0=96 fm^-3 
  G1=Scale(G1, qhat0/S0);

  double Ca=3, PI=3.1415926, Zeta3=1.2, Alphas=0.25, Alphas2=0.2,Energy=10000, Energy2=100000;
  double PreFactor=0.0;
  PreFactor=Ca*42.0*Zeta3/PI;
  TF1 *f1 = new TF1("f1","[0]*pow([1],2.0)*log(5.7*[2]/(24*3.1415926*[1]*x))",XMin,XMax);
  f1->FixParameter(0,PreFactor);
  f1->FixParameter(1,Alphas);
  f1->FixParameter(2,Energy);

  TF1 *f2 = new TF1("f2","[0]*pow([1],2.0)*log(5.7*[2]/(24*3.1415926*[1]*x))",XMin,XMax);
  f2->FixParameter(0,PreFactor);
  f2->FixParameter(1,Alphas);
  f2->FixParameter(2,Energy2);

  TF1 *f3 = new TF1("f3","[0]*pow([1],-3.0)",XMin,XMax);
  f3->FixParameter(0,T3Qhat0);
  f3->FixParameter(1,T0);

  TCanvas *c1 = new TCanvas("c1","c1",580,500);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(5);
  gPad->SetMargin(0.14, 0.05, 0.13, 0.03);
  CustomGlobal();//1=Black, 2=Red, 4=light Blue, 9=DarkBlue, 6=Pink, 8=Green                                                                   
  Custom(G1,9);
  G1->Draw("AL");f1->Draw("L same"); f2->Draw("L same");f3->Draw("L same");
  G1->SetLineWidth(9);f1->SetLineWidth(5);f2->SetLineWidth(5);f3->SetLineWidth(5);
  f1->SetLineColor(2); f2->SetLineColor(2);f3->SetLineColor(6);
  //f2->SetLineStyle(2);
  G1->GetYaxis()->SetTitle("#it{#hat{q}/T^{3}}");G1->GetXaxis()->SetTitle("#it{T} (MeV)");
  G1->GetYaxis()->SetRangeUser(YMin,YMax);
  G1->GetXaxis()->SetLimits(XMin, XMax);
  
  TLegend *L1 = new TLegend(0.2,0.85,0.7,0.9);
  //L1->AddEntry(G1,"#hat{q}/T^{3}=#hat{q}_{0}*s/(s_{0}*T^{3})","lp");
  //L1->AddEntry(f1,"Ca*42*#zeta*#alpha_{s}^{2}*log(9E/(8PI*#alpha*T))","lp");
  //L1->AddEntry(f2,"Ca*42*#zeta*#alpha_{s}^{2}*log(9E/(8PI*#alpha*T))","lp");
  //L1->AddEntry(f3, "q_{0}=2GeV^{2}/fm, #hat{q}=q_{0}T^{3}/T^{3}_{0}, T_{0}=0.362MeV","l");
  //L1->AddEntry(G1, "q_{0}=2GeV^{2}/fm, #hat{q}=q_{0}s/s_{0}, s_{0}=96fm^{-3}","l");

  L1->AddEntry(f1, "See Legends on #it{R}_{AA}'s plot","");
  //L1->AddEntry(f1, "#alpha_{s}=0.25,Q_{0}=2GeV, #hat{q}=#it{c}_{1}#alpha_{s}T^{3}log(#it{c}_{2}E/(#alpha_{s}T))","");
  //L1->AddEntry(f2, "E=100GeV","");
  L1->SetFillStyle(0);L1->SetTextFont(22);L1->SetTextSize(0.055);  L1->Draw();
    TLegend *L2 = new TLegend(0.65,0.65,0.8,0.76);
  L2->AddEntry(f2, "E=100GeV","");
  L2->SetFillStyle(0);L2->SetTextFont(22);L2->SetTextSize(0.055);  L2->SetTextColor(2);L2->Draw();
  TLegend *L3 = new TLegend(0.65,0.33,0.8,0.5);
  L3->AddEntry(f1, "E=10GeV","");
  L3->SetFillStyle(0);L3->SetTextFont(22);L3->SetTextSize(0.055); L3->SetTextColor(2); L3->Draw();
  
  char PDFFileOutput[800000],COMMAND[800000];
  //sprintf(PDFFileOutput,"~/Dropbox/QhatLattice/Plots/EntropyDensityOverT3VsT.pdf");
  //c1->SaveAs("MyFile1.eps");
  //system("epstopdf MyFile1.eps");
  //sprintf(COMMAND,"mv MyFile1.pdf %s",PDFFileOutput);
  //cout<<string(COMMAND)<<endl;

  //system(COMMAND);
}



TGraphErrors * Scale(TGraphErrors *g1, double K)
{
  const int N1=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N3; double X3[N1], Y3[N1], Ex3[N1], Ey3[N1];
  for (int i=0; i<N1; i++)
    {
      X3[i] = X1[i];
      N3=N1;
      cout<<"X1[i]="<<X1[i]<<"\tX3[i]="<<X3[i]<<endl;
      cout<<"Y1[i]="<<Y1[i]<<endl;
      Y3[i] = Y1[i]*K;
      Ey3[i] =0.0;
      Ex3[i] =0.0;
    }
  TGraphErrors *g3 = new TGraphErrors(N3,X3,Y3,Ex3,Ey3);
  return g3; 
}
