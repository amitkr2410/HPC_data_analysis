#include"iostream"
#include"Custom.h"

void PlotRe();
void PlotIm();
void PlotwrtBeta()
{
  PlotRe();
  //PlotIm();
  gPad->SetMargin(0.15,0.05,0.2,0.02);
}

void PlotRe()
{
  double Ymin = -0.2;
  double Ymax = 2.0;
  int Nt=6, Ns=4*Nt;
  //int Nt=4, Ns=Nt;
  char FileName[500];
  sprintf(FileName,"UnfinishedResult2/FFCorrelatorLO_NLO_Nt%i_Ns%i_wrtBeta.txt",Nt,Ns);
  TGraphErrors *Re1 = new TGraphErrors(FileName,"%lf %lf %lf"); 
  TGraphErrors *Re2 = new TGraphErrors(FileName,"%lf %*lf %*lf %*lf %*lf %lf %lf");
  TGraphErrors *Re3 = new TGraphErrors(FileName,"%lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %lf %lf");
  TGraphErrors *Re4 = new TGraphErrors(FileName,"%lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf  %lf %lf");
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(5);

  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  CustomGlobal();
  Custom(Re1,2);Custom(Re2,9);Custom(Re3,1);Custom(Re4,3);

  Re1->Draw("A*");Re2->Draw("*");Re3->Draw("*");Re4->Draw("*");
  Re1->GetYaxis()->SetTitle("<O>");Re1->GetXaxis()->SetTitle("#beta");
  Re1->GetYaxis()->SetRangeUser(Ymin, Ymax);

  TLegend *L = new TLegend(0.5,0.7,0.98,0.94);
  
  L->AddEntry(Re1,"F^{3i}F^{3i}-F^{4i}F^{4i}","el");
  L->AddEntry(Re2,"F^{3i}F^{4i}+F^{4i}F^{3i}","el");
  L->AddEntry(Re3,"F^{3i}DzF^{3i}-F^{4i}DzF^{4i}","el");
  L->AddEntry(Re4,"F^{3i}DzF^{4i}+F^{4i}DzF^{3i}","el");
  L->SetFillStyle(0);L->SetTextFont(22);L->SetTextSize(0.04);
  L->Draw();
}

