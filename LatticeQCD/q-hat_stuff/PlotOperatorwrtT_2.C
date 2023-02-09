//Plot of g^{2}*Tr[FF]/T^{4} vs T in MeV
//for Clover-definition, single-plaquette definition, Square+Rectangle definition, Square definition
//Options to pick Pure SU(3) Vs Full QCD
#include"Custom.h"
#include"../Script/StatTool.h"
// define a function of calculating the ratio between two root graphs
TGraphErrors * ScaleGraphBeta(TGraphErrors *g, TGraph *g2, int PowerOfU0);
TGraphErrors *ConvertToT(TGraphErrors *g, TGraphErrors *g2, int Nt);
void PlotOperatorwrtT_2()
{
  CustomGlobal();
  char FileName[10000];
  double XMin=150, XMax=1000, YMin=0, YMax=10;
  //string QCDTypeS="PureGauge";
  string QCDTypeS="FullQCD";
  double *MinMax; double Diff=0; int LW=2;
  gStyle->SetLineWidth(1);
  TCanvas *c0 = new TCanvas("c0","c0",1600,900);
  TGraphErrors *g, *G, *g2, *G2, *g3, *G3,  *g4, *G4;
  TGraphErrors *h, *H, *h2, *H2, *h3, *H3, *h4, *H4;
  int PowersOfU0Clover=4, PowersOfU0FourLinks=4; int Nt, Ns;
  
  Nt=4; Ns=16;  
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());  
  g = new TGraphErrors(FileName,"%lf %lf %lf ");
  sprintf(FileName,"%sSinglePlaquetteTracelessVacuumSubtractedOperatorNormalizedVsBeta_RealPartNt%d.txt", QCDTypeS.c_str(), Nt);
  g2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  g2 = ConvertToT(g,g2,Nt);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  g3 = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  g4 = new TGraphErrors(FileName,"%lf %*lf %*lf  %*lf %*lf  %lf %lf");

  Nt=6; Ns=24;
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  G = new TGraphErrors(FileName,"%lf %lf %lf ");
  sprintf(FileName,"%sSinglePlaquetteTracelessVacuumSubtractedOperatorNormalizedVsBeta_RealPartNt%d.txt", QCDTypeS.c_str(), Nt);
  G2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  G2 = ConvertToT(G,G2,Nt);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  G3 = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  G4 = new TGraphErrors(FileName,"%lf %*lf %*lf  %*lf %*lf  %lf %lf");


  Nt=8; Ns=32;
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  h = new TGraphErrors(FileName,"%lf %lf %lf ");
  //sprintf(FileName,"%sSinglePlaquetteTracelessVacuumSubtractedOperatorNormalizedVsBeta_RealPartNt%d.txt", QCDTypeS.c_str(), Nt);
  //h2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  //h2 = ConvertToT(h,h2,Nt);
  
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  h3 = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  h4 = new TGraphErrors(FileName,"%lf %*lf %*lf  %*lf %*lf  %lf %lf");
  
  
  Custom(g,1);Custom(g2,1); Custom(g3,1); Custom(g4,1);
  Custom(G,2);Custom(G2,2); Custom(G3,2); Custom(G4,2);
  Custom(h,4);//Custom(h2,4);
  Custom(h3,4); Custom(h4,4);

  g->Draw("Ap"); g2->Draw("Lp same");  g3->Draw("p same"); g4->Draw("Lp same");
  G->Draw("p");  G2->Draw("Lp same");	G3->Draw("p same"); G4->Draw("Lp same");
  h->Draw("p");  //h2->Draw("Lp");
  h3->Draw("p"); h4->Draw("Lp");

  g2->SetLineWidth(4);G2->SetLineWidth(4); //h2->SetLineWidth(4);
  g2->SetLineStyle(5);G2->SetLineStyle(5); //h2->SetLineStyle(5);  
  g2->SetMarkerStyle(24);G2->SetMarkerStyle(24); //h2->SetMarkerStyle(24);

  g3->SetMarkerStyle(21); G3->SetMarkerStyle(21);  h3->SetMarkerStyle(21);

  g4->SetLineWidth(LW); G4->SetLineWidth(LW); h4->SetLineWidth(LW);
  g4->SetMarkerStyle(25); G4->SetMarkerStyle(25);  h4->SetMarkerStyle(25);
  
  g->GetYaxis()->SetTitle("<#it{O}>");g->GetXaxis()->SetTitle("#it{T} (MeV)");
  MinMax=FindMinMax4(g, G2, G3, G4); Diff = (MinMax[1]-MinMax[0])/6.0;
  g->GetYaxis()->SetRangeUser(YMin,YMax);
  g->GetXaxis()->SetRangeUser(XMin, XMax);

  TLegend *leg1 = new TLegend(0.44,0.71,0.54,0.9);
  leg1->SetTextFont(22);leg1->SetTextSize(0.035);leg1->SetFillStyle(0);leg1->SetTextColor(0);
  leg1->AddEntry(g,"lp");  leg1->AddEntry(g2,"lp"); leg1->AddEntry(g3,"lp");   leg1->AddEntry(g4,"lp");
  leg1->Draw();

  TLegend *leg11 = new TLegend(0.54,0.71,0.64,0.9);
  leg11->SetTextFont(22);leg11->SetTextSize(0.035);leg11->SetFillStyle(0);leg11->SetTextColor(0); 
  leg11->AddEntry(G,"","lp");  leg11->AddEntry(G2,"","lp");  leg11->AddEntry(G3,"","lp"); leg11->AddEntry(G4,"","lp");
  leg11->Draw();
  
  TLegend *leg111 = new TLegend(0.68,0.71,0.9,0.9);
  leg111->SetTextFont(22);leg111->SetTextSize(0.035);leg111->SetFillStyle(0);leg111->SetTextColor(0);
  leg111->AddEntry(h,"","lp");  //leg111->AddEntry(h2,"","lp");
  leg111->AddEntry(h3,"","lp");
  leg111->AddEntry(h4,"","lp");
  leg111->Draw();

  TPaveText *Pt1 = new TPaveText(XMin+60,YMax/1.45,XMin+150,YMax/1.06,"NB");
  Pt1->SetTextFont(22);Pt1->SetTextSize(0.05);Pt1->SetFillColor(0); Pt1->SetTextAlign(12);
  Pt1->AddText("Clover:"); Pt1->AddText("Single Plaquette:"); Pt1->AddText("Square+Rectangle:"); Pt1->AddText("Square:");
  Pt1->Draw();

  TPaveText *Pt111 = new TPaveText(XMin+300,YMax/1.11,XMax-100,YMax/1.06,"NB");
  Pt111->SetTextFont(22);Pt111->SetTextSize(0.045);Pt111->SetFillColor(0);Pt111->SetTextAlign(11);
  Pt111->AddText("#it{n_{#tau}}=4             #it{n_{#tau}}=6              #it{n_{#tau}}=8");
  Pt111->Draw();
  
  //TLEGEND4(0.12,0.7,0.3,0.92, g,"Clover: #it{n}_{#tau}=4","lp", g2,"Single Plaquette:","lp",g3,"Square+Rectangle:","lp", g4,"Square:","lp", 0.04);
  //TLEGEND4(0.32,0.7,0.4,0.92,G,"#it{n}_{#tau}=6","lp", G2,"-","lp",G3,"-","lp", G4,"-","lp", 0.0);
  //TLEGEND4(0.475,0.7,0.65,0.92, h,"#it{n}_{#tau}=8","lp", h2,"-","lp", h3,"-","lp", h4,"-","lp", 0.04);
  //TLEGEND3(0.7,0.7,0.95,0.92, H,"(1x1+1x2)#it{g}^{2}_{0}FF/(T^{4}#it{u}^{4}_{0})","lp", H2, "(1x1+1x2)#it{g}^{2}_{0}FF/(T^{4}#it{u}^{4}_{0})","lp", H3,"(1x1+1x2)#it{g}^{2}_{0}FF/(T^{4}#it{u}^{4}_{0})","lp", 0.04);
  g->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetMargin(0.08,0.02,0.12,0.02);
  
  sprintf(FileName,"Graph%s_Plot.pdf", QCDTypeS.c_str());
  //c0->SaveAs("MyFile.eps");
  //system("epstopdf MyFile.eps");char COMMAND[100000];
  //sprintf(COMMAND,"mv MyFile.pdf %s",FileName);
  //cout<<string(COMMAND)<<endl;
  //system(COMMAND);

}

TGraphErrors * ScaleGraphBeta(TGraphErrors *g, TGraph *g2, int PowerOfU0)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();
  double x[N], y[N], ex[N], ey[N];
  cout<<"N="<<N<<", N2="<<N2<<endl;
  for(int i=0;i<N;i++)
    {
      x[i]=X[i];
      ex[i]=EX[i];
      y[i]  = Y[i]/pow(Y2[i],PowerOfU0);
      ey[i] = EY[i]/pow(Y2[i],PowerOfU0);
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

TGraphErrors *ConvertToT(TGraphErrors *g, TGraphErrors *g2, int Nt)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY(); double *EY2=g2->GetEY();
  double x[N], y[N], ex[N], ey[N];
  cout<<"N="<<N<<", N2="<<N2<<endl;
  for(int i=0;i<N;i++)
    {
      x[i]=X[i];
      ex[i]=EX[i];
      y[i]  = Y2[i]*pow(Nt,4.0);
      ey[i] = EY2[i]*pow(Nt,4.0);
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}  

