#include"Custom.h"
#include"../Script/StatTool.h"
// define a function of calculating the ratio between two root graphs
TGraphErrors * ScaleGraphBeta(TGraphErrors *g, TGraph *g2, int PowerOfU0);
TGraphErrors *Divideg0Square(TGraphErrors *g, TGraph *g2, string QCDTypeS);
double RenormalizationConstant(double g0Square, string QCDTypeS);
TGraphErrors *Renormalize(TGraphErrors *g,TGraph *g2, string QCDTypeS);
double FindErrorBar(TGraphErrors *g, double a);
TGraphErrors *FitMyGraph(TGraphErrors *g, TGraphErrors *g2, TGraphErrors *g3);
void ComputeQhat(TGraphErrors *CG, string QCDTypeS); 
double StrongCouplingOneLoop(string QCDTypeS, double T);
TGraphErrors *EntropyFromRawDataPureSU3();
void FullQCDFF()
{
  CustomGlobal();
  char FileName[10000];
  double XMin=100, XMax=1020, YMin=0, YMax=4;
  //string QCDTypeS="PureGauge";
  string QCDTypeS="FullQCD";
  double *MinMax; double Diff=0; int LW=2;
  gStyle->SetLineWidth(1);
  //TCanvas *c = new TCanvas("c","c",500,500);
  TGraphErrors *g, *G, *g2, *G2, *g3, *G3;
  TGraphErrors *h, *H, *h2, *H2, *h3, *H3;
  TGraphErrors *R, *R2, *R3;
  int PowersOfU0Clover=4, PowersOfU0FourLinks=4; int Nt, Ns;
  
  Nt=4; Ns=16;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());  
  g = new TGraphErrors(FileName,"%lf %lf %lf ");
  //G = ScaleGraphBeta(g,T , PowersOfU0Clover);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  TGraphErrors *hh = new TGraphErrors(FileName,"%lf %lf %lf ");
  //h = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  //H = ScaleGraphBeta(h, T, PowersOfU0FourLinks);
  g = Divideg0Square(g,T,QCDTypeS);
  //G = Divideg0Square(G,T,QCDTypeS);
  //h = Divideg0Square(h,T,QCDTypeS);
  //H = Divideg0Square(H,T,QCDTypeS);

      
  Nt=6;	Ns=24;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T2 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  g2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  //G2 = ScaleGraphBeta(g2,T2 , PowersOfU0Clover);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  TGraphErrors *hh2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  //h2 = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  //H2 = ScaleGraphBeta(h2, T2, PowersOfU0FourLinks);

  g2 = Divideg0Square(g2,T2,QCDTypeS);
  //G2 = Divideg0Square(G2,T2,QCDTypeS);
  //h2 = Divideg0Square(h2,T2,QCDTypeS);
  //H2 = Divideg0Square(H2,T2,QCDTypeS);

  

  Nt=8;	Ns=32;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T3 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  g3 = new TGraphErrors(FileName,"%lf %lf %lf ");
  //G3 = ScaleGraphBeta(g3,T3 , PowersOfU0Clover);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  TGraphErrors *hh3 = new TGraphErrors(FileName,"%lf %lf %lf ");
  //h3 = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  //H3 = ScaleGraphBeta(h3, T3, PowersOfU0FourLinks);

  g3 = Divideg0Square(g3,T3,QCDTypeS);
  //G3 = Divideg0Square(G3,T3,QCDTypeS);
  // h3 = Divideg0Square(h3,T3,QCDTypeS);
  //H3 = Divideg0Square(H3,T3,QCDTypeS);

    
  
  TCanvas *c0 = new TCanvas("c0","c0",450,400);
  TGraphErrors *S, *S2; S = EntropyFromRawDataPureSU3();
  S2 = ScaleTGraphErrors(S,0.1);
  S = ScaleTGraphErrors(S,0.15);
  
  Custom(g,1);Custom(g2,2); Custom(g3,4);
  // Custom(G,1);Custom(G2,2); Custom(G3,4);
  //Custom(h,1);Custom(h2,2);Custom(h3,4);
  //Custom(H,1);Custom(H2,2);Custom(H3,4);
  
  g->Draw("Ap"); g2->Draw("p same");  g3->Draw("p same");  
  //G->Draw("Lp");  G2->Draw("Lp same");	G3->Draw("Lp same");
  //h->Draw("p");  h2->Draw("p"); h3->Draw("p");
  //H->Draw("Lp");  H2->Draw("Lp"); H3->Draw("Lp");  
    
  Custom(S,2); S->Draw("L");
  Custom(S2,2); S2->Draw("L");
  
  g->SetLineWidth(LW);g2->SetLineWidth(LW); g3->SetLineWidth(LW);

  S->SetLineWidth(4);
  S2->SetLineWidth(4); S2->SetLineStyle(7);
  
  g->GetYaxis()->SetTitle("<#it{O}>");g->GetXaxis()->SetTitle("#it{T} (MeV)");
  g->GetYaxis()->SetRangeUser(YMin,YMax);
  g->GetXaxis()->SetLimits(XMin, XMax);
  g->GetXaxis()->SetNdivisions(306);
  TPaveText *Pt1 = new TPaveText(XMin+60,YMax-1.4,XMin+150,YMax-0.43,"NB");
  Pt1->AddText("#it{n}_{#tau}=4");Pt1->AddText("#it{n}_{#tau}=6");Pt1->AddText("#it{n}_{#tau}=8");
  Pt1->SetTextAlign(21);Pt1->SetFillColor(0);  Pt1->SetTextFont(22);Pt1->SetTextSize(0.05);  Pt1->Draw();

  TPaveText *Pt2 = new TPaveText(XMin+170,YMax-0.5,XMin+330,YMax-0.1,"NB");
  Pt2->AddText("#it{FF}/#it{T}^{4}");
  Pt2->SetTextAlign(21);Pt2->SetFillColor(0);  Pt2->SetTextFont(22);Pt2->SetTextSize(0.05);  Pt2->Draw();

  //TPaveText *Pt3 = new TPaveText(XMin+370,YMax-0.8,XMin+450,YMax-0.1,"NB");  
  //Pt3->AddText("#it{Z}^{(3)}_{T}#it{FF}/#it{T}^{4}");
  //Pt3->SetTextAlign(21);Pt3->SetFillColor(0);  Pt3->SetTextFont(22);Pt3->SetTextSize(0.05);  Pt3->Draw();

  TLEGEND3(0.32,0.7,0.42,0.92, g,"","lp",  g2,"","lp", g3,"","lp", 0.04);
  //TLEGEND3(0.32,0.7,0.4,0.92,  G,"FF/(T^{4}#it{u}^{4}_{0})","lp",  G2,"FF/(T^{4}#it{u}^{4}_{0})","lp", G3,"FF/(T^{4}#it{u}^{4}_{0})","lp", 0.04);
  //TLEGEND3(0.475,0.7,0.65,0.92, R,"","lp", R2,"","lp",R3,"","lp",   0.04);
  //TLEGEND1(0.26,0.1,0.83,0.32, S, "0.5#it{s}/#it{T}^{3} [LNP 583 209 (2002)]","l", 0.05);
    TLEGEND1(0.26,0.22,0.83,0.3, S, "0.15#it{s}/#it{T}^{3} [JHEP 1011 077 (2010)]","l", 0.05);
    TLEGEND1(0.26,0.15,0.83,0.25, S2, "0.10#it{s}/#it{T}^{3} [JHEP 1011 077 (2010)]","l", 0.05);
  //TLEGEND3(0.475,0.7,0.65,0.92, h,"(1x1+1x2)#it{g}^{2}_{0}FF/T^{4};","lp", h2,"(1x1+1x2)#it{g}^{2}_{0}FF/T^{4};","lp",h3,"(1x1+1x2)#it{g}^{2}_{0}FF/T^{4};","lp",   0.04);
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
if(N!=N2){ cout<<"Error: N="<<N<<", N2="<<N2<<endl;}
  for(int i=0;i<N;i++)
    {
      x[i]=X[i];
      ex[i]=EX[i];//cout<<"X[i] = "<<endl;
      y[i]  = Y[i]/pow(Y2[i],PowerOfU0);
      ey[i] = EY[i]/pow(Y2[i],PowerOfU0);
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

TGraphErrors *Divideg0Square(TGraphErrors *g, TGraph *g2, string QCDTypeS)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();
  double x[N], y[N], ex[N], ey[N]; double g0Square=0;
  if(N!=N2){ cout<<"Error: N="<<N<<", N2="<<N2<<endl;}
  for(int i=0;i<N;i++)
    {
      x[i]=X[i];
      ex[i]=EX[i];
      if(QCDTypeS =="FullQCD")   { g0Square = 10.0/X2[i]; }
      if(QCDTypeS =="PureGauge") { g0Square =  6.0/X2[i]; }
      y[i]  = Y[i]/g0Square ;
      ey[i] = EY[i]/g0Square ;
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}



double StrongCouplingOneLoop(string QCDTypeS, double T)
{
  int Nf=0; 
  double QSquare = pow(3.1415926*2.0*T, 2.0);// MeV^2
  double LambdaSquare = pow(340, 2.0); // in MS-bar scheme, MeV^2 https://arxiv.org/pdf/1604.08082.pdf
  if(QCDTypeS == "FullQCD" ){Nf=3;}
  else {Nf=0; LambdaSquare= pow(243,2);} //https://inspirehep.net/literature/554365
  double beta0= 11 - (2.0*Nf/3.0); 
  double alphas = 4*3.1415926/(beta0*log(QSquare/LambdaSquare));cout<<"alphas="<<alphas<<"\t";
  return alphas;
}

TGraphErrors * EntropyFromRawDataPureSU3()
{
  TGraph *g1 = new TGraph("../EntropyDensity/EntropyDensity_RawData_QCD.txt","%lf %lf ");
  //TGraphErrors *g1 = new TGraphErrors("../EntropyDensity/EntropyDensity_RawData_SU3_MovingFrame.txt","%lf %lf %lf");
  const int N1=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(); //double *EY1=g1->GetEY();
  int N3; double X3[N1], Y3[N1], Ex3[N1], Ey3[N1];
  for (int i=0; i<N1; i++)
    {
      X3[i] = X1[i];
      N3=N1;
      cout<<"X1[i]="<<X1[i]<<"\tX3[i]="<<X3[i]<<endl;
      cout<<"Y1[i]="<<Y1[i]<<endl;
      Y3[i] = Y1[i];      Ey3[i] =0.0;      Ex3[i] =0.0;
    }
  TGraphErrors *g3 = new TGraphErrors(N3,X3,Y3,Ex3,Ey3);
  return g3;
}
