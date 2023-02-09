#include"Custom.h"
#include"../Script/StatTool.h"
// define a function of calculating the ratio between two root graphs
TGraphErrors *Divideg0Square(TGraphErrors *g, TGraph *g2, string QCDTypeS);
TGraphErrors *RenormalizeFFoT4(TGraphErrors *g, TGraph *g2, string QCDTypeS);
double RenormalizationConstant(double g0Square, string QCDTypeS);
TGraphErrors *DiscretizeHISQEntropy();
TGraphErrors *ConvertEQCD();
  
void FullQCD_FF_s()
{
  CustomGlobal();
  char FileName[10000];
  double XMin=120, XMax=1000, YMin=0, YMax=5;
  //string QCDTypeS="PureGauge";
  string QCDTypeS="FullQCD";
  double *MinMax; double Diff=0; int LW=2;
  gStyle->SetLineWidth(2);
  //TCanvas *c = new TCanvas("c","c",500,500);
  TGraphErrors *g, *G, *g2, *G2, *g3, *G3;
  int Nt, Ns;    

    Nt=4; Ns=16; QCDTypeS="FullQCD";
    sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
    TGraph *T = new  TGraph(FileName,"%lf %*lf %*lf %lf");
    sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
    g = new TGraphErrors(FileName,"%lf %lf %lf");
    G = new TGraphErrors(FileName,"%lf %*lf %*lf %lf %lf ");
    g = Divideg0Square(g,T,QCDTypeS);
    G = Divideg0Square(G,T,QCDTypeS);
    
    Nt=6;    Ns=24; QCDTypeS="FullQCD";
    sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
    TGraph *T2 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
    sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
    g2 = new TGraphErrors(FileName,"%lf %lf %lf");
    G2 = new TGraphErrors(FileName,"%lf %*lf %*lf %lf %lf ");
    g2 = Divideg0Square(g2,T2,QCDTypeS);
    G2 = Divideg0Square(G2,T2,QCDTypeS);

    Nt=8;    Ns=32; QCDTypeS="FullQCD";
    sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
    TGraph *T3 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
    sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
    g3 = new TGraphErrors(FileName,"%lf %lf %lf");
    G3 = new TGraphErrors(FileName,"%lf %*lf %*lf %lf %lf ");
    g3 = Divideg0Square(g3,T3,QCDTypeS);
    G3 = Divideg0Square(G3,T3,QCDTypeS);
    
    TGraph *R, *R2, *R3;

    TCanvas *c = new TCanvas("c","c",500,400);
    //Renormalize FF/T^4
    g  = RenormalizeFFoT4(g, T,QCDTypeS);
    g2 = RenormalizeFFoT4(g2,T2,QCDTypeS);
    g3 = RenormalizeFFoT4(g3,T3,QCDTypeS);
    
    //HISQ
    TGraphErrors *H;
    //= new TGraphErrors("../EntropyDensity/EntropyDensity_HISQ_nf3.txt","%lf %*lf %*lf %*lf %*lf %*lf %*lf %lf %lf");
    H = DiscretizeHISQEntropy();
    H = ScaleGraph(H, 0.33/2.0); //FullQCD/PureSu3
    
    //EQCD
    TGraphErrors *EQCD;
    EQCD = ConvertEQCD();
    EQCD = ScaleGraph(EQCD, 1.0/2.0);
    
    
    TGraph *gEQCD;
    gEQCD = MakeBand(EQCD);
    gEQCD->SetFillStyle(1001);gEQCD->SetLineColorAlpha(4,0.7);gEQCD->SetFillColorAlpha(4,0.7);

    TGraph *DummyH = new TGraph(FileName,"%lf %lf");
    TGraph *DummyEQCD = new TGraph(FileName,"%lf %lf");
    TCanvas *c0 = new TCanvas("c0","c0",650,450);
    Custom(G,2); Custom(G2,4); Custom(G3,1);
    Custom(g,2); Custom(g2,4); Custom(g3,1);
    Custom(H,2); Custom(EQCD,4);
    G->Draw("Ap"); G2->Draw("p same"); G3->Draw("p same");

    gEQCD->Draw("f");EQCD->SetLineWidth(14);
    gStyle->SetLineStyleString(11,"100 100");
    H->Draw("L "); H->SetLineStyle(11); H->SetLineWidth(14);
    H->SetLineColorAlpha(2,0.6);
    DummyEQCD->SetLineColorAlpha(4,0.7); DummyEQCD->SetLineWidth(4);DummyEQCD->SetLineStyle(1);
    DummyH->SetLineColorAlpha(2,0.6); DummyH->SetLineWidth(4); DummyH->SetLineStyle(11);
    G->GetYaxis()->SetTitle("<#it{O}>");
    G->GetXaxis()->SetTitle("#it{T} (MeV)");
    G->GetYaxis()->SetRangeUser(YMin,YMax);
    G->GetXaxis()->SetRangeUser(XMin,XMax);
    G->GetXaxis()->SetNdivisions(307);
    G->GetYaxis()->SetTitleOffset(0.7);
    G->GetXaxis()->SetTitleOffset(1.0);
    gPad->SetMargin(0.1,0.02,0.14,0.02);    
    G->SetMarkerSize(2);G2->SetMarkerSize(2);G3->SetMarkerSize(2);
    G->SetMarkerColorAlpha(2,0.6); G2->SetMarkerColorAlpha(4,0.6); G3->SetMarkerColorAlpha(1,0.6);
    G->SetMarkerStyle(20); G2->SetMarkerStyle(20);G3->SetMarkerStyle(20);
	
     TPaveText *Pt1,  *Pt22, *Pt3;
     Pt1 = new TPaveText(XMin+200,YMin+0.16,XMin+250,YMin+1.3,"NB");
     Pt1->AddText("#it{n}_{#tau}=4");Pt1->AddText("#it{n}_{#tau}=6");Pt1->AddText("#it{n}_{#tau}=8");
     Pt1->SetTextAlign(21);Pt1->SetBorderSize(0);Pt1->SetFillStyle(0);  Pt1->SetTextFont(22);Pt1->SetTextSize(0.06);  Pt1->Draw();

     Pt22 = new TPaveText(XMin+250,YMin+1.2,XMin+350,YMin+2.,"NB");
     Pt22->AddText("#it{FF}/(#it{T}^{4}#it{u}^{4}_{0})");
     Pt22->SetTextAlign(21);Pt22->SetBorderSize(0);Pt22->SetFillStyle(0);  Pt22->SetTextFont(22); Pt22->SetTextSize(0.06);
     Pt22->Draw();
     //Pt3 = new TPaveText(XMin+220,YMax-1.2,XMin+330,YMax-0.2,"NB");  
     //Pt3->AddText("#it{Z}^{(3)}_{T}#it{FF}/#it{T}^{4}");
     //Pt3->SetTextAlign(21);Pt3->SetBorderSize(0);Pt3->SetFillStyle(0);  Pt3->SetTextFont(22);Pt3->SetTextSize(0.048);
     //Pt3->Draw();
  
    TLEGEND3(0.42,0.18,0.47,0.38, G,"","p",  G2,"","p", G3,"","p", 0.06);
    //TLEGEND3(0.4,0.7,0.5,0.9, g,"","p",  g2,"","p", g3,"","p", 0.06);
    TLEGEND1(0.35,0.27,0.98,0.34, DummyH,"#it{R}_{SB} #times #it{s}/2#it{T}^{3} [HISQ]","l", 0.06);
    TLEGEND1(0.35,0.18,0.98,0.26, DummyEQCD,"#it{R}_{EQCD} #times #it{s}/2#it{T}^{3} [HISQ]","l", 0.06);
    
    
  sprintf(FileName,"Graph%s_Plot.pdf", QCDTypeS.c_str());
  //c0->SaveAs("MyFile.eps");
  //system("epstopdf MyFile.eps");char COMMAND[100000];
  //sprintf(COMMAND,"mv MyFile.pdf %s",FileName);
  //cout<<string(COMMAND)<<endl;
  //system(COMMAND);

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


TGraphErrors *RenormalizeFFoT4(TGraphErrors *g, TGraph *g2, string QCDTypeS)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();
  double x[N], y[N], ex[N], ey[N]; double g0Square=0;
  if(N!=N2){ cout<<"Error: N="<<N<<", N2="<<N2<<endl;}
  double Rx[N],Ry[N];
  for(int i=0;i<N;i++)
    {
      x[i]=X[i];
      ex[i]=EX[i];
      if(QCDTypeS =="FullQCD")   { g0Square = 10.0/X2[i]; }
      if(QCDTypeS =="PureGauge") { g0Square =  6.0/X2[i]; }
      Rx[i] = X[i];
      Ry[i] = RenormalizationConstant(g0Square,QCDTypeS);
      y[i]  = Y[i]*Ry[i];
      ey[i] = EY[i]*Ry[i];
      cout<<"T="<<x[i]<<", g0*g0="<<g0Square <<", y="<<y[i]<<endl;
    }
  TGraph *Rg = new TGraph(N,Rx,Ry);
  if(X[0]==201){Custom(Rg,2); Rg->Draw("AL"); Rg->GetYaxis()->SetTitle("Z^{(3)}_{T}"); Rg->GetXaxis()->SetTitle("#it{T} (MeV)");}
  if(X[0]==147){Custom(Rg,4); Rg->Draw("L p same");}
  if(X[0]==182){Custom(Rg,1); Rg->Draw("L p same");}
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

double RenormalizationConstant(double g0Square, string QCDTypeS)
{
  double Z6 = ((1.0 - 0.4367*g0Square)/(1.0-0.7074*g0Square)) - 0.0971*pow(g0Square,2) + 0.0886*pow(g0Square,3) - 0.2909*pow(g0Square,4);
  double zT = (1.0-0.509*g0Square)/(1.0-0.4789*g0Square);
  double Z3 =0.0;
  cout<<"Z6 = "<<Z6<<", zT="<<zT<<endl;
  Z3 = zT*Z6;
  return Z3;
}

TGraphErrors *DiscretizeHISQEntropy()
{
  TGraph *H = new TGraph("../EntropyDensity/EntropyDensity_HISQ_nf3.txt","%lf %*lf %*lf %*lf %*lf %*lf %*lf %lf %*lf");
  const int NH=12;
  double T[NH]={130,150,170,200, 250, 300,  400,  500,  600,  700,  800,  850};
  double sT[NH], EX[NH], EY[NH];

  for(int i=0; i<NH; i++)
    {
      sT[i] = H->Eval(T[i]);
      EX[i]=0;
      EY[i]=0;
    }
  TGraphErrors *ans = new TGraphErrors(NH,T,sT,EX,EY);
  return ans;
}

TGraphErrors *ConvertEQCD()
{  
  double TcFullQCD = 156;//MeV Prog.Part.Nucl.Phys. 70 (2013) 55-107
  double TcPureGauge = 265;//MeV Prog.Part.Nucl.Phys. 70 (2013) 55-107
  TGraph *H = new TGraph("../EntropyDensity/EntropyDensity_HISQ_nf3.txt","%lf %*lf %*lf %*lf %*lf %*lf %*lf %lf %*lf");
  TGraph *Herror = new TGraph("../EntropyDensity/EntropyDensity_HISQ_nf3.txt","%lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %lf");

  const int NH=12;
  double T[NH]={130,150,170,200, 250, 300,  400,  500,  600,  700,  800,  850};
  double TTc[NH], sT[NH], EXH[NH], EYH[NH];
  
  for(int i=0; i<NH; i++)
    {
      TTc[i] = T[i]/TcFullQCD;
      sT[i] = H->Eval(T[i]);
      EXH[i]=0;
      EYH[i]=Herror->Eval(T[i]);
    }

  string FileNameNf3="../EntropyDensity/EntropyDensity_EQCD_nf3.txt";
  string FileNameNf0="../EntropyDensity/EntropyDensity_EQCD_nf0.txt";
  TGraph *g  = new TGraph(FileNameNf0.c_str(), "%lf %lf");
  TGraph *g2 = new TGraph(FileNameNf0.c_str(), "%lf %*lf %lf");
  TGraph *g3 = new TGraph(FileNameNf0.c_str(), "%lf %*lf %*lf %lf");

  TGraph *G  = new TGraph(FileNameNf3.c_str(), "%lf %lf");
  TGraph *G2 = new TGraph(FileNameNf3.c_str(), "%lf %*lf %lf");
  TGraph *G3 = new TGraph(FileNameNf3.c_str(), "%lf %*lf %*lf %lf");
  double RY[NH], ERY[NH], ratio=0, low=0, high=0, temp=0;
  for(int i=0; i<NH; i++)
    {
      high = sT[i]*g->Eval(TTc[i])/G->Eval(TTc[i]);
      low =  sT[i]*g3->Eval(TTc[i])/G3->Eval(TTc[i]); 
      if(high < low){temp=low; low=high; high=temp;}
      cout<<"XH[i] = "<<T[i]<<",high="<<high<<", low="<<low<<endl;
      RY[i] = (high+low)/2.0;
      ERY[i] =  sqrt(pow(RY[i]*EYH[i]/sT[i],2) + pow((high-low)/2.0,2) );
      cout<<"RY="<<RY[i]<<", ERY="<<ERY[i]<<endl;
    }
  
  TGraphErrors *ans = new TGraphErrors(NH,T,RY,EXH, ERY);
  return ans;
}




