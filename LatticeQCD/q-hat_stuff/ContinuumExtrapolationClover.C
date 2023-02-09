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
double StrongCouplingOneLoop(string QCDTypeS, double muSquare);
TGraphErrors *EntropyFromRawDataPureSU3();
TGraphErrors *EntropyFromRawDataFullQCD();
void ContinuumExtrapolationClover()
{
  CustomGlobal();
  char FileName[10000];
  double XMin=80, XMax=1020, YMin=0, YMax=6;
  string QCDTypeS="PureGauge";
  //string QCDTypeS="FullQCD";
  if(QCDTypeS=="PureGauge"){XMin=150; XMax=1020;}
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
  g = Divideg0Square(g,T,QCDTypeS);
  R = Renormalize(g,T,QCDTypeS);
  G = ScaleGraphBeta(g,T , PowersOfU0Clover);
  
  Nt=6;	Ns=24;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T2 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  g2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  g2 = Divideg0Square(g2,T2,QCDTypeS);
  R2 = Renormalize(g2,T2,QCDTypeS);
  G2 = ScaleGraphBeta(g2,T2 , PowersOfU0Clover);
 
  Nt=8;	Ns=32;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T3 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  g3 = new TGraphErrors(FileName,"%lf %lf %lf ");
  g3 = Divideg0Square(g3,T3,QCDTypeS);
  R3 = Renormalize(g3,T3,QCDTypeS);
  G3 = ScaleGraphBeta(g3,T3 , PowersOfU0Clover);
  
  TGraphErrors *CG; if(QCDTypeS=="PureGauge"){ CG=FitMyGraph(R,R2,R3);}
  if(QCDTypeS=="PureGauge"){ ComputeQhat(CG,QCDTypeS);}
  if(QCDTypeS=="FullQCD"){   ComputeQhat(G2,QCDTypeS);}
  
  TCanvas *c0 = new TCanvas("c0","c0",350,300);
  TGraphErrors *S, *S2;
  if(QCDTypeS=="PureGauge"){ S = EntropyFromRawDataPureSU3();}
  if(QCDTypeS=="FullQCD"){ S = EntropyFromRawDataFullQCD();}
  S2 = ScaleTGraphErrors(S,0.35);
  S = ScaleTGraphErrors(S,0.5);
  
  Custom(g,2);Custom(g2,4); Custom(g3,1);  
  Custom(R,2);Custom(R2,4);Custom(R3,1);
  Custom(G,2);Custom(G2,4); Custom(G3,1);
  g->Draw("Ap"); g2->Draw("p same");  g3->Draw("p same");  
  R->Draw("p");  R2->Draw("p"); R3->Draw("p");
  G->Draw("p"); G2->Draw("p same");  G3->Draw("p same");
  if(QCDTypeS=="PureGauge"){ Custom(CG,8);CG->Draw("L"); 
  Custom(S,6); S->Draw("L");
}
  
  g->SetLineWidth(LW);g2->SetLineWidth(LW); g3->SetLineWidth(LW);
  R->SetMarkerStyle(29);R2->SetMarkerStyle(29); R3->SetMarkerStyle(29);
  R->SetMarkerSize(2);R2->SetMarkerSize(2);R3->SetMarkerSize(2);
  G->SetMarkerStyle(24);G2->SetMarkerStyle(24); G3->SetMarkerStyle(24);
  S->SetLineWidth(5);
  if(QCDTypeS=="PureGauge"){ CG->SetLineWidth(6); }

  
  g->GetYaxis()->SetTitle("<#it{O}>");g->GetXaxis()->SetTitle("#it{T} (MeV)");
  MinMax=FindMinMax4(g, g2, g3, g2); Diff = (MinMax[1]-MinMax[0])/6.0;
  g->GetYaxis()->SetRangeUser(YMin,YMax);
  g->GetXaxis()->SetLimits(XMin, XMax);
  TPaveText *Pt1, *Pt2, *Pt22, *Pt3; 
  if(QCDTypeS=="FullQCD"){Pt1 = new TPaveText(XMin+60,YMax-2.0,XMin+150,YMax-0.53,"NB");}
  if(QCDTypeS=="PureGauge"){Pt1 = new TPaveText(XMin+40,YMax-2.0,XMin+100,YMax-0.53,"NB");}
  Pt1->AddText("#it{n}_{#tau}=4");Pt1->AddText("#it{n}_{#tau}=6");Pt1->AddText("#it{n}_{#tau}=8");
  Pt1->SetTextAlign(21);Pt1->SetBorderSize(0);Pt1->SetFillStyle(0);  Pt1->SetTextFont(22);Pt1->SetTextSize(0.05);  Pt1->Draw();

  if(QCDTypeS=="FullQCD"){Pt2 = new TPaveText(XMin+170,YMax-0.8,XMin+330,YMax-0.1,"NB");}
  if(QCDTypeS=="PureGauge"){Pt2 = new TPaveText(XMin+80,YMax-0.8,XMin+230,YMax-0.1,"NB");}
  Pt2->AddText("#it{FF}/#it{T}^{4}");
  Pt2->SetTextAlign(21);Pt2->SetBorderSize(0);Pt2->SetFillStyle(0);  Pt2->SetTextFont(22);Pt2->SetTextSize(0.042);  Pt2->Draw();

  if(QCDTypeS=="FullQCD"){ Pt22 = new TPaveText(XMin+370,YMax-0.8,XMin+450,YMax-0.1,"NB");}
  if(QCDTypeS=="PureGauge"){Pt22 = new TPaveText(XMin+260,YMax-0.8,XMin+320,YMax-0.1,"NB");}
  Pt22->AddText("#it{FF}/(#it{T}^{4}#it{u}^{4}_{0})");
  Pt22->SetTextAlign(21);Pt22->SetBorderSize(0);Pt22->SetFillStyle(0);  Pt22->SetTextFont(22); Pt22->SetTextSize(0.042);  Pt22->Draw();
  
  if(QCDTypeS=="FullQCD"){  Pt3 = new TPaveText(XMin+370,YMax-0.8,XMin+450,YMax-0.1,"NB");  }
  if(QCDTypeS=="PureGauge"){ Pt3 = new TPaveText(XMin+420,YMax-0.8,XMin+500,YMax-0.1,"NB");  }
  Pt3->AddText("#it{Z}^{(3)}_{T}#it{FF}/#it{T}^{4}");
  Pt3->SetTextAlign(21);Pt3->SetBorderSize(0);Pt3->SetFillStyle(0);  Pt3->SetTextFont(22);Pt3->SetTextSize(0.042);
  if(QCDTypeS=="PureGauge"){ Pt3->Draw();}

  if(QCDTypeS=="FullQCD"){
    TLEGEND3(0.32,0.7,0.42,0.92, g,"","p",  g2,"","p", g3,"","p", 0.04);
    TLEGEND3(0.475,0.7,0.65,0.92, G,"","p", G2,"","p",G3,"","p",   0.04);
  }
  if(QCDTypeS=="PureGauge"){
    TLEGEND3(0.23,0.7,0.33,0.92, g,"","p",  g2,"","p", g3,"","p", 0.04);
    TLEGEND3(0.36,0.7,0.46,0.92, G,"","p", G2,"","p",G3,"","p",   0.04);
    TLEGEND3(0.50,0.7,0.54,0.92, R,"","p", R2,"","p",R3,"","p",   0.04);
    TLEGEND1(0.56,0.7,0.83,0.92, CG, "Continuum [#it{Z}^{(3)}_{T}#it{FF}/#it{T}^{4}]","lp", 0.045);
  //TLEGEND1(0.26,0.1,0.83,0.32, S, "0.5#it{s}/#it{T}^{3} [LNP 583 209 (2002)]","l", 0.05);
    TLEGEND1(0.26,0.2,0.83,0.28, S, "#it{s}/2#it{T}^{3} [PLB 769 385 (2017)]","l", 0.05);
  }
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

TGraphErrors *Renormalize(TGraphErrors *g,TGraph *g2, string QCDTypeS)
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
      
      y[i]  = Y[i]*RenormalizationConstant(g0Square,QCDTypeS);
      ey[i] = EY[i]*RenormalizationConstant(g0Square,QCDTypeS);
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

double RenormalizationConstant(double g0Square, string QCDTypeS)
{
  double Z6 = ((1.0 - 0.4367*g0Square)/(1.0-0.7074*g0Square)) - 0.0971*pow(g0Square,2) + 0.0886*pow(g0Square,3) - 0.2909*pow(g0Square,4);
  double zT = (1.0-0.509*g0Square)/(1.0-0.4789*g0Square);
  double Z3 =0.0;

  Z3 = zT*Z6;
  return Z3;
}

TGraphErrors *FitMyGraph(TGraphErrors *g, TGraphErrors *g2, TGraphErrors *g3)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  const int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  TGraph *IG, *IG2, *IG3;
  IG  = g;//new (TGraph)g->Clone();
  IG2 = g2;//new (TGraph)g2->Clone();
  double RX[N], RY[N], REX[N], REY[N]; double CX[1], CY[1], CEX[1], CEY[1];
  double x[3], y[3], ex[3], ey[3], Lx[2], Ly[2], Lex[2], Ley[2]; int Size[3]={N,N2,N3}, count=0;
  const int Nt[3]={4,6,8}; char VarName[1000];
  TCanvas *CFit = new TCanvas("CFit","CFit",1000,700);
  CFit->Divide(3,3);
  for(int i=0; i<N3; i++)
    {      
      for(int k=0;k<3;k++)
	{ x[k]=0; y[k]=0; ex[k]=0; ey[k]=0;
	}
      for(int k=0;k<3;k++)
        {
	  x[k]=1.0/(Nt[k]*Nt[k]); 
	  if(k==0){ y[k]=IG->Eval(X3[i]);  ey[k]=FindErrorBar(g,X3[i]); }
	  if(k==1){ y[k]=IG2->Eval(X3[i]); ey[k]=FindErrorBar(g2,X3[i]);}
	  if(k==2){ y[k]=Y3[i]; ey[k]=EY3[k];}	  
	  cout<<"(T,y,ey)="<<X3[i]<<","<<y[k]<<","<<ey[k]<<endl;
	}    
      CFit->cd(i+1);
      TGraphErrors *gr = new TGraphErrors(3,x,y,ex,ey);
      Lx[0]=x[1];Lex[0]=ex[1];Ly[0]=y[1];Ley[0]=ey[1];
      Lx[1]=x[2];Lex[1]=ex[2];Ly[1]=y[2];Ley[1]=ey[2];
      TGraphErrors *Lgr = new TGraphErrors(2,Lx,Ly,Lex,Ley);
      Custom(gr,2);gr->Draw("Ap"); Custom(Lgr,4); Lgr->Draw("p");
      gr->SetMarkerSize(1);  Lgr->SetMarkerSize(1);
      gr->GetXaxis()->SetLimits(-0.01,0.08); gr->GetYaxis()->SetRangeUser(y[0]-1.4,y[0]+0.5);
      TF1 *f = new TF1("f","[0]+[1]*x +[2]*x*x",-0.01,1.08);
      TF1 *Lf = new TF1("Lf","[0]+[1]*x",-0.01,1.08); Lf->SetLineColor(2);
      f->SetParameter(0,1.0);Lf->SetParameter(0,1.0);
      f->SetParameter(1,1.0);Lf->SetParameter(1,1.0);
      f->SetParameter(2,1.0);
      gr->Fit(f,"R"); Lgr->Fit(Lf,"R");

      RX[i] = X3[i];
      REX[i] = 0;
      RY[i] = Lf->GetParameter(0); // (f->GetParameter(0) + Lf->GetParameter(0))/2.0;
      REY[i] = std::sqrt(  pow(Lf->GetParError(0),2.0) + pow(0.5*(f->GetParameter(0) - Lf->GetParameter(0)),2.0) );
      //REY[i] = std::sqrt(  pow(f->GetParError(0),2.0)  );
      CX[0]=0; CY[0]=RY[i]; CEX[0]=REX[i]; CEY[0]=REY[i];
      TGraphErrors *Cg = new TGraphErrors(1,CX,CY,CEX,CEY);
      Custom(Cg,1);Cg->Draw("p same");
      sprintf(VarName,"#it{T}=%.1f MeV",X3[i]);
      TLEGEND1(0.12,0.7,0.3,0.92, gr,VarName,"l", 0.07);
      cout<<"(RX,RY,REY)="<<RX[i]<<","<< RY[i] <<"," <<REY[i]<<") : ("<<Lf->GetParameter(0)<<","<<Lf->GetParError(0) <<"\n"<<endl;
      cout<<"(i,NDF,ChiSquare)="<<i<<","<<f->GetNDF()<<","<<f->GetChisquare()<<"\n"<<endl;

      
    }
  
  TGraphErrors *Ans = new TGraphErrors(N,RX,RY,REX,REY);
  return Ans;
}

double FindErrorBar(TGraphErrors *g, double a)
{
 const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
 double ey=0;
 double YHigh[N];
 for(int i=0; i<N; i++)
   {
     YHigh[i] = Y[i] + EY[i]; 
   }
 TGraph *GHigh = new TGraph(N,X,YHigh);
 TGraph *G = new TGraph(N,X,Y);
 ey = GHigh->Eval(a) - G->Eval(a);
 if(ey <0) {ey=0;cout<<"ey = 0 ........"<<endl;}
 return ey;
}

void ComputeQhat(TGraphErrors *g, string QCDTypeS)
{
   const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
   //PreFactor = 8*\sqrt(2)*pi*alphas/(Nc*2*sqrt(2)*2); for <M| FF/T^4 |M>
   double x[N], y[N], ex[N], ey[N]; double Nc=3, PreFactor=0, AlphaSMid=0, AlphasErr=0;
   double muSquare1=0, muSquare2=0, qhatHigh=0;
   ofstream fout;
   if(QCDTypeS=="PureGauge"){ fout.open("Qhat_Pure_SU3_Clover.txt",ios::out); }
   if(QCDTypeS=="FullQCD"){ fout.open("Qhat_Full_QCD_Clover.txt",ios::out); }
   fout<<"#T(GeV) \t q-hat/T^3 \t ErrX    \t ErrY   ErrYAlphas"<<endl;
   for(int i=0; i<N; i++)
     {
       x[i] = X[i];
       ex[i]=0;
       muSquare1=pow(3.1415926*2.0*X[i],2);
       muSquare2=pow(3.1415926*4.0*X[i],2);
       AlphaSMid = (StrongCouplingOneLoop(QCDTypeS,muSquare1) + StrongCouplingOneLoop(QCDTypeS,muSquare2))/2.0;
       AlphasErr = fabs(StrongCouplingOneLoop(QCDTypeS,muSquare1) - StrongCouplingOneLoop(QCDTypeS,muSquare2))/2.0;
       PreFactor = 8.0*std::sqrt(2)*3.1415*AlphaSMid/(Nc*2*std::sqrt(2)*2.0);
       y[i]= Y[i]*PreFactor;
       qhatHigh = Y[i]*PreFactor*(AlphaSMid+AlphasErr)/AlphaSMid;
       
       ey[i] = sqrt( pow(EY[i]*PreFactor,2) +  pow(qhatHigh - y[i],2) );
       cout<<"\n \n PreFactor="<<PreFactor<<endl<<"T=";
       cout<<x[i]/1000.0<<"\t"<<y[i]<<"\t"<<ex[i]<<"\t"<<ey[i]<<endl;
       //fout<<x[i]/1000.0<<"\t"<<y[i]<<"\t"<<ex[i]<<"\t"<<ey[i]<<endl;}
        fout<<x[i]/1000.0<<"\t"<<y[i]<<"\t"<<ex[i]<<"\t"<<EY[i]*PreFactor<<"\t"<<fabs(qhatHigh - y[i])<<endl;
     }
   fout.close();
   TGraphErrors *Ans = new TGraphErrors(N,x,y,ex,ey);
   TCanvas *qhat = new TCanvas("qhat","qhat",500,450);
   Custom(Ans,2);Ans->Draw("AL");Ans->SetLineWidth(3);
   Ans->GetYaxis()->SetTitle("#hat{#it{q}}/#it{T}^{3}");
   Ans->GetXaxis()->SetTitle("#it{T} (MeV)");
   Ans->GetYaxis()->SetRangeUser(0.0,6.0);
   Ans->GetXaxis()->SetRangeUser(150.0, 1000);
   Ans->GetYaxis()->SetTitleOffset(1.0);
   Ans->GetXaxis()->SetNdivisions(310);
   gPad->SetMargin(0.12,0.02,0.12,0.02);
}

double StrongCouplingOneLoop(string QCDTypeS, double muSquare)
{
  int Nf=0; 
  double QSquare = muSquare;//pow(3.1415926*2.0*T, 2.0);// MeV^2
  double LambdaSquare = pow(340, 2.0); // in MS-bar scheme, MeV^2 https://arxiv.org/pdf/1604.08082.pdf
  if(QCDTypeS == "FullQCD" ){Nf=3;}
  else {Nf=0; LambdaSquare= pow(243,2);} //https://inspirehep.net/literature/554365
  double beta0= 11 - (2.0*Nf/3.0); 
  double alphas = 4*3.1415926/(beta0*log(QSquare/LambdaSquare));cout<<"alphas="<<alphas<<"\t";
  return alphas;
}

TGraphErrors * EntropyFromRawDataPureSU3()
{
  //TGraph *g1 = new TGraph("../EntropyDensity/EntropyDensity_RawData_SU3.txt","%lf %lf ");
  TGraphErrors *g1 = new TGraphErrors("../EntropyDensity/EntropyDensity_RawData_SU3_MovingFrame.txt","%lf %lf %lf");
  const int N1=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(); double *EY1=g1->GetEY();
  int N3; double X3[N1], Y3[N1], Ex3[N1], Ey3[N1];
  for (int i=0; i<N1; i++)
    {
      X3[i] = X1[i]*265;
      N3=N1;
      cout<<"X1[i]="<<X1[i]<<"\tX3[i]="<<X3[i]<<endl;
      cout<<"Y1[i]="<<Y1[i]<<endl;
      //Y3[i] = Y1[i]*4.0/3.0;      Ey3[i] =0.0;      Ex3[i] =0.0;
      Y3[i] = Y1[i];      Ey3[i] =EY1[i];      Ex3[i] =0.0; 
    }
  TGraphErrors *g3 = new TGraphErrors(N3,X3,Y3,Ex3,Ey3);
  return g3;
}

TGraphErrors *EntropyFromRawDataFullQCD()
{
  TGraphErrors *g = new TGraphErrors("../EntropyDensity/EntropyDensity_RawData_QCD.txt","%lf %lf");
  return g;
}
