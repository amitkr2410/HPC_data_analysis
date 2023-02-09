#include"Custom.h"
#include"../Script/StatTool.h"
// define a function of calculating the ratio between two root graphs
TGraphErrors * ScaleGraphBeta(TGraphErrors *g, TGraph *g2, int PowerOfU0);
TGraphErrors *Divideg0Square(TGraphErrors *g, TGraph *g2, string QCDTypeS);
void FitMyGraph(TGraphErrors *g, TGraphErrors *g2, TGraphErrors *g3,double *Value);
void FitMyGraphCombined(TGraphErrors *g, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, double *Value);
void ContinuumExtrapolation()
{
  CustomGlobal();
  char FileName[10000];
  double XMin=150, XMax=1000, YMin=0, YMax=20;
  string QCDTypeS="PureGauge";
  //string QCDTypeS="FullQCD";
  double *MinMax; double Diff=0; int LW=2;
  gStyle->SetLineWidth(1);
  //TCanvas *c = new TCanvas("c","c",500,500);
  TGraphErrors *g, *G, *g2, *G2, *g3, *G3;
  TGraphErrors *h, *H, *h2, *H2, *h3, *H3;
  int PowersOfU0Clover=8, PowersOfU0FourLinks=4; int Nt, Ns;
  
  Nt=4; Ns=16;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());  
  g = new TGraphErrors(FileName,"%lf %lf %lf ");
  G = ScaleGraphBeta(g,T , PowersOfU0Clover);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  TGraphErrors *hh = new TGraphErrors(FileName,"%lf %lf %lf ");
  h = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  H = ScaleGraphBeta(h, T, PowersOfU0FourLinks);
  g = Divideg0Square(g,T,QCDTypeS);
  G = Divideg0Square(G,T,QCDTypeS);
  h = Divideg0Square(h,T,QCDTypeS);
  H = Divideg0Square(H,T,QCDTypeS);
  
  Nt=6;	Ns=24;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T2 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  g2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  G2 = ScaleGraphBeta(g2,T2 , PowersOfU0Clover);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  TGraphErrors *hh2 = new TGraphErrors(FileName,"%lf %lf %lf ");
  h2 = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  H2 = ScaleGraphBeta(h2, T2, PowersOfU0FourLinks);

  g2 = Divideg0Square(g2,T2,QCDTypeS);
  G2 = Divideg0Square(G2,T2,QCDTypeS);
  h2 = Divideg0Square(h2,T2,QCDTypeS);
  H2 = Divideg0Square(H2,T2,QCDTypeS);
  

  Nt=8;	Ns=32;
  sprintf(FileName,"%sDataPloopNt%d_Ns%d.txt", QCDTypeS.c_str(), Nt, Ns);
  TGraph *T3 = new  TGraph(FileName,"%lf %*lf %*lf %lf");
  sprintf(FileName,"Nt%d_%s.txt",Nt, QCDTypeS.c_str());
  g3 = new TGraphErrors(FileName,"%lf %lf %lf ");
  G3 = ScaleGraphBeta(g3,T3 , PowersOfU0Clover);
  sprintf(FileName,"%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);
  TGraphErrors *hh3 = new TGraphErrors(FileName,"%lf %lf %lf ");
  h3 = new TGraphErrors(FileName,"%lf %*lf %*lf  %lf %lf ");
  H3 = ScaleGraphBeta(h3, T3, PowersOfU0FourLinks);

  g3 = Divideg0Square(g3,T3,QCDTypeS);
  G3 = Divideg0Square(G3,T3,QCDTypeS);
  h3 = Divideg0Square(h3,T3,QCDTypeS);
  H3 = Divideg0Square(H3,T3,QCDTypeS);

  TCanvas *c = new TCanvas("c","c",500,500);
  c->Divide(3,2);
  double Value[2]; Value[0]=0; Value[1]=0;
  c->cd(1); FitMyGraph(g,g2,g3,Value);
  c->cd(2); FitMyGraph(G,G2,G3,Value);
  c->cd(4); FitMyGraph(h,h2,h3,Value);
  c->cd(5); FitMyGraph(H,H2,H3,Value);
  c->cd(6); FitMyGraphCombined(G,G2,G3,H,H2,H3,Value);
  /*
  TCanvas *c0 = new TCanvas("c0","c0",1600,900);
  Custom(g,1);Custom(g2,2); Custom(g3,4);
  Custom(G,1);Custom(G2,2); Custom(G3,4);
  Custom(h,1);Custom(h2,2);Custom(h3,4);
  Custom(H,1);Custom(H2,2);Custom(H3,4);
  g->Draw("Ap"); g2->Draw("p same");  g3->Draw("p same");
  G->Draw("Lp");  G2->Draw("Lp same");	G3->Draw("Lp same");
  h->Draw("p");  h2->Draw("p"); h3->Draw("p");
  H->Draw("Lp");  H2->Draw("Lp"); H3->Draw("Lp");
  g->SetLineWidth(LW);g2->SetLineWidth(LW); g3->SetLineWidth(LW);
  G->SetLineWidth(4);G2->SetLineWidth(4); G3->SetLineWidth(4);
  G->SetLineStyle(5);G2->SetLineStyle(5); G3->SetLineStyle(5);
  G->SetMarkerStyle(24);G2->SetMarkerStyle(24); G3->SetMarkerStyle(24);

  h->SetLineWidth(LW); h2->SetLineWidth(LW); h3->SetLineWidth(LW);
  h->SetMarkerStyle(21); h2->SetMarkerStyle(21);  h3->SetMarkerStyle(21);
  H->SetLineWidth(4); H2->SetLineWidth(4);  H3->SetLineWidth(4);
  H->SetMarkerStyle(25); H2->SetMarkerStyle(25); H3->SetMarkerStyle(25);
  H->SetLineStyle(5); H2->SetLineStyle(5); H3->SetLineStyle(5);

  
  g->GetYaxis()->SetTitle("<#it{O}>");g->GetXaxis()->SetTitle("#it{T} (MeV)");
  MinMax=FindMinMax4(G, G2, G3, H2); Diff = (MinMax[1]-MinMax[0])/6.0;
  g->GetYaxis()->SetRangeUser(YMin,YMax);
  g->GetXaxis()->SetRangeUser(XMin, XMax);
  TLEGEND3(0.12,0.7,0.3,0.92, g,"#it{g}^{2}_{0}FF/T^{4}; #it{n}_{#tau}=4","lp",  g2,"#it{g}^{2}_{0}FF/T^{4}; #it{n}_{#tau}=6","lp", g3,"#it{g}^{2}_{0}FF/T^{4}; #it{n}_{#tau}=8","lp", 0.04);
  TLEGEND3(0.32,0.7,0.4,0.92,  G,"#it{g}^{2}_{0}FF/(T^{4}#it{u}^{8}_{0})","lp",  G2,"#it{g}^{2}_{0}FF/(T^{4}#it{u}^{8}_{0})","lp", G3,"#it{g}^{2}_{0}FF/(T^{4}#it{u}^{8}_{0})","lp", 0.04);
  TLEGEND3(0.475,0.7,0.65,0.92, h,"(1x1+1x2)#it{g}^{2}_{0}FF/T^{4};","lp", h2,"(1x1+1x2)#it{g}^{2}_{0}FF/T^{4};","lp",h3,"(1x1+1x2)#it{g}^{2}_{0}FF/T^{4};","lp",   0.04);
  TLEGEND3(0.7,0.7,0.95,0.92, H,"(1x1+1x2)#it{g}^{2}_{0}FF/(T^{4}#it{u}^{4}_{0})","lp", H2, "(1x1+1x2)#it{g}^{2}_{0}FF/(T^{4}#it{u}^{4}_{0})","lp", H3,"(1x1+1x2)#it{g}^{2}_{0}FF/(T^{4}#it{u}^{4}_{0})","lp", 0.04);
  g->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetMargin(0.08,0.02,0.12,0.02);
  */
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

void FitMyGraph(TGraphErrors *g, TGraphErrors *g2, TGraphErrors *g3,double *Value)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  const int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  double x[3], y[3], ex[3], ey[3]; int Size[3]={N,N2,N3}, count=0;
  const int Nt[3]={4,6,8};
  for(int k=0;k<3;k++)
    {
      x[k]=1.0/(Nt[k]*Nt[k]); 
      y[k]=0; ex[k]=0; ey[k]=0;
      count=0;
      for(int i=0;i<Size[k];i++)
	{
	  if(k==0 && X[i] > 400.0){y[k]= y[k]+Y[i];  ey[k]= ey[k] + pow(EY[i],2);  count++;}
	  if(k==1 && X[i] > 400.0){y[k]= y[k]+Y2[i]; ey[k]= ey[k] + pow(EY2[i],2); count++;}
	  if(k==2 && X[i] > 400.0){y[k]= y[k]+Y3[i]; ey[k]= ey[k] + pow(EY3[i],2); count++;}
	}
      y[k]=y[k]/count; ey[k] = sqrt(ey[k]/count);
      cout<<"(y,ey)="<<y[k]<<","<<ey[k]<<endl;
    }
  
  TGraphErrors *gr = new TGraphErrors(3,x,y,ex,ey);
  Custom(gr,2);gr->Draw("Ap");
  gr->GetXaxis()->SetLimits(-0.01,0.08);
  TF1 *f = new TF1("f","[0]+[1]*x +[2]*x*x",-0.01,1.08);
  f->SetParameter(0,1.0);
  f->SetParameter(1,1.0);
  f->SetParameter(2,1.0);
  gr->Fit(f,"R");
  cout<<"(NDF,ChiSquare)="<<f->GetNDF()<<","<<f->GetChisquare()<<"\n"<<endl;
}


void FitMyGraphCombined(TGraphErrors *g, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, double *Value)
{
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  const int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  const int N4=g4->GetN(); double *X4=g4->GetX();double *Y4=g4->GetY();double *EY4=g4->GetEY(); double *EX4=g4->GetEX();
  const int N5=g5->GetN(); double *X5=g5->GetX();double *Y5=g5->GetY();double *EY5=g5->GetEY(); double *EX5=g5->GetEX();
  const int N6=g6->GetN(); double *X6=g6->GetX();double *Y6=g6->GetY();double *EY6=g6->GetEY(); double *EX6=g6->GetEX();
  double x[6], y[6], ex[6], ey[6]; int Size[6]={N,N2,N3,N4, N5,N6}, count=0;
  const double Nt[6]={4,6,8,4.001,6.001,8.001};
  for(int k=0;k<6;k++)
    {
      x[k]=1.0/(Nt[k]*Nt[k]);
      y[k]=0; ex[k]=0; ey[k]=0;
      count=0;
      for(int i=0;i<Size[k];i++)
        {
          if(k==0 && X[i] > 400.0){y[k]= y[k]+Y[i];  ey[k]= ey[k] + pow(EY[i],2);  count++;}
          if(k==1 && X[i] > 400.0){y[k]= y[k]+Y2[i]; ey[k]= ey[k] + pow(EY2[i],2); count++;}
          if(k==2 && X[i] > 400.0){y[k]= y[k]+Y3[i]; ey[k]= ey[k] + pow(EY3[i],2); count++;}
	  if(k==3 && X[i] > 400.0){y[k]= y[k]+Y4[i]; ey[k]= ey[k] + pow(EY4[i],2);  count++;}
          if(k==4 && X[i] > 400.0){y[k]= y[k]+Y5[i]; ey[k]= ey[k] + pow(EY5[i],2); count++;}
          if(k==5 && X[i] > 400.0){y[k]= y[k]+Y6[i]; ey[k]= ey[k] + pow(EY6[i],2); count++;}
        }
      y[k]=y[k]/count; ey[k] = sqrt(ey[k]/count);
      cout<<"(y,ey)="<<y[k]<<","<<ey[k]<<endl;
    }

  TGraphErrors *gr = new TGraphErrors(6,x,y,ex,ey);
  Custom(gr,2);gr->Draw("Ap");
  gr->GetXaxis()->SetLimits(-0.01,0.08);
  TF1 *f = new TF1("f","[0] + [1]*x + [2]*x*x",-0.01,1.08);
  f->SetParameter(0,1.0);
  f->SetParameter(1,1.0);
  f->SetParameter(2,1.0); 
  gr->Fit(f,"R");
  cout<<"(NDF,ChiSquare)="<<f->GetNDF()<<","<<f->GetChisquare()<<"\n\n";
}

