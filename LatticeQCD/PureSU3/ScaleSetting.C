#include "Custom.h"

TGraphErrors * Subtract(TGraphErrors *g1,TGraphErrors *g2);
double AsymptoticBetaFunction(double beta);
double KarshRBetaFunction(double beta);
double KarschLatticeSpacing(double beta);
double AmitLatticeSpacing(double beta, int Nt);
double AmitLatticeSpacingNEW(double beta);
double LatticeSpacing(double beta);
TGraphErrors * ConvertToT(TGraphErrors *g1, int Nt, double Norm);
void ScaleSetting()
{
  int SU3_Polyakov = 1;
  CustomGlobal();
  int ANS=AmitLatticeSpacingNEW(6);
  
      TGraphErrors *gT1 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_2_Ns8_I2100.txt","%lf %lf %lf");
      TGraphErrors *gT2 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_4_Ns16_I2100.txt","%lf %lf %lf");
      TGraphErrors *gT3 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_6_Ns24_I2100.txt","%lf %lf %lf");
      TGraphErrors *gT4 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_8_Ns32_I2100.txt","%lf %lf %lf");
      TGraphErrors *gT5 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_10_Ns40_I2100.txt","%lf %lf %lf");

      TGraphErrors *gV1 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_8_Ns8_I2100.txt","%lf %lf %lf");
      TGraphErrors *gV2 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_16_Ns16_I2100.txt","%lf %lf %lf");
      TGraphErrors *gV3 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_24_Ns24_I2100.txt","%lf %lf %lf");
      TGraphErrors *gV4 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_32_Ns32_I2100.txt","%lf %lf %lf");
      TGraphErrors *gV5 = new TGraphErrors("~/Dropbox/Lattice/PolyakovLoopSU3/result/Data_AbsAvgPolyakovLoopSU3_Nt_40_Ns40_I2100.txt","%lf %lf %lf");
      
      TCanvas *c = new TCanvas("c","c",600,600);

      Custom(gT2,4);Custom(gT3,1);Custom(gT4,6);Custom(gT5,2);
      Custom(gV2,7);Custom(gV3,40);Custom(gV4,42);Custom(gV5,3);
      gT2->Draw("APL");gT3->Draw("PL Same");gT4->Draw("PL Same");gT5->Draw("PL Same");
      gV2->Draw("PL");gV3->Draw("PL Same");gV4->Draw("PL Same");gV5->Draw("PL Same");

      gT2->GetXaxis()->SetLimits(2.,10.0);
      gT2->GetYaxis()->SetRangeUser(-0.5,3);
      gT2->GetYaxis()->SetTitle("< |P| >");
      gT2->GetXaxis()->SetTitle("#beta_{0}=6/g^{2}");

      gPad->SetMargin(0.14,0.05,0.15,0.1);
      TLegend *leg = new TLegend(0.20,0.4,0.45,0.89);
      leg->SetFillStyle(0);
      //leg->AddEntry(gT1," 2x8^3","leq");
      leg->AddEntry(gT2," 4x16^3","leq");
      leg->AddEntry(gT3," 6x24^3","leq");
      leg->AddEntry(gT4," 8x32^3","leq");
      leg->AddEntry(gT5," 10x40^3","leq");

      //leg->AddEntry(gV1," 8x8^3","leq");
      leg->AddEntry(gV2," 16x16^3","leq");
      leg->AddEntry(gV3," 24x24^3","leq");
      leg->AddEntry(gV4," 32x32^3","leq");
      leg->AddEntry(gV5," 40x40^3","leq");

      leg->SetTextFont(22);leg->SetTextSize(0.05);
      leg->Draw();
      //c->SaveAs("/home/amit/Dropbox/ConferencePlotsHardProbes2018/SU3_PloyakovNt4_Nt6_Nt8_Nt10I2100.png");
      //Subtracting the two graph (Thermal -Vacuum)
      TCanvas *c0 = new TCanvas("c0","c0",600,600);
      TGraphErrors *g1, *g2, *g3, *g4, *g5;
      //g1=Subtract(gT1, gV1);
      g2=Subtract(gT2, gV2);
      g3=Subtract(gT3, gV3);
      g4=Subtract(gT4, gV4);
      g5=Subtract(gT5, gV5);

      Custom(g2,4);Custom(g3,1); Custom(g4,6);Custom(g5,2);
      g2->Draw("APL");g3->Draw("PL Same");g4->Draw("PL Same");g5->Draw("PL Same");
      
      g2->GetXaxis()->SetLimits(3.,10.0);
      g2->GetYaxis()->SetRangeUser(-0.5,3);
      g2->GetYaxis()->SetTitle("< |P| >");
      g2->GetXaxis()->SetTitle("#beta_{0}=6/g^{2}");

      gPad->SetMargin(0.2,0.05,0.15,0.1);
      TLegend *leg1 = new TLegend(0.23,0.61,0.5,0.85);
      leg1->SetFillStyle(0);
      //leg1->AddEntry(g1," 2x8^3 (T-V)","leq");
      leg1->AddEntry(g2," 4x16^3(T-V)","leq");
      leg1->AddEntry(g3," 6x24^3(T-V)","leq");
      leg1->AddEntry(g4," 8x32^3(T-V)","leq");
      leg1->AddEntry(g5," 10x40^3(T-V)","leq");

      leg1->SetTextFont(22);leg1->SetTextSize(0.05);
      leg1->Draw();

      //c0->SetGrid();
      //c0->SaveAs("/home/amit/Dropbox/ConferencePlotsHardProbes2018/SU3_PloyakovTminusVNt4_Nt6_Nt8_Nt10I2100.png");
      //c0->SaveAs("SU3_PloyakovNt2_Ns8_I2100.pdf");
            
      //Perturbative conversion to T
      TCanvas *c1 = new TCanvas("c1","c1",600,500);
      double Norm =1.0;
      TGraphErrors *G1, *G2, *G3, *G4;
      //G1=ConvertToT(g1,  2,  1.26);//1.26
      G2=ConvertToT(g2,  4,  1.0);//1.28, 1.38
      G3=ConvertToT(g3,  6,  1.0);//1.5, 1.65
      G4=ConvertToT(g4,  8,  1.0);//1.76, 1.7
      //G5=ConvertToT(g5,  10,  1.8);//1.75
      Custom(G2,4);Custom(G3,1);Custom(G4,2);//Custom(G5,2);
      G2->Draw("ALP");G3->Draw("LP");G4->Draw("LP");//G5->Draw("LP");
      G2->SetMarkerSize(1.6);G3->SetMarkerSize(2);G4->SetMarkerSize(1.6);
      G2->SetMarkerStyle(20);G3->SetMarkerStyle(22);G4->SetMarkerStyle(21);
      G2->GetXaxis()->SetLimits(100.0,800.0);
      G2->GetYaxis()->SetRangeUser(-0.5,2.0);
      G2->GetYaxis()->SetTitle("< |P| >");
      G2->GetXaxis()->SetTitle("T (MeV)");

      gPad->SetMargin(0.2,0.05,0.15,0.1);
      TLegend *leg2 = new TLegend(0.53,0.61,0.8,0.85);
      leg2->SetFillStyle(0);
      //leg2->AddEntry(g1," 2x8^3","leq");
      leg2->AddEntry(G2," 4x16^3","lp");
      leg2->AddEntry(G3," 6x24^3","lp");
      leg2->AddEntry(G4," 8x32^3","lp");
      //leg2->AddEntry(G5," 10x40^3","leq");

      leg2->SetTextFont(22);leg2->SetTextSize(0.05);
      leg2->Draw();
      TLine *Line1 =new TLine(265.0,-0.5,265.0,2);
      Line1->SetLineColor(8); Line1->SetLineStyle(10);Line1->SetLineWidth(4);
      Line1->Draw();

      //c1->SaveAs("/home/amit/Dropbox/ConferencePlotsHardProbes2018/SU3_PloyakovNt4_Nt6_Nt8wrt_nonperturbativeRenormalization_I2100.png");
      //c0->SaveAs("SU3_PloyakovNt2_Ns8_I2100.pdf"); 
      
      //Lattice beta function Karsch and Amit
      const int Num=20; double Beta[Num], BETA[Num-2];
      double a0Karsch[Num], RbetaKarsch[Num-2];
      double a0Amit[Num],   RbetaAmit[Num-2];
      double a0AmitNEW[Num],   RbetaAmitNEW[Num-2]; 
      Beta[0]=5.0;Beta[Num-1]=7.0; int NT=8;
      for(int i=0;i<Num;i++)
	{ Beta[i]=Beta[0] + (i*(Beta[Num-1]-Beta[0])/(Num-1));
	  a0Karsch[i] =  KarschLatticeSpacing(Beta[i]);
	  a0Amit[i] =  AmitLatticeSpacing(Beta[i],NT);
	  a0AmitNEW[i]= AmitLatticeSpacingNEW(Beta[i]);
	}
      for(int i=1;i<Num-1;i++)
	{
	  BETA[i-1]= Beta[i];
	  RbetaKarsch[i-1]=   KarshRBetaFunction(BETA[i-1]);
	  RbetaAmit[i-1]=   -a0Amit[i]*(Beta[i+1]-Beta[i-1])/( a0Amit[i+1]-a0Amit[i-1] );
	  RbetaAmitNEW[i-1]=   -a0AmitNEW[i]*(Beta[i+1]-Beta[i-1])/( a0AmitNEW[i+1]-a0AmitNEW[i-1] );
	  //cout<<"beta ="<<BETA[i-1]<<", RbetaKarsch="<<RbetaKarsch[i-1]<<", RbetaAmit="<<RbetaAmit[i-1]<<endl;
	}
      TCanvas *c2 = new TCanvas("c2","c2",600,600);
      c2->Divide(2,1);
      c2->cd(1);
      TGraph *Gk = new TGraph(Num-2,BETA,RbetaKarsch);
      TGraph *Ga = new TGraph(Num-2,BETA,RbetaAmit);
      TGraph *GaNEW = new TGraph(Num-2,BETA,RbetaAmitNEW);
      Custom(Ga,2);Custom(GaNEW,4);Custom(Gk,1);
      Ga->Draw("ALP");GaNEW->Draw("LP same"); Gk->Draw("LP same");
      c2->cd(2);
      TGraph *Gka = new TGraph(Num,Beta,a0Karsch);
      TGraph *Gaa = new TGraph(Num,Beta,a0Amit);
      TGraph *GaaNEW = new TGraph(Num,Beta,a0AmitNEW);
      Custom(Gaa,2);Custom(GaaNEW,4);Custom(Gka,1);
      Gaa->Draw("ALP");GaaNEW->Draw("LP same"); Gka->Draw("LP same");
  
}

TGraphErrors * Subtract(TGraphErrors *g1,TGraphErrors *g2)
{
  const int N1=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2=g2->GetN(); double *X2=g2->GetX(), *Y2=g2->GetY(), *Ex2=g2->GetEX(), *Ey2=g2->GetEY();
  int N3; double X3[N1], Y3[N1], Ex3[N1], Ey3[N1];

  for(int i=0;i<N1;i++)
    {
      cout<<"x1[i] ="<<X1[i]<<"\t x2[i] = "<<X2[i];
      X3[i]=X1[i];   N3=N1;cout<<"X3[i]="<<X3[i]<<endl;
      double Factor = 1.0;
      Y3[i]=Factor*(Y1[i] - Y2[i]);//cout<<"Y3[i]="<<Y3[i]<<"\n"<<endl;
      Ex3[i]=0.0; Ey3[i]=Factor*TMath::Sqrt(pow(Ey1[i],2.0) + pow(Ey2[i],2.0));
    }
  TGraphErrors *g3 = new TGraphErrors(N3, X3,Y3,Ex3, Ey3);
  return g3;
}

TGraphErrors * ConvertToT(TGraphErrors *g1, int Nt, double Norm)
{
  const int N1=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2; double X2[N1], Y2[N1], Ex2[N1], Ey2[N1];

  for(int i=0;i<N1;i++)
    {
      //double a0 = LatticeSpacing(X1[i])/Norm;//cout<<"g0^2*a0^4="<<g0Square*pow(a0,4.0)<<endl;
      //double a0 = KarschLatticeSpacing(X1[i])/Norm;
      double a0 = AmitLatticeSpacingNEW(X1[i])/Norm;
      X2[i]=1.0/(Nt*a0);   N2=N1;cout<<"X2[i]="<<X2[i];
      double Factor = 1.0;
      Y2[i]=Factor*Y1[i];cout<<"Y2[i]="<<Y2[i]<<endl;                                                                                                                                               
      Ex2[i]=0.0; Ey2[i]=Factor*TMath::Sqrt(pow(Ey1[i],2.0));
    }
  TGraphErrors *g2 = new TGraphErrors(N2, X2,Y2,Ex2, Ey2);
  return g2;
}

double AsymptoticBetaFunction(double beta)
{
  double g0 = sqrt(6.0/beta);
  double Beta0 = 11.0/(16*pow(3.1415926,2.0));
  double Beta1 =102/pow(16.0*3.14159*3.14159 , 2.0);
  double B = 11.0*g0*g0/(24.0*3.141*3.141);
  double Result = pow(Beta0*g0*g0 , (-Beta1/(2*Beta0*Beta0)) )*exp(-1.0/(2*Beta0*g0*g0) );
  //cout<<"beta = "<<beta<<"\t means perturbative a*LambdaL ="<<Result<<" \n"<<endl;
  return Result;
}

double KarshRBetaFunction(double beta)
{
  double Beta[31]={5.7, 5.75, 5.8,  5.85,  5.9, 5.95,   6,      6.05,   6.1,    6.15,   6.2,  6.25, 6.3, 6.35, 6.4, 6.45,                    6.5,    6.55,   6.6,    6.65,   6.7,    6.75,   6.8,   6.85,  6.9,  6.95, 7,  7.05, 7.1, 7.15, 7.2};
  double BetaFunctionLattice[31]={0.46713,     0.475614,        0.490164,       0.510318,    0.5349,   0.562194, 0.589032,                                     0.614058,    0.63579,        0.654882,       0.672762,       0.689244,       0.704238,     
                                 0.718008,       0.730992,       0.742794,       0.753006,       0.76209,        0.770652, 
                                 0.7785,    0.785448,       0.791328,       0.795984,       0.799782,       0.80313,                                          0.806088,     0.808752,       0.811176,       0.813444,       0.815526,       0.817404};
  TGraph *g = new TGraph(31,Beta, BetaFunctionLattice);
  return g->Eval(beta);
}
double KarschLatticeSpacing(double beta)
{
  //double LambdaLattice=243/83.5; //Lambda(MSbar)=243 MeV , Karsch Lambda0=Tc/34.38 MeV
  double Beta[31]={5.7,	5.75, 5.8,  5.85,  5.9,	5.95,	6,	6.05,	6.1,	6.15,	6.2,  6.25, 6.3, 6.35, 6.4, 6.45,                       6.5,    6.55,   6.6,    6.65,   6.7,    6.75,   6.8,   6.85,  6.9,  6.95, 7,  7.05, 7.1, 7.15, 7.2};
  double lambdaFit[31]={2.1947,	2.087578,	1.99072,	1.905342,	1.831481,	1.768642,	1.715383,	                 1.670045,	1.631051,	1.596995,	1.567082,	1.540692,	1.517309,    1.496477,	              1.477892,	1.461275,	1.446329,	1.432777,	1.420468,	1.409274,	1.399057,	                   1.389672,	1.380963,	1.372792,	1.36508,	1.357773,	1.350826,	1.344204,	                1.337881,	1.331835,	1.326041};																														  
  double BetaFunctionLattice[31]={0.46713,     0.475614,	0.490164,	0.510318,    0.5349,   0.562194, 0.589032,	                           0.614058,	0.63579,	0.654882,	0.672762,	0.689244,	0.704238,                                       0.718008,	0.730992,	0.742794,	0.753006,	0.76209,	0.770652,	                             0.7785,	0.785448,	0.791328,	0.795984,	0.799782,	0.80313,	                                  0.806088,	0.808752,	0.811176,	0.813444,	0.815526,	0.817404};
  TGraph *g = new TGraph(31,Beta,lambdaFit);
  double Tc=265;                  //in MeV Critical temperature 265MeV, Pure SU(3)
  double LambdaLattice=Tc/34.38;  //in MeV By karsch paper Pure SU(3) scale setting
  double a0 =  AsymptoticBetaFunction(beta)*(g->Eval(beta))/LambdaLattice;
  return a0;
}

double AmitLatticeSpacing(double beta, int Nt)
{
  int FileCount=0; double *Beta, *Temperature; string *BetaS;
  int QCDType = 0;
if(QCDType<=0 && Nt==4)
    {
      double BETA[10] = {5.35,     5.6,      5.7,      5.8,      5.9,      6,        6.2,      6.35,     6.5,      6.6};
      string BETAS[10]= {"5.3500", "5.6000", "5.7000", "5.8000", "5.9000", "6.0000", "6.2000", "6.3500", "6.5000", "6.6000"};
      double T[10]   =  {191.423,  253.235,  283.295,  316.964,  354.68,   396.931,  497.306,  589.089,  697.975,  781.629};
      //Karsch                                                                                              
      //double T[10]   =  {137.753,  220.99,   270.21,   333.316,  405.588,  485.75,   667.11,   824.04,   1010.12,  1151.89};
      FileCount=10; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }

  if(QCDType<=0 && Nt==6)
    {
      double BETA[10] ={5.6,      5.85,     5.9,      6,        6.1,      6.25,     6.45,     6.6,      6.75,     6.85};
      string BETAS[10]={"5.6000", "5.8500", "5.9000", "6.0000", "6.1000", "6.2500", "6.4500", "6.6000", "6.7500", "6.8500"};
      double T[10]   ={197.84,    261.943,  277.094,  310.103,  347.084,  411.075,  515.307,  610.648,  723.786,  810.725};
      FileCount=10; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }

  if(QCDType<=0 && Nt==8)
    {
      double BETA[10] = {5.7,       5.95,      6,       6.1,      6.2,       6.35,     6.55,     6.7,      6.85,     6.95};
      string BETAS[10]= {"5.7000", "5.9500",  "6.0000","6.1000", "6.2000",  "6.3500", "6.5500", "6.7000", "6.8500", "6.9500"};
      double T[10]   = {194.765,  257.954,  272.89,  305.434,  341.898,   404.999,  507.794,     601.833,  713.438,  799.206};

      //Karsch                                                                                                      
      //double T[10]   = {135.10,   221.80,   242.87,   284.31,    333.55,    412.02,    539.30,   657.06,   801.33,   915.22};
      FileCount=10; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }

  TGraph *g = new TGraph(FileCount,Beta,Temperature);
  double a0 = 1.0/(Nt*g->Eval(beta));
  return a0;
  
}

double AmitLatticeSpacingNEW(double beta)
{ const int N=11;
  double NT[N]={4,5,   6,8,   10,12,  14,16,   18,20, 22};
  double BetaCritical[N]={5.69, 5.8,   5.89, 6.06,  6.2, 6.33,  6.44, 6.54,  6.63, 6.71,  6.79}; //5.2
  double LambdaCritical[N]={0,0,  0,0, 0,0, 0,0, 0,0, 0 };
  double Tc=265;                  //in MeV Critical temperature 265MeV, Pure SU(3)                                           
  double LambdaLattice=Tc/34.38;  //in MeV By karsch paper Pure SU(3) scale setting Tc/Lambda=34.38
  
  for(int i=0; i<N; i++)
    {LambdaCritical[i]=1.0/(NT[i]*AsymptoticBetaFunction(BetaCritical[i])*Tc/LambdaLattice);
      cout<<"Beta="<<BetaCritical[i]<<", LambdaCritical="<<LambdaCritical[i]<<endl;
    }
  TCanvas *can=new TCanvas("can","can",600,600);
  TH1D *h = new TH1D("h","h",100,5,7);
  TAxis *Ax=h->GetXaxis();
  for(int i=0; i<N; i++)
    {
      int BinNum = Ax->FindBin(BetaCritical[i]);
      h->SetBinContent(BinNum, LambdaCritical[i]);
    }
  TF1 *f = new TF1("f","([0]  +[1]*x + [2]*x*x)/([3]  + [4]*x + [5]*x*x)",BetaCritical[0]-0.2, BetaCritical[N-1]+0.2);
  f->SetParameter(0,1.0); f->FixParameter(1,0.0); f->SetParameter(2,1.0) ;f->SetParameter(3,1.0);f->FixParameter(4,0.0); f->SetParameter(5,1.0);
  h->Fit("f","R");
  double FitParameter0 = f->GetParameter(0), FitParameter1=f->GetParameter(1), FitParameter2=f->GetParameter(2), FitParameter3=f->GetParameter(3), FitParameter4=f->GetParameter(4), FitParameter5=f->GetParameter(5);
  //FitParameter0 = 1.0; FitParameter1=0; FitParameter2=0;
  //y = 0.3944x4 - 10.679x3 + 108.49x2 - 490.42x + 833.87

  //double a0 =  AsymptoticBetaFunction(beta)*( 0.3944*pow(beta,4) - 10.679*pow(beta,3) + 108.49*pow(beta,2) - 490.42*beta + 833.87)/LambdaLattice;
  double a0 =  AsymptoticBetaFunction(beta)*( (FitParameter2*pow(beta,2) + FitParameter1*pow(beta,1) + FitParameter0)/(FitParameter5*pow(beta,2) + FitParameter4*pow(beta,1) + FitParameter3 ) )/LambdaLattice;
  return a0;
}
