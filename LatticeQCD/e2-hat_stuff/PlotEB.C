#include "../Script/Custom.h"
#include "../Script/StatTool.h"
int Color(int Nt);
double RenormalizationConstant(double g0Square, string QCDTypeS);
TGraph *ScaleByZ3(TGraph *g);
TGraph *ScaleByRBeta(TGraph *g);
void MakePlot(int nt, int ns,int  QCDType,TCanvas *c,TCanvas *c2,TCanvas *c3,TCanvas *c4);
double AsymptoticBetaFunction(double beta);
double LatticeSpacing(int Nt, int QCDType, double beta);
TGraph *ConvertToObservable(TGraph *g, int Nt, int QCDType);
void PlotEB()
{
  int Nt[4]={4,6,8,10}, Ns[4]={16,24,32,40},  QCDType=0;
   TCanvas *c  = new TCanvas("c","c", 500, 600);   
   TCanvas *c2 = new TCanvas("c2","c2", 500, 500);
   TCanvas *c3 = new TCanvas("c3","c3", 500, 500);
   TCanvas *c4 = new TCanvas("c4","c4", 500, 500);

   for(int i=0;i<4;i++)
     {
       MakePlot(Nt[i],Ns[i], QCDType,c,c2,c3,c4);
     }
   
}
void MakePlot(int Nt, int Ns, int  QCDType, TCanvas *c,TCanvas *c2,TCanvas *c3,TCanvas *c4)
{
  int NT=4;
  CustomGlobal();
  string QCDTypeS;
  if(QCDType<=0){QCDTypeS="PureGauge";}
  if(QCDType > 0){QCDTypeS="FullQCD";}

  char FileName1[1000], FileName2[1000]; 
  sprintf(FileName1,"DataEB_Nt%i_Ns%i.txt",Nt,Ns);	
  sprintf(FileName2,"DataEB_Nt%i_Ns%i.txt",Ns,Ns);		  
  TGraph *g = new TGraph(FileName1,"%lf %lf %*lf %*lf");
  TGraph *g2 = new TGraph(FileName1,"%lf %*lf %lf %*lf");
  TGraph *g3 = new TGraph(FileName1,"%lf %*lf %*lf %lf");

  TGraph *gv = new TGraph(FileName2,"%lf %lf %*lf %*lf");
  TGraph *gv2 = new TGraph(FileName2,"%lf %*lf %lf %*lf");
  TGraph *gv3 = new TGraph(FileName2,"%lf %*lf %*lf %lf");

  TGraph *G; G = SubtractTGraph(g,gv);
  TGraph *G2; G2 = SubtractTGraph(g2,gv2);
  TGraph *G3; G3 = SubtractTGraph(g3,gv3);
  TGraph *H;  H	 = SubtractTGraph(G,G2);
  TGraph *H2; H2 = AddTGraph(G,G2);

  TGraph *UnRenE2hatT; UnRenE2hatT = ConvertToObservable(G2, Nt, QCDType);
  string OutputFileS0="UnRenormalizedNt" + to_string(Nt) + "_" + QCDTypeS + ".txt";
  WriteToFile(UnRenE2hatT,OutputFileS0);
  
  TGraph *RMinus; RMinus = ScaleByZ3(H);
  TGraph *RPlus; RPlus =  ScaleByRBeta(H2);

  TGraph *F3iS; F3iS = AddTGraph(RPlus, RMinus); F3iS = ScaleGraph(F3iS,0.5);
  TGraph *F4iS; F4iS = SubtractTGraph(RPlus, RMinus); F4iS =	ScaleGraph(F4iS,0.5);
  TGraph *E2hat; E2hat = ScaleGraph(F4iS,-0.5);

  TGraph *RMinusT, *E2hatT;
  RMinusT = ConvertToObservable(RMinus, Nt, QCDType);
  E2hatT = ConvertToObservable(E2hat, Nt, QCDType);
  string OutputFileS="RenormalizedNt" + to_string(Nt) + "_" + QCDTypeS + ".txt";
  WriteToFile(E2hatT,OutputFileS);

  //TCanvas *c= new TCanvas("c","c",500,600);
  c->cd();
  Custom(g,0);Custom(g2,2);Custom(g3,4);
  Custom(gv,1);Custom(gv2,2);Custom(gv3,4);
  Custom(G,1);Custom(G2,2);Custom(G3,4);
  Custom(H,1);Custom(H2,2);
 if(Nt==NT){  H->Draw("ALp");}
  //gv->Draw("ALp");//  g2->Draw("Lp same");  g3->Draw("Lp same");
  //gv->Draw("Lp same"); gv2->Draw("Lp same");  gv3->Draw("Lp same");
  G->Draw("Lp same"); G2->Draw("Lp same");//  G3->Draw("Lp same");
  H->Draw("Lp same"); H2->Draw("Lp same");

  G->SetLineStyle(6); G2->SetLineStyle(8); G3->SetLineStyle(8);
  G->SetLineWidth(12); G2->SetLineWidth(10); G3->SetLineWidth(10);
  
  H->GetYaxis()->SetTitle("<O>");
  H->GetXaxis()->SetTitle("#beta_{0}=6/g^{2}_{0}");
  H->GetYaxis()->SetRangeUser(-0.02,0.015);
  gPad->SetMargin(0.15,0.01,0.15,0.01);


  //TCanvas *c2 = new TCanvas("c2","c2", 500, 500);
  c2->cd();
  Custom(H,1);Custom(H2,2); Custom(RMinus,1);   Custom(RPlus,2);
  Custom(F3iS, 4); Custom(F4iS,6); RMinus->SetLineStyle(8); RPlus->SetLineStyle(8);
  if(Nt==NT){H->Draw("ALp same");}
  H2->Draw("Lp same"); RMinus->Draw("Lp same"); RPlus->Draw("Lp same");
  F3iS->Draw("Lp same");F4iS->Draw("Lp same");
  H->GetYaxis()->SetTitle("<O>");
  H->GetXaxis()->SetTitle("#beta_{0}=6/g^{2}_{0}");
  H->GetYaxis()->SetRangeUser(-0.004,0.0025);
  gPad->SetMargin(0.15,0.01,0.15,0.01);

  //TCanvas *c3 = new TCanvas("c3","c3", 500, 500);
  c3->cd();
  Custom(RMinus,Color(Nt));Custom(E2hat,Color(Nt));  
  RMinus->SetLineStyle(8); RMinus->SetMarkerStyle(26);
  RMinus->SetLineWidth(6); E2hat->SetLineWidth(6);
  if(Nt==NT){RMinus->Draw("ALp same");}
  RMinus->Draw("Lp same"); E2hat->Draw("Lp same");
  RMinus->GetYaxis()->SetTitle("<O>");
  RMinus->GetXaxis()->SetTitle("#beta_{0}=6/g^{2}_{0}");
  RMinus->GetYaxis()->SetRangeUser(-0.0005,0.02);
  gPad->SetMargin(0.15,0.01,0.15,0.01);

  c4->cd();
  Custom(RMinusT,Color(Nt));Custom(E2hatT,Color(Nt));
  RMinusT->SetLineStyle(8); RMinusT->SetMarkerStyle(26);
  RMinusT->SetLineWidth(6); E2hatT->SetLineWidth(6);

  if(Nt==NT){RMinusT->Draw("ALp same");}
  RMinusT->Draw("Lp same"); E2hatT->Draw("Lp same");
  RMinusT->GetYaxis()->SetTitle("<O>");
  RMinusT->GetXaxis()->SetTitle("T (MeV)");
  RMinusT->GetYaxis()->SetRangeUser(-0.0005,4.5);
  gPad->SetMargin(0.15,0.01,0.15,0.01);
}
int Color(int Nt)
{ int Ans=1;
  if(Nt==4){Ans=2;}
  if(Nt==6){Ans=3;}
  if(Nt==8){Ans=1;}
   if(Nt==10){Ans=4;}
  return Ans;
}
TGraph *ScaleByZ3(TGraph *g)
{
    const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();
    double x[N], y[N], ex[N], ey[N]; double g0Square=0;
    for(int i=0;i<N;i++)
      {
	x[i]=X[i]; g0Square = 6.0/X[i];
	y[i]  = Y[i]*RenormalizationConstant(g0Square,"PureGauge");
      }
    TGraph *gr = new TGraphErrors(N,x,y);
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

TGraph *ScaleByRBeta(TGraph *g)
{
    const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();
    double x[N], y[N], ex[N], ey[N]; double g0Square=0;
    TGraph *RBeta = new TGraph("R_beta_PureSU3GaugeAction.txt","%lf %lf");
    for(int i=0;i<N;i++)
      {
        x[i]  = X[i];
        y[i]  = Y[i]*RBeta->Eval(X[i]);
      }
    TGraph *gr = new TGraphErrors(N,x,y);
    return gr;
}

TGraph *ConvertToObservable(TGraph *g, int Nt, int QCDType)
{
    const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();
    double x[N], y[N], ex[N], ey[N]; double Beta=0, g0Square=0, T=0;
    for(int i=0;i<N;i++)
      {
        Beta=X[i]; if(QCDType<=0){ g0Square = 6.0/X[i];}
	T=LatticeSpacing(Nt, QCDType, Beta);
	x[i] = T;
        y[i]  = Y[i]*pow(Nt,4)/g0Square;
      }
    TGraph *gr = new TGraphErrors(N,x,y);
    return gr;
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
  
double LatticeSpacing(int Nt, int QCDType, double beta)
{
  double *Beta;string *BetaS; int FileCount;
  double *Temperature;
if(QCDType<=0 && Nt==4)
    {
      double BETA[9] = {5.6,      5.7,      5.8,      5.9,      6,        6.2,      6.35,     6.5,      6.6};
      string BETAS[9]= {"5.6000", "5.7000", "5.8000", "5.9000", "6.0000", "6.2000", "6.3500", "6.5000", "6.6000"};
      double T[9] =    {208.61,   270.65,   335.85,    405.79,   481.89,   657.89,   816.19,   1002.76,  1145.77};
      //Old: 191.423                                                                          
      //double T[9]   =  {  253.235,  283.295,  316.964,  354.68,   396.931,  497.306,  589.089,  697.975,  781.629};  
      //Karsch:  137.753                                                                                     
      //double T[9]   =  {  220.99,   270.21,   333.316,  405.588,  485.75,   667.11,   824.04,   1010.12,  1151.89}; 
      FileCount=9; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }
if(QCDType<=0 && Nt==6)
    {
      double BETA[10] ={5.6,      5.85,     5.9,      6,        6.1,      6.25,     6.45,     6.6,      6.75,     6.85};
      string BETAS[10]={"5.6000", "5.8500", "5.9000", "6.0000", "6.1000", "6.2500", "6.4500", "6.6000", "6.7500", "6.8500"};
      double T[10] =   {139.075,  246.76,   270.532,   321.264,  376.988,  471.904,  624.72,  763.846,   928.821,  1055.68};
      //Old: double T[10]   ={197.84,    261.943,  277.094,  310.103,  347.084,  411.075,  515.307,  610.648,  723.786,  810.725};                                                                                                    
      FileCount=10; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }

if(QCDType<=0 && Nt==8)
    {
      double BETA[10] = {5.7,       5.95,      6,       6.1,      6.2,       6.35,     6.55,     6.7,      6.85,     6.95};
      string BETAS[10]= {"5.7000", "5.9500",  "6.0000","6.1000", "6.2000",  "6.3500", "6.5500", "6.7000", "6.8500", "6.9500"};
      double T[10]    = {135.326,   221.497,   240.948,  282.741, 328.949,   408.096,  536.123,  652.992,   791.763, 898.55};
      //Old: double T[10]   = {194.765,  257.954,  272.89,  305.434,  341.898,   404.999,  507.794,     601.833,  713.438, 799.206};                                                                                                          
      //Karsch                                                                                             
      //double T[10]   = {135.10,   221.80,   242.87,   284.31,    333.55,    412.02,    539.30,   657.06,   801.33,   915.22};                                                                                                    
      FileCount=10; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }
//TGraph *g = new TGraph(FileCount,Beta,Temperature);
// return g->Eval(beta);
 
 double Tc=265; //in MeV Critical temperature 265MeV, Pure SU(3)
 double LambdaLattice=Tc/34.38;  //in MeV By karsch paper Pure SU(3) scale setting Tc/Lambda=34.38
 double FitParameter0 = -23.2088, FitParameter1 = 0.993504, FitParameter2 = -25.2466, FitParameter3 = 0.905950;
 double a0 =  AsymptoticBetaFunction(beta)*( (FitParameter0 + FitParameter1*pow(beta,2))/( FitParameter2 + FitParameter3*pow(beta,2)) )/LambdaLattice;
 double AnsT=1.0/(Nt*a0); cout<<"Nt="<<Nt<<", beta0="<<beta<<", T="<<AnsT<<endl;  
 return AnsT;
 
}
