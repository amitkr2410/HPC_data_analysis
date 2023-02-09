#include "StatTool.h"
#include "Custom.h"
using namespace std;
void Mean(TGraph *g, double *Ymean, double *Yerr, int FirstPoint, int SkipN);
void AnalyzeNtNs(int Nt, int Ns, string OutputFile);
void QhatRe(int Nt, int Ns, string File1Output, string File2Output);
void QhatIm(int Nt, int Ns, string File1Output, string File2Output);
double LatticeSpacing(double beta);
void StatisticalAverage()
{
  int Nt=4;
  int Ns=4*Nt;
  
  char FileNameOutput1[10000];sprintf(FileNameOutput1,"Test/MyFileNt%d_Ns%d.txt",Nt,Ns);
  char FileNameOutput2[10000];sprintf(FileNameOutput2,"Test/MyFileNt%d_Ns%d.txt",Ns,Ns);
  AnalyzeNtNs(Nt, Ns, FileNameOutput1);
  AnalyzeNtNs(Ns, Ns, FileNameOutput2);
  QhatRe(Nt,Ns, FileNameOutput1, FileNameOutput2);
  QhatIm(Nt, Ns, FileNameOutput1, FileNameOutput2);  

}

void AnalyzeNtNs(int Nt, int Ns, string OutputFile)
{
  CustomGlobal();
  string File[100];
  //  int Nt=4, Ns=16;
  string NtNs; 
  if(Nt==4){NtNs="Nt4_Ns16_";}
  if(Nt==16 && Ns==16){NtNs="Nt16_Ns16_";}
  int FileCount=7;  
  double Beta[7]={6.05, 6.125,  6.354, 6.608,  6.800,  7.150,  7.373};
  int FirstRun=200; int SkipN=3;
  int i=0;
  File[i]="DataTraceFmunu" + NtNs + "Beta6.0500_ml0.005320_ms0.106400_u0_1.000.txt"; i=i+1;
  File[i]="DataTraceFmunu" + NtNs + "Beta6.1250_ml0.004830_ms0.096600_u0_1.000.txt"; i=i+1;
  File[i]="DataTraceFmunu" + NtNs + "Beta6.3540_ml0.003640_ms0.072800_u0_1.000.txt"; i=i+1;
  File[i]="DataTraceFmunu" + NtNs + "Beta6.6080_ml0.002710_ms0.054200_u0_1.000.txt"; i=i+1;
  File[i]="DataTraceFmunu" + NtNs + "Beta6.8000_ml0.002240_ms0.044800_u0_1.000.txt"; i=i+1;
  File[i]="DataTraceFmunu" + NtNs + "Beta7.1500_ml0.001600_ms0.032000_u0_1.000.txt"; i=i+1;
  File[i]="DataTraceFmunu" + NtNs + "Beta7.3730_ml0.001250_ms0.025000_u0_1.000.txt"; i=i+1;
  
  
  ofstream fout; char FileNameOutput[10000];sprintf(FileNameOutput,"%s",OutputFile.c_str());
  fout.open(FileNameOutput, ios::out);
  TGraph *G;
  vector<double> Ymean , Yerr;

 for(int i=0; i<FileCount; i++)   
   { char FileName[10000]; Ymean.resize(0); Yerr.resize(0);
     double mean, err; mean=0; err=0;
     if(Nt!=0){sprintf(FileName,"RawDataFullQCDWrong/Nt%d_Ns%d/%s",Nt,Ns,File[i].c_str());}
     //if(Nt==Ns){sprintf(FileName,"/wsu/home/fy/fy41/fy4125/Lattice/MILC/MyFFMeasurement/Output/%s",File[i].c_str());}

     cout<<"\n F3iF3i.Re (Nt,Ns)="<<Nt<<","<<Ns<<endl;
     G = new TGraph(FileName,"%lf %lf");
     Mean(G, &mean, &err, FirstRun, SkipN);
     Ymean.push_back(mean); Yerr.push_back(err);

     cout<<"\n F3iF3i.Im (Nt,Ns)="<<Nt<<","<<Ns<<endl;
     G = new TGraph(FileName,"%lf %*lf %lf");
     Mean(G, &mean, &err, FirstRun, SkipN);
     Ymean.push_back(mean); Yerr.push_back(err);

     cout<<"\n F3iF4i.Re (Nt,Ns)="<<Nt<<","<<Ns<<endl;
     G = new TGraph(FileName,"%lf %*lf %*lf     %*lf %*lf     %lf");
     Mean(G, &mean, &err, FirstRun, SkipN);
     Ymean.push_back(mean); Yerr.push_back(err);

     cout<<"\n F3iF4i.Im (Nt,Ns)="<<Nt<<","<<Ns<<endl;     
     G = new TGraph(FileName,"%lf %*lf %*lf     %*lf %*lf     %*lf %lf ");
     Mean(G, &mean, &err, FirstRun, SkipN);
     Ymean.push_back(mean); Yerr.push_back(err);
     
     fout<<Beta[i]<<"\t"<<Ymean[0]<<"\t"<<Yerr[0]<<"\t"<<Ymean[1]<<"\t"<<Yerr[1]<<"\t"
                        <<Ymean[2]<<"\t"<<Yerr[2]<<"\t"<<Ymean[3]<<"\t"<<Yerr[3]<<endl;
   }

 TCanvas *c1 = new TCanvas(NtNs.c_str(),NtNs.c_str(),900,400);
 c1->Divide(2,1);
 c1->cd(1);
 TGraphErrors *g = new TGraphErrors(FileNameOutput,"%lf %lf %lf");
 TGraphErrors *g2 = new TGraphErrors(FileNameOutput,"%lf %*lf %*lf %lf %lf");
 Custom(g,2); Custom(g2,1);
 g->Draw("AL *");g2->Draw("L* same ");
 //g->GetYaxis()->SetRangeUser(0.0,0.034);
 g->GetYaxis()->SetTitle("g^{2}a^{4}< F^{3i}F^{3i} - F^{4i}F^{4i} >");
 g->GetXaxis()->SetTitle("#beta=10/g^{2}");
 gPad->SetMargin(0.2,0.02,0.2,0.05);
 g->GetYaxis()->SetTitleOffset(1.6);

 c1->cd(2);
 TGraphErrors *G1 = new TGraphErrors(FileNameOutput,"%lf   %*lf %*lf %*lf %*lf %lf %lf ");
 TGraphErrors *G2 = new TGraphErrors(FileNameOutput,"%lf   %*lf %*lf %*lf %*lf %*lf %*lf %lf %lf");
 Custom(G1,2); Custom(G2,1);
 G1->Draw("AL* ");G2->Draw(" L* same");
 //G1->GetYaxis()->SetRangeUser(0.1,4.0);
 G1->GetYaxis()->SetTitle("2g^{2}a^{4}< F^{3i}F^{4i} >");
 G1->GetXaxis()->SetTitle("#beta=10/g^{2}"); 
 gPad->SetMargin(0.2,0.02,0.2,0.05);G1->GetYaxis()->SetTitleOffset(1.8);
  fout.close();

}

void QhatRe(int Nt, int Ns,string File1Output, string File2Output)
{
  double PreFactor, Qminus=20000, MuSquare=20000.0,alphas=0.2, NC=3;
  double TOnePlusTtwo = 2.0*16*2.0;
  PreFactor=(8*sqrt(2)*pow(3.1415926,1.0)*alphas/(NC*TOnePlusTtwo));

  char FileName[50000];
  sprintf(FileName,"%s",File1Output.c_str());
  TGraphErrors *g1 = new TGraphErrors(FileName,"%lf %lf %lf");
  sprintf(FileName,"%s",File2Output.c_str());
  TGraphErrors *g2 = new TGraphErrors(FileName,"%lf %lf %lf");
  const int N=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2=g2->GetN(); double *X2=g2->GetX(), *Y2=g2->GetY(), *Ex2=g2->GetEX(), *Ey2=g2->GetEY();
  double X[N], Y[N], Ex[N], Ey[N]; double YY[N], Eyy[N];
  
  for(int i=0;i<N;i++)
    {
      double Beta=X1[i];
      double g0Square=10.0/Beta;
      double T = LatticeSpacing(Beta);
      double a0= 1.0/(Nt*T);
      double Factor = 1/((4.0*g0Square*pow(a0,4.0)));
      double TDepPart = pow(T,-1.0);
      Y[i]=PreFactor*TDepPart*( (Y1[i] - Y2[i])*Factor)/pow(T,3.0) ;
      Ey[i]=PreFactor*TDepPart*TMath::Sqrt(pow(Ey1[i],2.0) + pow(Ey2[i],2.0))*Factor/pow(T,3.0);
      Ex[i]=0;
      X[i] = T;
      cout<<"(Beta1,Beta2)=("<<X1[i]<<","<<X2[i]<<") \t (Y1,Y2)=("<<Y1[i]<<","<<Y2[i]<<")"<<endl;
      cout<<"Y1-Y2 = "<<Y1[i]-Y2[i]<<endl;
      cout<<"T="<<T<<"MeV\t"<<Y[i]<<endl; YY[i]=Y[i]/2.0; Eyy[i]=Ey[i]/2.0;
    }
  cout<<endl;
  TCanvas *cc1 = new TCanvas("QhatReal","QhatReal",600,500);
  CustomGlobal();
  TGraphErrors *G = new TGraphErrors(N,X,Y,Ex,Ey);
  TGraphErrors *GG = new TGraphErrors(N,X,YY,Ex,Eyy);
  Custom(G,2);Custom(GG,2) ;
  G->Draw("AL* ");//GG->Draw("L same");
  //G->GetYaxis()->SetRangeUser(-0.1,6.0);
  G->GetYaxis()->SetTitle("Real(#hat{q}/T^{3})");
  G->GetXaxis()->SetTitle("T(MeV)");
  gPad->SetMargin(0.12,0.02,0.2,0.05);G->GetYaxis()->SetTitleOffset(0.8);

}

void QhatIm(int Nt, int Ns,string File1Output, string File2Output)
{
  double PreFactor, Qminus=20000, MuSquare=20000.0,alphas=0.2, NC=3;
  double TOnePlusTtwo = 2.0*16*2.0;
  PreFactor=(8*sqrt(2)*pow(3.1415926,1.0)*alphas/(NC*TOnePlusTtwo));

  char FileName[50000];
  sprintf(FileName,"%s",File1Output.c_str());
  TGraphErrors *g1 = new TGraphErrors(FileName,"%lf %*lf %*lf %*lf %*lf %lf %lf");
  sprintf(FileName,"%s",File2Output.c_str());
  TGraphErrors *g2 = new TGraphErrors(FileName,"%lf %*lf %*lf %*lf %*lf %lf %lf");
  const int N=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2=g2->GetN(); double *X2=g2->GetX(), *Y2=g2->GetY(), *Ex2=g2->GetEX(), *Ey2=g2->GetEY();
  double X[N], Y[N], Ex[N], Ey[N]; double YY[N], Eyy[N];
  
  for(int i=0;i<N;i++)
    {
      double Beta=X1[i];
      double g0Square=10.0/Beta;
      double T = LatticeSpacing(Beta);
      double a0= 1.0/(Nt*T);
      double Factor = 1/((4.0*g0Square*pow(a0,4.0)));
      double TDepPart = pow(T,-1.0);
      Y[i]=PreFactor*TDepPart*( (Y1[i] - Y2[i])*Factor)/pow(T,3.0) ;
      Ey[i]=PreFactor*TDepPart*TMath::Sqrt(pow(Ey1[i],2.0) + pow(Ey2[i],2.0))*Factor/pow(T,3.0);
      Ex[i]=0;
      X[i] = T;
      cout<<"(Beta1,Beta2)=("<<X1[i]<<","<<X2[i]<<") \t (Y1,Y2)=("<<Y1[i]<<","<<Y2[i]<<")"<<endl;
      cout<<"Y1-Y2 = "<<Y1[i]-Y2[i]<<endl;
      cout<<"T="<<T<<"MeV\t"<<Y[i]<<endl; YY[i]=Y[i]/2.0; Eyy[i]=Ey[i]/2.0;
    }

  TCanvas *cc2 = new TCanvas("QhatIm","QhatIm",600,500);
  CustomGlobal();
  TGraphErrors *G = new TGraphErrors(N,X,Y,Ex,Ey);
  TGraphErrors *GG = new TGraphErrors(N,X,YY,Ex,Eyy);
  Custom(G,2);Custom(GG,1) ;
  G->Draw("AL* ");//GG->Draw("L same");
  //G->GetYaxis()->SetRangeUser(-0.1,6.0);
  G->GetYaxis()->SetTitle("Im(#hat{q}/T^{3})");
  G->GetXaxis()->SetTitle("T(MeV)");
  gPad->SetMargin(0.12,0.02,0.2,0.05);G->GetYaxis()->SetTitleOffset(0.8);

}

double LatticeSpacing(double beta)
{
  double Beta[7]={6.05, 6.125,  6.354, 6.608,  6.800,  7.150,  7.373};
  double    T[7]={232,  249,    311,   399,    480,    669,    819};
  TGraph *g = new TGraph(7,Beta,T);
  return g->Eval(beta);
}
