//Fmunu in data files is (Umunu - Udaggermunu)/(2i)
//Here, We divide Fmunu by barecoupling constant g0=sqrt(6/Beta), and lattice spacing square a^2; See function Subtract()
//Thus, F^2_{munu} should be divided by g0^2 and a^4;
//We use perturbative two-loop function result to set the scale.
//Plot is for Correlator w.r.t Temperature.   
#include"iostream"
#include"Custom.h"
class EntropyPlot
{
public:
  TGraphErrors * Scale(TGraphErrors *g1);

};

void Plot()
{
  EntropyPlot LC;

  char FileName[3000];
  TGraphErrors *g1 = new TGraphErrors("EntropyDensity_RawData_SU3.txt","%lf %lf ");
  TGraphErrors *G2 = new TGraphErrors("EntropyDensity_RawData_QCD.txt","%lf %lf f");
  TGraphErrors *G1;
  G1 = LC.Scale(g1);

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(5);

  TCanvas *c1 = new TCanvas("c1","c1",420,400);
  gPad->SetMargin(0.14, 0.05, 0.13, 0.03);
  CustomGlobal();//1=Black, 2=Red, 4=light Blue, 9=DarkBlue, 6=Pink, 8=Green                                                                   
  Custom(G1,9);Custom(G2,2);
  G1->Draw("AL");
  G2->Draw("L");
  G1->SetLineWidth(7);G2->SetLineWidth(7);
  G1->GetYaxis()->SetTitle("#it{s/T^{3}}");G1->GetXaxis()->SetTitle("#it{T} (MeV)");
  G1->GetYaxis()->SetRangeUser(0, 25);
  G1->GetXaxis()->SetLimits(100, 1000);
  
  TLegend *L1 = new TLegend(0.2,0.68,0.6,0.96);
  L1->AddEntry(G2,"2+1 flavour","lp");
  L1->AddEntry(G1,"Pure gauge SU(3)","lp");
  L1->SetFillStyle(0);L1->SetTextFont(22);L1->SetTextSize(0.055);
  L1->Draw();
  
  char PDFFileOutput[800000],COMMAND[800000];
  sprintf(PDFFileOutput,"~/Dropbox/QhatLattice/Plots/EntropyDensityOverT3VsT.pdf");
  c1->SaveAs("MyFile1.eps");
  system("epstopdf MyFile1.eps");
  sprintf(COMMAND,"mv MyFile1.pdf %s",PDFFileOutput);
  cout<<string(COMMAND)<<endl;

  system(COMMAND);
}



TGraphErrors * EntropyPlot::Scale(TGraphErrors *g1)
{

  const int N1=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N3; double X3[N1], Y3[N1], Ex3[N1], Ey3[N1];
  for (int i=0; i<N1; i++)
    {
      X3[i] = X1[i]*265;
      N3=N1;
      cout<<"X1[i]="<<X1[i]<<"\tX3[i]="<<X3[i]<<endl;
      cout<<"Y1[i]="<<Y1[i]<<endl;
      Y3[i] = Y1[i]*4.0/3.0;
      Ey3[i] =0.0;
      Ex3[i] =0.0;
    }
  TGraphErrors *g3 = new TGraphErrors(N3,X3,Y3,Ex3,Ey3);
  return g3;
 
}
