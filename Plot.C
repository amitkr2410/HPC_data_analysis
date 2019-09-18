#include"Custom.h"

// define a function of calculating the ratio between two root graphs

TH1D * GetRatioAsHist(string Name, TGraphErrors* Numerator, TGraphErrors* Denominator)
{
  int N1=Numerator->GetN(); double *X1=Numerator->GetX();double *Y1=Numerator->GetY();double *EY1=Numerator->GetEY(); double *EX1=Numerator->GetEX();
  int N2=Denominator->GetN(); double *X2=Denominator->GetX();double *Y2=Denominator->GetY();double *EY2=Denominator->GetEY(); double *EX2=Denominator->GetEX();
  double Ratio=0.0, RatioError=0.0;
  double * XBinArray = new double [N1+1]; // boundaries of bins of the new histogram
  for(int i=0;i<N1;i++)
    { XBinArray[i] = X1[i] - EX1[i];
    }
  XBinArray[N1] = X1[N1-1] + EX1[N1-1];
  TH1D *G = new TH1D(Name.c_str(), Name.c_str(), N1, XBinArray);
  if(N1==N2)
    {   for(int i=0;i<N1;i++)
	{ //if(Y1[i] > pow(10,-15) && Y2[i] > pow(10,-15) && fabs(Y1[i]) < pow(10,5) && fabs(Y2[i]) < pow(10,5))
	  {cout<<"Amit X1,X2 = "<<X1[i]<<","<<X2[i]<<endl;
	      Ratio = Y1[i]/Y2[i];cout<<"Amit Y1,Y2 = "<<Y1[i]<<","<<Y2[i]<<endl;
	      RatioError = Ratio*TMath::Sqrt( pow(EY1[i]/Y1[i],2.0) + pow(EY2[i]/Y2[i],2.0) );
	    }
	    //else
	    {//Ratio = 0.0; RatioError=0.0;
	    }
	      G->SetBinContent(i+1,Ratio);
	      G->SetBinError(i+1,RatioError);
	      cout<<"XBinNum = "<<i+1<<"\t"<<"X[i] ="<<X1[i]<<"\t"<<Ratio<<endl;
	}
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
  cout<<"\n\n";
  return G;
}

void Plot(int pTHatMin, int pTHatMax)
{
  CustomGlobal();
  char FileName[5000]; double RMax=pow(10,0.4), RMin=pow(10,-2);
    //int pTHatMin=110;
    //int pTHatMax=120;
  TFile *Spectra1  = new TFile("Files/Jetscape_Spectrum_7000GeV.root");
  TFile *Spectra2  = new TFile("Files/Pythia_Spectrum_7000GeV.root");
  char Name[1000];
  //sprintf(Name, "DifferentialJetCrossSectionBin%i_%i",pTHatMin,pTHatMax);
  //sprintf(Name, "DifferentialSingleHadronYieldBin%i_%i",pTHatMin,pTHatMax);
  sprintf(Name, "TotalDifferentialSingleHadronYield");
  //  sprintf(Name,"");
  TGraphErrors *R1         = (TGraphErrors *)Spectra1->FindObjectAny(Name)->Clone();
  TGraphErrors *g2        = (TGraphErrors *)Spectra2->FindObjectAny(Name)->Clone();
  

  TH1D *R         = GetRatioAsHist("Ratio", R1, g2);
  
  TCanvas *canvas = new TCanvas("canvas","canvas",600,450);
  //R->Draw("AL ");Custom(R,2);
  R->Draw("AXIS"); Custom(R,2);
  R->Draw("Hist E X0 ][  same");R->Draw("Hist ][ same");
  R->GetXaxis()->SetRangeUser(5,120);R->GetYaxis()->SetRangeUser(RMin, RMax);
  R->GetXaxis()->SetTitle(" parton #it{p}_{T}  (GeV)");R->GetYaxis()->SetTitle("diff cross section");
  gPad->SetMargin(0.22,0.03,0.2,0.1);
  //g1->Draw("L same");Custom(g1,2);
  //g2->Draw("L same");Custom(g2,1);
  //R->GetYaxis()->SetTitleOffset(1.2);
  //canvas->SetLogy();
  TLegend *leg1 = new TLegend(0.7,0.65,0.92,0.8);
  leg1->SetTextFont(22);leg1->SetTextSize(0.04);
  leg1->AddEntry(R, "JETSCAPE/PYTHIA","l");
  //leg1->AddEntry(R, "JETSCAPE","l");
  //leg1->AddEntry(g2, "PYHTIA","l");
  leg1->SetFillStyle(0);  leg1->Draw();
  TPaveText *Pt1 = new TPaveText(20,RMax*1.01,100,RMax*2.2,"NB");
  //sprintf(Name, "b-quark only, pTHat=(%i,%i)GeV, |#eta|<3.0",pTHatMin, pTHatMax);
  sprintf(Name, "Full spectrum, b-quark only,|#eta|<3.0");
  //sprintf(Name, "All particles, R=0.7, |#eta_{jet}|<2.0");  
  Pt1->AddText(Name);Pt1->SetTextAlign(12);Pt1->SetFillColor(0);
  Pt1->SetTextFont(22);Pt1->SetTextSize(0.04);Pt1->Draw();
  //sprintf(Name, "BQuark/Spectrum_b-quark_pTHat%i_%i.pdf",pTHatMin, pTHatMax);
  //sprintf(Name, "BQuark/FullSpectrum_b-quark.pdf");
  canvas->SaveAs(Name);

}
