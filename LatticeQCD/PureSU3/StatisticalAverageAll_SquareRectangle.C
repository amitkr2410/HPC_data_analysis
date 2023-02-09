#include <stdlib.h>
#include "StatTool.h"
#include "Custom.h"

//void Mean(TGraph *g, double *Ymean, double *Yerr, int FirstPoint, int SkipN);
void AnalyzeNtNs(int Nt,int Ns,int QCDType, int WithOrWithoutTadpoleCorrection,  string File0, string File00, string File1,  string File11, string File20); 
TGraphErrors * Compute_AlphasFFOverT4_Vs_T(TGraphErrors *g,TGraphErrors *g2, int Nt, int QCDType, TGraphErrors *gU0, int FlagU0 );
TGraphErrors * Compute_1x1_AlphasFFOverT4_Vs_T(TGraphErrors *g,TGraphErrors *g2, int Nt, int QCDType);
double BareCouplingSquare(int QCDType, double Beta);
double StrongCouplingOneLoop(int QCDType, int Nt, double T);
double RbetaFmunuSquare(int QCDType, double beta);
double LatticeSpacing(int Nt, int QCDType, double beta);
double TadpoleFactorTUMQCD(double beta);
void VacuumResultNt32(double beta, int j, int WithOrWithoutTadpoleCorrection, double *MeanAndErrorFit);

void StatisticalAverageAll_SquareRectangle(int NT)
{
  //gSystem->Exec("rm Test");
  //gSystem->Exec("mkdir Test");
  int Nt=NT;
  int Ns=Nt*4;
  int QCDType=1, WithOrWithoutTadpoleCorrection=0;
  double Qminus[5]={100000, 20000, 50000, 80000, 100000}; //in MeV
  //Plaquette, tadpolefactor, and bare free energy
  char Output0[10000], Output00[10000];
  //   T+V ,            V (w.r.t beta)
  char Output1[10000], Output2[10000];
  char Output11[10000], Output20[10000], Output21[10000]; 
  string QCDTypeS, FmunuTypeS;
  if(QCDType<=0){QCDTypeS="PureGauge";}
  if(QCDType>0){QCDTypeS="FullQCD";}

  sprintf(Output0, "Test/%sDataPloopNt%d_Ns%d.txt",QCDTypeS.c_str(),Nt,Ns);
  sprintf(Output00,"Test/%sDataPloopNt%d_Ns%d.txt",QCDTypeS.c_str(),Ns,Ns);  
  sprintf(Output1, "Test/%sSquareRectangleVsBeta",QCDTypeS.c_str());
  //sprintf(Output2, "Test/%s%sOperatorVsBeta_NLO_Nt%d_Ns%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt,Ns);
  //sprintf(Output3, "Test/%s%sOperatorVsBeta_LO_Nt%d_Ns%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Ns,Ns);
  //sprintf(Output4, "Test/%s%sOperatorVsBeta_NLO_Nt%d_Ns%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Ns,Ns);

  sprintf(Output11,"Test/%s%sVacuumSubtractedOperatorVsBeta_RealPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt); //Vacuum subtracted Square and Rectangle  w.r.t beta
  sprintf(Output20,"Test/%sVacuumSubtracted_AlphasFFOverT4_Vs_T_RealPartNt%d.txt",QCDTypeS.c_str(), Nt);

  AnalyzeNtNs(Nt, Ns, QCDType, WithOrWithoutTadpoleCorrection,  Output0, Output00, Output1,  Output11, Output20);
  
  //DimensionlessOperatorVsT(Nt, Qminus, Output11, Output12,   Output21, Output22);


}

void AnalyzeNtNs(int Nt, int Ns, int QCDType, int WithOrWithoutTadpoleCorrection, string OutputFile0, string OutputFile00, string OutputFile1, string OutputFile11, string OutputFile20)
{
  CustomGlobal();
  int FirstRun=400; int SkipN=1; if(QCDType>0){SkipN=1; }
  int Norm;
  string  NtNs, QCDTypeS, QCDTypeS1, QCDTypeS2; 

  if(QCDType<=0){QCDTypeS="PureGauge";Norm=1.0;}
  if(QCDType>0){QCDTypeS="FullQCD";Norm=1.0;}

  int FileCount,NtTemp; double *Beta; string *BetaS; double *Temperature; char NtTempS[100], NsTempS[100];  

  if(QCDType<=0 && Nt==4)
    {
      double BETA[8] = {5.6,      5.7,      5.8,      5.9,      6,        6.2,      6.35,     6.5    };
      string BETAS[8]= {"5.6000", "5.7000", "5.8000", "5.9000", "6.0000", "6.2000", "6.3500", "6.5000"};
      double T[8] =    {208.61,   270.65,   335.85,    405.79,   481.89,   657.89,   816.19,   1002.76};
      //Old: 191.423 
      //double T[9]   =  { 253.235,  283.295,  316.964,  354.68,   396.931,  497.306,  589.089,  697.975,  781.629};
      //Karsch: 137.753,
      //double T[9]   =  220.99,   270.21,   333.316,  405.588,  485.75,   667.11,   824.04,   1010.12,  1151.89};
      FileCount=8; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
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
      double BETA[8] = {     6,       6.1,      6.2,       6.35,     6.55,     6.7,      6.85,     6.95};
      string BETAS[8]= {  "6.0000","6.1000", "6.2000",  "6.3500", "6.5500", "6.7000", "6.8500", "6.9500"};
      double T[8]    = {   240.948,  282.741, 328.949,   408.096,  536.123,  652.992,   791.763, 898.55};
      //Old: double T[10]   = {194.765,  257.954,  272.89,  305.434,  341.898,   404.999,  507.794,     601.833,  713.438,  799.206};
      //Karsch
      //double T[10]   = {135.10,   221.80,   242.87,   284.31,    333.55,    412.02,    539.30,   657.06,   801.33,   915.22};
      FileCount=8; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }

  if(QCDType>0 && Nt==4)
    {
      double BETA[8] ={5.9,      6.0,       6.285,     6.515,    6.664,    6.95,      7.15,     7.373};
      string BETAS[8]={"5.9000", "6.0000",  "6.2850",  "6.5150", "6.6640", "6.9500",  "7.1500", "7.3730"};
      double MS[8]  ={0.132,     0.1138,    0.079,     0.0603,   0.0514,   0.0386,    0.032,    0.025};
      double T[8]   ={201,    221,     291,     364,     421,     554,     669,     819};
      FileCount=8; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }


  if(QCDType>0 && Nt==6)
    {
      double BETA[10] ={6,        6.215,    6.285,    6.423,    6.664,    6.95,     7.15,      7.373,    7.596,       7.825};
      string BETAS[10]={"6.0000", "6.2150", "6.2850", "6.4230", "6.6640", "6.9500", "7.1500",  "7.3730", "7.5960",   "7.8250"};
      double MS[10]  ={0.1138,  0.0862,   0.079,  0.067,   0.0514,  0.0386,  0.032,   0.025,   0.0202,  0.0164};
      double T[10]   ={147,      181,      194,    222,     281,     370,     446,     547,     667,     815};
      FileCount=10; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }
  
  if(QCDType>0 && Nt==8)
    { // 6.2850, 7.8250
      double BETA[9] ={  6.515,   6.575,   6.664,   6.95,    7.28,    7.5,     7.596, 7.8250,  8.2};
      string BETAS[9]={ "6.5150", "6.5750",  "6.6640",  "6.9500", "7.2800",  "7.5000", "7.5960", "7.8250",  "8.2000"};
      //0.079, 0.0164
      double MS[9]  ={  0.0604,  0.0564,  0.0514,  0.0386,  0.0284,  0.0222,  0.0202, 0.0164,  0.01167};
      //146, 611
      double T[9]   ={ 182,     193,     211,     277,     377,     459,     500,  611,     843};
      FileCount=9; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }

  ofstream fout; char FileNameOutput[10000];
  string File[3],File2; char FileName[3][10000], FileName2[10000];  string LOorNLO[3], TadpoleOrder;
  char FileNameOutputLO[10000], FileNameOutputNLO[10000], FileNameOutputNNLO[10000], FileNameOutputNNNNLO[10000];
  char FileNameOutputVacLO[10000],FileNameOutputVacNLO[10000], FileNameOutputVacNNLO[10000], FileNameOutputVacNNNNLO[10000];
  TGraph *GM0, *GM1, *GM2;
  double Ymean[2], Yerr[2];
  double mean, err; double MeanAndErrorFit[2]={0,0};
  for(int k=0; k<2;k++)                   //Thermal or Vacuum
    {
      for(int i=0; i<FileCount; i++)      // different beta
	{ 
	  NtTemp=Nt;
	  if(k==1){NtTemp=Ns;}//NtNs= "Nt16_Ns16_"                                                                                          
	  sprintf(NtTempS,"%d",NtTemp);
	  sprintf(NsTempS,"%d",Ns);

	  for(int j=0; j<1;j++)           // Tadpole0 or Tadpole2 or Tadpole4 D4z derivatives terms
	    {
	      File[j]=std::string("DataTraceSquareRectangleLoop") + "Nt" + NtTempS + "_Ns"+ NsTempS  +"_" + "Beta"+ BetaS[i] + ".txt";
	      sprintf(FileName[j],"RawData%sSquareRectangle/Nt%d_Ns%d/%s",QCDTypeS.c_str(), NtTemp, Ns, File[j].c_str());
	      //RawDataFullQCDCloverTraceless/Nt4_Ns16/DataTraceFmunuTadpole0_Clover_Traceless_Nt4_Ns16_Beta5.9000.txt 
              File2=std::string("DataPloopNt") + NtTempS + "_Ns"+ NsTempS  + "_" + "Beta"+ BetaS[i] + ".txt";
              sprintf(FileName2,"RawData%s_Plaquette_PLoopRenormalizationU0/Nt%d_Ns%d/%s",QCDTypeS.c_str(), NtTemp, Ns, File2.c_str());
              //RawDataFullQCD_Plaquette_PLoop/Nt4_Ns16/DataPloopNt4_Ns16_Beta5.9000.txt
	      for(int m=0;m<2;m++){ Ymean[m]=0; Yerr[m]=0;}
	      mean=0; err=0;
	    }
	  double MEAN=0, ERROR=0;
	  cout<<"\n Computing Avg plaq Tr(U1234)/3.0, tadpole factor u0, and bare free energy=-log(|P|) for "<<NtTemp<<","<<Ns<<endl;
	  GM0 = new TGraph(FileName2,"%lf %lf");
	  JackknifeMean(GM0, &MEAN, &ERROR, FirstRun, SkipN);
	  double AveragePlaquette = MEAN;  double AveragePlaquetteError = ERROR;
	  double TadpoleFactorU0= pow(AveragePlaquette, 1.0/4.0); double TadpoleFactorU0Error=TadpoleFactorU0*ERROR/(4.0*MEAN);
	  GM0 = new TGraph(FileName2,"%lf %*lf %*lf %*lf %*lf   %*lf %*lf %*lf %*lf   %*lf %lf %*lf %*lf %*lf");
	  JackknifeMean(GM0, &MEAN, &ERROR, FirstRun, SkipN);
	  double BareFreeEnergy = MEAN; double BareFreeEnergyError= ERROR;
	  if(k==0){sprintf(FileNameOutput,"%s",OutputFile0.c_str());}  //Thermal                                                             
	  if(k==1){sprintf(FileNameOutput,"%s",OutputFile00.c_str());} //Vacuum                                            

	  if(QCDType>0 && Beta[i]>=6.488)
	    {TadpoleFactorU0= TadpoleFactorTUMQCD( Beta[i]); AveragePlaquette = pow(TadpoleFactorU0,4.0);        
	    }
	  ofstream fout0; fout0.open(FileNameOutput, ios::app);
	  fout0<<Beta[i]<<"\t"<<AveragePlaquette<<"\t"<<AveragePlaquetteError<<"\t"
	       <<TadpoleFactorU0<<"\t"<<TadpoleFactorU0Error<<"\t"
	       <<BareFreeEnergy<<"\t"<<BareFreeEnergyError<<endl;
	  cout<<Beta[i]<<"\t"<<AveragePlaquette<<"\t"<<AveragePlaquetteError<<"\t"
	      <<TadpoleFactorU0<<"\t"<<TadpoleFactorU0Error<<"\t"
	      <<BareFreeEnergy<<"\t"<<BareFreeEnergyError<<endl;
	  fout0.close();
	  
	    sprintf(FileNameOutput,"%s_Nt%d_Ns%d.txt",OutputFile1.c_str(), NtTemp, Ns );  //LO  T+V
	     if(k==0 ){ sprintf(FileNameOutputLO,    "%s",FileNameOutput);}
	     if(k==1 ){ sprintf(FileNameOutputVacLO,  "%s",FileNameOutput);}

	     fout.open(FileNameOutput, ios::app); double Small=pow(10,-22);
	     if(  k==0 )
	       {
		 cout<<"\n Square (Nt,Ns)="<<NtTemp<<","<<Ns<<": Beta="<<Beta[i]<<endl;
		 mean=0; err=0;
		 GM0 = new TGraph(FileName[0],"%lf %lf");
		 JackknifeMean(GM0,  &mean, &err, FirstRun, SkipN);
		 Ymean[0]=mean; Yerr[0]=err;
		 
		 cout<<" Rectangle :"<<endl;
		 mean=0; err=0;
		 GM0 = new TGraph(FileName[0],"%lf %*lf %*lf %lf");
		 JackknifeMean(GM0, &mean, &err, FirstRun, SkipN);
		 Ymean[1]=mean; Yerr[1]=err;		

		 fout<<Beta[i]<<"\t"<<Ymean[0]<<"\t"<<Yerr[0]<<"\t"<<Ymean[1]<<"\t"<<Yerr[1]<<endl;
		 cout<<Beta[i]<<"\t"<<Ymean[0]<<"\t"<<Yerr[0]<<"\t"<<Ymean[1]<<"\t"<<Yerr[1]<<endl;	
	       }
	     else {
	       fout<<Beta[i]<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
	       cout<<Beta[i]<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
	     }
	     fout.close();
	}   //end of LO or  NNLO or NNNNLO Fmunu loop 0-1
    }	//end of different beta files loop 0-9
           //end of thermal or vacuum loop 0-1
  
  
  //Plot for Average Plaquette, TadpoleFactor, Bare free energy
  TCanvas *c0 = new TCanvas("c0","Average Plaquette and Bare Free energy",800,400);
  c0->Divide(2,1);c0->cd(1);
  TGraphErrors *gAP  = new TGraphErrors(OutputFile0.c_str(), "%lf  %lf  %lf");
  TGraphErrors *gAP2 = new TGraphErrors(OutputFile00.c_str(),"%lf  %lf  %lf"); //Vacuum
  TGraphErrors *gU0  = new TGraphErrors(OutputFile0.c_str(), "%lf  %*lf %*lf %lf %lf");
  TGraphErrors *gU02 = new TGraphErrors(OutputFile00.c_str(),"%lf  %*lf %*lf %lf %lf"); //Vacuum
  Custom(gAP,2); Custom(gAP2,1); Custom(gU0,8); Custom(gU02,4);
  gAP->SetMarkerStyle(22);  gAP2->SetMarkerStyle(25);  gU0->SetMarkerStyle(33);  gU02->SetMarkerStyle(24);
  gAP->Draw("AL p");gAP2->Draw("p same ");gU0->Draw("Lp same "); gU02->Draw("p same ");
  gAP->GetYaxis()->SetRangeUser(0.0,2.0);
  gAP->GetYaxis()->SetTitle("Re<#it{O}> "); gAP->GetXaxis()->SetTitle("#beta_{0}");
  gPad->SetMargin(0.2,0.02,0.2,0.1); gAP->GetYaxis()->SetTitleOffset(1.6);
  TLEGEND4(0.2,0.54,0.88,0.92, gAP,"Average Plaquette, #it{n_{#tau}.n_{s}^{3}}","elp", gAP2,"Average Plaquette, #it{n_{s}.n_{s}^{3}}","ep", gU0,"#it{u}_{0}, #it{n_{#tau}.n_{s}^{3}}","elp", gU02,"#it{u}_{0}, #it{n_{s}.n_{s}^{3}}","ep", 0.05);
  c0->cd(2);
  TGraphErrors *gBE  = new TGraphErrors(OutputFile0.c_str(), "%lf  %*lf %*lf %*lf %*lf  %lf %lf");
  TGraphErrors *gBE2 = new TGraphErrors(OutputFile00.c_str(),"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
  Custom(gBE,2); Custom(gBE2,1);  gBE->Draw("AL p");gBE2->Draw("p same ");
  gBE->SetMarkerStyle(22);  gBE2->SetMarkerStyle(23);
  double *MinMax=FindMinMax4(gBE,gBE2,gBE,gBE); double Diff = (MinMax[1]-MinMax[0])/6.0;
  gBE->GetYaxis()->SetRangeUser(MinMax[0]-2*Diff,MinMax[1]+4*Diff);
  gBE->GetYaxis()->SetTitle("Bare free energy (-log(|#it{P}|)) "); gBE->GetXaxis()->SetTitle("#beta_{0}");
  gPad->SetMargin(0.2,0.02,0.2,0.1); gBE->GetYaxis()->SetTitleOffset(1.6);
  TLEGEND2(0.5,0.7,0.88,0.88, gBE,"#it{n_{#tau}.n_{s}^{3}}","elp", gBE2,"#it{n_{s}.n_{s}^{3}}","ep", 0.05);
  sprintf(FileNameOutput,"Test/AP_u0_FreeEnergy%s_Nt%d.pdf", QCDTypeS1.c_str(),Nt);
//c0->SaveAs(FileNameOutput);
  
 //(Nt,Ns) //Re-part   LO1, LO2, NNLO1, NNLO2, NNNNLO1, NNNNLO2
 TGraphErrors *g  = new TGraphErrors(FileNameOutputLO,    "%lf  %lf  %lf");  
 TGraphErrors *g2 = new TGraphErrors(FileNameOutputLO,    "%lf  %*lf %*lf  %lf  %lf");

 //(Ns,Ns) //Re-part   LO1, LO2, NNLO1, NNLO2, NNNNLO1, NNNNLO2 
 TGraphErrors *gv  = new TGraphErrors(FileNameOutputVacLO,    "%lf  %lf  %lf");
 TGraphErrors *gv2 = new TGraphErrors(FileNameOutputVacLO,    "%lf  %*lf %*lf   %lf  %lf");

 TGraphErrors *stad, *stad2, *stad3, *stad4, *stad5, *stad6, *stad7, *stad8; //Real with tad-pole factor corrected

 stad=SubtractTGraphErrors(g,gv); 
 stad2=SubtractTGraphErrors(g2,gv2); 

 WriteToFile(stad, stad2,  OutputFile11);//Vacuum subtracted Real-part LO1 LO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t beta

 
 //convert to temperature dependent form and write to a file
 TGraphErrors *Gct, *Gct2, *Gct3, *Gct4, *Gct5, *Gct6; //Real, 1x1 + 1x2
 TGraphErrors *Gct1x1; //Real 1x1
 Gct    = Compute_AlphasFFOverT4_Vs_T(stad, stad2, Nt, QCDType,gU0, 0 );  //Real-part LO1, without 1/u0^4
 Gct2   = Compute_AlphasFFOverT4_Vs_T(stad, stad2, Nt, QCDType, gU0, 1);  // with 1/u0^4
 Gct1x1 = Compute_1x1_AlphasFFOverT4_Vs_T(stad, stad2, Nt, QCDType);
 WriteToFile(Gct, Gct2, Gct1x1, OutputFile20);  //Vacuum subtracted Real-part 1x1+1x2, 1x1+u0*1x2,  1x1 w.r.t T

 NtNs="Nt"+std::to_string(Nt)+"Ns"+std::to_string(Ns);
 TCanvas *c1 = new TCanvas(NtNs.c_str(),NtNs.c_str(),900,800);
 double SmallPad=pow(10,-54);
 c1->Divide(2,2,SmallPad,SmallPad);
 c1->SetFillStyle(4000);
 c1->cd(1);
 Custom(g,2); Custom(g2,4);
 g->SetMarkerStyle(22); g2->SetMarkerStyle(23); 
 g->Draw("AL p");g2->Draw("p same ");
 MinMax=FindMinMax2(g,g2); Diff = (MinMax[1]-MinMax[0])/6.0;
 g->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 g->GetYaxis()->SetTitle("T+V: Real<#it{O}1x1 and 1x2>"); g->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.2,0.02,0.2,0.1); g->GetYaxis()->SetTitleOffset(1.6);

 c1->cd(2);
 Custom(gv,2); Custom(gv2,4);
 gv->SetMarkerStyle(22);  gv2->SetMarkerStyle(23);
 gv->Draw("AL p");gv2->Draw("Lp same ");
 MinMax=FindMinMax2(gv,gv2); Diff = (MinMax[1]-MinMax[0])/6.0;
 gv->GetYaxis()->SetRangeUser(MinMax[0] -Diff, MinMax[1]+Diff);
 gv->GetYaxis()->SetTitle("Vacuum: Real<#it{O} 1x1 and 1x2> "); gv->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.2,0.02,0.12,0.1); gv->GetYaxis()->SetTitleOffset(1.6);
//sprintf(FileNameOutput,"Test/BareFmunuFmunuVsBeta%s_Nt%d.pdf", QCDTypeS1.c_str(),Nt); 
//c1->SaveAs(FileNameOutput); 
 
 c1->cd(3);
 Custom(stad,2); Custom(stad2,4);
 stad->SetMarkerStyle(22);  stad2->SetMarkerStyle(23);
 stad->Draw("AL p");stad2->Draw("Lp same ");
 MinMax=FindMinMax2(stad,stad2); Diff = (MinMax[1]-MinMax[0])/6.0;
 stad->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 stad->GetYaxis()->SetTitle("(Vacuum subtracted): Real<#it{O} 1x1 1x2>"); stad->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.2,0.02,0.2,0.1); stad->GetYaxis()->SetTitleOffset(1.6);
  
 c1->cd(4);
 Custom(Gct,2); Custom(Gct2,4);
 Gct->SetMarkerStyle(22);  Gct2->SetMarkerStyle(23);
 Gct->Draw("AL p");Gct2->Draw("Lp same ");
 MinMax=FindMinMax2(Gct,Gct2); Diff = (MinMax[1]-MinMax[0])/6.0;
 Gct->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 Gct->GetYaxis()->SetTitle("(Vacuum subtracted): Real<#it{O} >"); Gct->GetXaxis()->SetTitle("#it{T}(MeV)");
 gPad->SetMargin(0.2,0.02,0.2,0.1); Gct->GetYaxis()->SetTitleOffset(1.6);
 
}

TGraphErrors * Compute_AlphasFFOverT4_Vs_T(TGraphErrors *g,TGraphErrors *g2, int Nt, int QCDType, TGraphErrors *gU0, int FlagU0 )
{
  cout<<"1x1+1x2:Started Convert to temperature  Calculation:: for Nt="<<Nt<<endl;
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  const int NU0=gU0->GetN(); double *XU0=gU0->GetX();double *YU0=gU0->GetY();double *EYU0=gU0->GetEY(); double *EXU0=gU0->GetEX();
  double x[N], y[N], ex[N], ey[N];
  double Beta, T, a0, a0T, g0Square,  PreFactor=0, FactorU0=0;
  for(int i=0; i<N; i++)
    {
      Beta = X[i];
      T  = LatticeSpacing(Nt, QCDType, Beta);
      a0 = 1.0/(Nt*T); g0Square = BareCouplingSquare(QCDType, Beta);
      a0T = 1.0/Nt;
      x[i] = T;
      ex[i]=0;
      if(FlagU0==0)  { FactorU0=1.0;}
      if(FlagU0==1) { FactorU0=1.0/pow(YU0[i],2);}

      PreFactor=10.0/pow(a0T,4);
      y[i]  = PreFactor*(  (Y[i]/3.0)  + (FactorU0*Y2[i]/60.0) );
      ey[i] = PreFactor*sqrt( pow( EY[i]/3.0,2)  +  pow(FactorU0*EY2[i]/60.0,2)   );
      cout<<"1x1+1x2:"<<FlagU0<<"\t,Nt="<<Nt<<",Three Beta="<<X[i]<<","<<X2[i]<<","<<XU0[i]<<",T="<<T<<endl;
      cout<<"ConvertToT: PreFactor="<<PreFactor<<", U0="<<YU0[i]<<", Y="<<y[i]<<", ey"<<ey[i]<<endl;
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

TGraphErrors * Compute_1x1_AlphasFFOverT4_Vs_T(TGraphErrors *g,TGraphErrors *g2, int Nt, int QCDType)
{					   
  cout<<"1x1: Started Convert to temperature  Calculation:: for Nt="<<Nt<<endl;
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  const int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
						 
  double x[N], y[N], ex[N], ey[N];
  double Beta, T, a0, a0T, g0Square,  PreFactor=0, FactorU0=0;
  for(int i=0; i<N; i++)
    {
      Beta = X[i];
      T  = LatticeSpacing(Nt, QCDType, Beta);
      a0 = 1.0/(Nt*T); g0Square = BareCouplingSquare(QCDType, Beta);
      a0T = 1.0/Nt;
      x[i] = T;
      ex[i]=0;      
      PreFactor=6.0/pow(a0T,4);
      y[i]  = PreFactor*(Y[i]/3.0);
      ey[i] = PreFactor*( EY[i]/3.0);
    cout<<"1x1: \t,Nt="<<Nt<<", Beta="<<X[i]<<", T="<<T<<endl;
    cout<<"ConvertToT: PreFactor="<<PreFactor<<", Y="<<y[i]<<", ey"<<ey[i]<<endl;
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}


//...............................

double BareCouplingSquare(int QCDType, double Beta)
{
  double g0Square=0.0;
  if(QCDType>0) {g0Square=10.0/Beta;}
  if(QCDType<=0){g0Square=6.0/Beta;}
  return g0Square;
}

double StrongCouplingOneLoop(int QCDType, int Nt, double T)
{
  int Nf=0; 
  double QSquare = pow(3.1415926*Nt*T, 2.0);// MeV^2
  double LambdaSquare = pow(340, 2.0); // in MS-bar scheme, MeV^2
  if(QCDType>0 ){Nf=3;}
  else {Nf=3;}
  double beta0= 11 - (2.0*Nf/3.0); 
  double alphas = 4*3.1415926/(beta0*log(QSquare/LambdaSquare));
  return alphas;
}

double RbetaFmunuSquare(int QCDType, double beta)
{ double RBETA=1.0;
  TGraph *g;
  if(QCDType>0){ g= new TGraph("HISQ_RBeta/R_beta_HISQActionHOTQCD.txt","%lf %lf");}
  else { g= new TGraph("HISQ_RBeta/R_beta_PureSU3GaugeAction.txt","%lf %lf"); }
  
  RBETA = g->Eval(beta);
  return RBETA;
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
      //Old: double T[10]   = {194.765,  257.954,  272.89,  305.434,  341.898,   404.999,  507.794,     601.833,  713.438,  799.206};
      //Karsch
      //double T[10]   = {135.10,   221.80,   242.87,   284.31,    333.55,    412.02,    539.30,   657.06,   801.33,   915.22};
      FileCount=10; Beta = new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] = BETAS[i]; Temperature[i] = T[i];
        }
    }  

  if(QCDType>0 && Nt==4)
    {
      double BETA[8] ={5.9,      6.0,       6.285,     6.515,    6.664,    6.95,      7.15,     7.373};
      string BETAS[8]={"5.9000", "6.0000",  "6.2850",  "6.5150", "6.6640", "6.9500",  "7.1500", "7.3730"};
      double MS[8]  ={0.132,     0.1138,    0.079,     0.0603,   0.0514,   0.0386,    0.032,    0.025};
      double T[8]   ={201,    221,     291,     364,     421,     554,     669,     819};
      FileCount=8; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }

  if(QCDType>0 && Nt==6)
    {
      double BETA[10] ={6,        6.215,    6.285,    6.423,    6.664,    6.95,     7.15,      7.373,    7.596,      7.825};
      string BETAS[10]={"6.0000", "6.2150", "6.2850", "6.4230", "6.6640", "6.9500", "7.1500",  "7.3730", "7.5960",   "7.8250"};
      double MS[10]  ={0.1138,  0.0862,   0.079,  0.067,   0.0514,  0.0386,  0.032,   0.025,   0.0202,  0.0164};
      double T[10]   ={147,      181,      194,    222,     281,     370,     446,     547,     667,     815};
      FileCount=10; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }
  
  if(QCDType>0 && Nt==8)
    {
      //6.2850,  7.8250                                                                       
      double BETA[9] ={  6.515,   6.575,   6.664,   6.95,    7.28,    7.5,     7.596, 7.8250,    8.2};
      string BETAS[9]={ "6.5150", "6.5750",  "6.6640",  "6.9500", "7.2800",  "7.5000", "7.5960", "7.8250",    "8.2000"};
      //0.0790,  0.0164
      double MS[9]  ={  0.0604,  0.0564,  0.0514,  0.0386,  0.0284,  0.0222,  0.0202, 0.0164,   0.01167};
      //146, 611
      double T[9]   ={    182,     193,     211,     277,     377,     459,     500, 611,       843};
      FileCount=9; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }

  TGraph *g = new TGraph(FileCount,Beta,Temperature);
  return g->Eval(beta);
}

double TadpoleFactorTUMQCD(double beta)
{
  double BETA[28]={6.488, 6.515, 6.608, 6.664,  6.7, 6.74, 6.77, 6.8, 6.84, 6.88, 6.91, 6.95, 7.03, 7.1, 7.15, 7.2, 7.28,                     7.373, 7.596, 7.65,  7.825,  8.0, 8.2,  8.4,  8.57, 8.71, 8.85, 9.06};
  double u0[28]  ={0.86000, 0.86453, 0.86823, 0.87026, 0.87152, 0.87288, 0.87388, 0.87485, 0.87613, 0.87736, 0.87827,                         0.87945, 0.88173, 0.88363, 0.88493, 0.88621, 0.88817, 0.89035, 0.89517, 0.89627, 0.89962, 0.90274,                         0.90604, 0.90909, 0.91152, 0.91341, 0.91522, 0.91778};
  TGraph *g = new TGraph(28,BETA,u0);
  return g->Eval(beta);  
}

void VacuumResultNt32(double beta, int j, int WithOrWithoutTadpoleCorrection, double *MeanAndErrorFit)
{
  double BETA[9]={6.515,  6.575, 6.664, 6.95, 7.28, 7.5, 7.596, 7.825, 8.2};
  double Y[9];
  double EY[9];
  
  if(j==0)
    {
      double y[9]={-8.52176e-06, -8.52176e-06, -8.52176e-06, -8.52176e-06, -8.52176e-06, -8.52176e-06, -8.52176e-06, -8.52176e-06, -8.52176e-06};
      double ey[9]={9.98077e-05, 9.98077e-05, 9.98077e-05, 9.98077e-05, 9.98077e-05, 9.98077e-05, 9.98077e-05, 9.98077e-05, 9.98077e-05};
      for(int i=0;i<9;i++)
	{ Y[i]=y[i]; EY[i]=ey[i]*0.01;}
    }

  if(WithOrWithoutTadpoleCorrection==0)
    {
      if(j==1)
	{
	  double y[9]={0.598544, 0.587548, 0.571237, 0.526876, 0.484372, 0.459823, 0.449666, 0.427772, 0.391919};
	  double ey[9]={0.000480284, 0.000334801, 0.000119, 4.35e-05, 8.15314e-05, 0.000760121, 0.001405, 0.00411071, 0.00854147};
	  for(int i=0;i<9;i++)
	    { Y[i]=y[i]; EY[i]=ey[i];}
	}
      if(j==2)
	{ double y[9] ={-2.44367, -2.40279, -2.34214, -2.17658, -2.01713, -1.92438, -1.88585,  -1.80264,  -1.66638};
	  double ey[9]={0.0018001, 0.00122624, 0.000375, 6e-05, 0.000181256, 0.00287471, 0.00518, 0.0150281, 0.0311549}; 
	  for(int i=0;i<9;i++)
	    { Y[i]=y[i]; EY[i]=ey[i];}
	}            
    }
  

  if(WithOrWithoutTadpoleCorrection==1)
    {
      if(j==1)
        {
          double y[9] ={0.801303, 0.782359, 0.754259, 0.681215, 0.614286, 0.576623, 0.561164, 0.528605, 0.475288};
          double ey[9]={0.000905716, 0.000604824, 0.0001585, 2.7e-05, 7.65022e-05, 0.00145683, 0.002614, 0.0075324, 0.0155866};
	  for(int i=0;i<9;i++)
	    { Y[i]=y[i]; EY[i]=ey[i];}
        }
      if(j==2)
	{ double y[9] ={-3.29361, -3.22238, -3.11673, -2.84034, -2.58505, -2.44012, -2.38041, -2.2541, -2.04726};
          double ey[9]={0.00324861, 0.00212367, 0.000455, 1.5e-05, 0.0001363, 0.0055154, 0.009745, 0.0279004, 0.0576309};
	  for(int i=0;i<9;i++)
	    { Y[i]=y[i]; EY[i]=ey[i];}
        }
    }
  TGraph *g = new TGraph(9,BETA, Y);
  TGraph *G = new TGraph(9,BETA, EY);
  MeanAndErrorFit[0] = g->Eval(beta); MeanAndErrorFit[1] = G->Eval(beta);

}
