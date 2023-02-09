#include <stdlib.h>
#include "StatTool.h"
#include "Custom.h"

void Mean(TGraph *g, double *Ymean, double *Yerr, int FirstPoint, int SkipN);
void AnalyzeNtNs(int Nt,int Ns,int QCDType,int FmunuType, int NNLOFlag, double Qminus,  string File0, string File00, string File1,  string File11, string File12, string File13, string File14, string File21, string File22); 
TGraphErrors * ConvertToT(TGraphErrors *g , int Nt, int QCDType, double Qminus, string LOorNLO);
void QhatRe(int Nt, int QCDType, double *Qminus, string File1Output, string File3Output);
void QhatIm(int Nt, int QCDType, double *Qminus, string File1Output, string File3Output);
double BareCouplingSquare(int QCDType, double Beta);
double RbetaFmunuSquare(int QCDType, double beta);
double LatticeSpacing(int Nt, int QCDType, double beta);

void StatisticalAverageAll()
{
  gSystem->Exec("./Clean.sh");
  int Nt=4;
  int Ns=Nt*4;
  int QCDType=1, FmunuType=4, NNLOFlag=1;
  double Qminus[5]={10000, 20000, 50000, 80000, 100000}; //in MeV
  //Plaquette, tadpolefactor, and bare free energy
  char Output0[10000], Output00[10000];
  //   T+V LO,         T+V NLO,        V LO            V NLO (w.r.t beta)
  char Output1[10000], Output2[10000], Output3[10000], Output4[10000];
  //   Vacuum Sub_LO    Vacuum Sub_NLO (w.r.t beta)
  char Output11[10000], Output12[10000], Output13[10000], Output14[10000];
  //   Vacuum Sub_LO    Vacuum Sub_NLO (Operators dimensionless w.r.t T) 
  char Output21[10000], Output22[10000];
  //   qhatRe           qhatIm         (w.r.t T)
  char Output31[10000], Output32[10000];
  string QCDTypeS, FmunuTypeS;
  if(QCDType<=0){QCDTypeS="PureGauge";}
  if(QCDType>0){QCDTypeS="FullQCD";}
  if(FmunuType==2){FmunuTypeS="SinglePlaquetteTraceless";}
  if(FmunuType==4){FmunuTypeS="CloverTraceless";}

  sprintf(Output0, "Test/%sDataPloopNt%d_Ns%d.txt",QCDTypeS.c_str(),Nt,Ns);
  sprintf(Output00,"Test/%sDataPloopNt%d_Ns%d.txt",QCDTypeS.c_str(),Ns,Ns);  
  sprintf(Output1, "Test/%s%sOperatorVsBeta",QCDTypeS.c_str(), FmunuTypeS.c_str());
  //sprintf(Output2, "Test/%s%sOperatorVsBeta_NLO_Nt%d_Ns%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt,Ns);
  //sprintf(Output3, "Test/%s%sOperatorVsBeta_LO_Nt%d_Ns%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Ns,Ns);
  //sprintf(Output4, "Test/%s%sOperatorVsBeta_NLO_Nt%d_Ns%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Ns,Ns);

  sprintf(Output11,"Test/%s%sVacuumSubtractedOperatorVsBeta_RealPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt); //Vacuum subtracted Real-part LO1 LO2 NLO1 NLO2  w.r.t beta
  sprintf(Output12,"Test/%s%sVacuumSubtractedOperatorVsBeta_ImagPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt);//Vacuum subtracted Im-part LO1 LO2 NLO1 NLO2  w.r.t beta
  sprintf(Output13,"Test/%s%sVacuumSubtractedOperatorTadpoleCorrectedVsBeta_RealPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt); //Vacuum subtracted Real-part LO1 LO2 NLO1 NLO2  w.r.t beta 
  sprintf(Output14,"Test/%s%sVacuumSubtractedOperatorTadpoleCorrectedVsBeta_ImagPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt);//Vacuum subtracted Im-part LO1 LO2 NLO1 NLO2  w.r.t beta
  sprintf(Output21,"Test/%s%sVacuumSubtractedOperatorReNormalizedVsT_RealPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt);
  sprintf(Output22,"Test/%s%sVacuumSubtractedOperatorReNormalizedVsT_ImagPartNt%d.txt", QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt);
  AnalyzeNtNs(Nt, Ns, QCDType,FmunuType, NNLOFlag, Qminus[0], Output0, Output00, Output1,  Output11, Output12, Output13, Output14,  Output21, Output22);
  
  //DimensionlessOperatorVsT(Nt, Qminus, Output11, Output12,   Output21, Output22);
  sprintf(Output31,"Test/QhatRealVsT_Nt%d_Type%s%s",Nt,QCDTypeS.c_str(),FmunuTypeS.c_str());
  sprintf(Output32,"Test/QhatImVsT_Nt%d_Type%s%s",Nt, QCDTypeS.c_str(),FmunuTypeS.c_str());
  QhatRe(Nt, QCDType, &Qminus[0], Output13,  Output31);
  QhatIm(Nt, QCDType, &Qminus[0], Output13,  Output32);  

}

void AnalyzeNtNs(int Nt, int Ns, int QCDType, int FmunuType, int NNLOFlag, double Qminus, string OutputFile0, string OutputFile00, string OutputFile1, string OutputFile11, string OutputFile12, string OutputFile13, string OutputFile14,string OutputFile21, string OutputFile22)
{
  CustomGlobal();
  int FirstRun=500; int SkipN=1;if(QCDType<=0){SkipN=0;}
  string File[1000],File2[1000];int Nf;
  string LOorNLO;
  string  NtNs, QCDTypeS, QCDTypeS1, QCDTypeS2; 

  if(QCDType<=0){QCDTypeS="PureGauge";Nf=1.0;}
  if(QCDType>0){QCDTypeS="FullQCD";Nf=1.0;}
  if(FmunuType==1) { QCDTypeS1=QCDTypeS+"SinglePlaquette";         QCDTypeS2="SinglePlaquette_";}
  if(FmunuType==2) { QCDTypeS1=QCDTypeS+"SinglePlaquetteTraceless"; QCDTypeS2="SinglePlaquette_Traceless_";}
  if(FmunuType==3) { QCDTypeS1=QCDTypeS+"Clover";                  QCDTypeS2="Clover_";}
  if(FmunuType==4) { QCDTypeS1=QCDTypeS+"CloverTraceless";         QCDTypeS2="Clover_Traceless_";}
  int FileCount,NtTemp; double *Beta; string *BetaS; double *Temperature; char NtTempS[100], NsTempS[100];  

  if(QCDType<=0 && Nt==4)
    {
      double BETA[10] = {5.35,     5.6,      5.7,      5.8,      5.9,      6,        6.2,      6.35,     6.5,      6.6};
      string BETAS[10]= {"5.3500", "5.6000", "5.7000", "5.8000", "5.9000", "6.0000", "6.2000", "6.3500", "6.5000", "6.6000"};
      double T[10]   =  {191.423,  253.235,  283.295,  316.964,  354.68,   396.931,  497.306,  589.089,  697.975,  781.629};
      //Karsch
      //double T[10]   ={137.753,  220.99,   270.21,   333.316,  405.588,  485.75,   667.11,   824.04,   1010.12,  1151.89};
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
    {
      double BETA[10] ={6.285,  6.515,   6.575,   6.664,   6.95,    7.28,    7.5,     7.596,   7.825,   8.2};
      string BETAS[10]={"6.2850", "6.5150", "6.5750",  "6.6640",  "6.9500", "7.2800",  "7.5000", "7.5960", "7.8250",   "8.2000"};
      double MS[10]  ={0.079,  0.0604,  0.0564,  0.0514,  0.0386,  0.0284,  0.0222,  0.0202,  0.0164,  0.01167};
      double T[10]   ={146,    182,     193,     211,     277,     377,     459,     500,     611,     843};
      FileCount=10; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }

  ofstream fout; char FileNameOutput[10000];char FileName[10000], FileName2[10000];
  char FileNameOutputLO[10000], FileNameOutputNLO[10000], FileNameOutputNNLO[10000], FileNameOutputNNNNLO[10000];
  char FileNameOutputVacLO[10000],FileNameOutputVacNLO[10000], FileNameOutputVacNNLO[10000], FileNameOutputVacNNNNLO[10000];
  TGraph *GM;
  vector<double> Ymean , Yerr;
  
  for(int k=0; k<2;k++)                   //Thermal or Vacuum
    {
      for(int i=0; i<FileCount; i++)      // different beta
	{ 
	  for(int j=0; j<4;j++)           // LO or NLO or NNLO or NNNNLO Fmunu
	    { if(j==0) {LOorNLO="LO_";}	    
	      if(j==1) {LOorNLO="NLO_";}
	      if(j==2) {LOorNLO="NNLO_";}
	      if(j==3) {LOorNLO="NNNNLO_";}
	      if(QCDType<1 && j>1) {LOorNLO="NLO_";}
	      NtTemp=Nt;
	      if(k==1){NtTemp=Ns;}//NtNs= "Nt16_Ns16_"
	      sprintf(NtTempS,"%d",NtTemp);
	      sprintf(NsTempS,"%d",Ns);
	      File[i]=std::string("DataTraceFmunu") + LOorNLO + QCDTypeS2 + "Nt" + NtTempS + "_Ns"+ NsTempS  +"_" + "Beta"+ BetaS[i] + ".txt";
	      sprintf(FileName,"RawData%s/Nt%d_Ns%d/%s",QCDTypeS1.c_str(), NtTemp, Ns, File[i].c_str());
	      //RawDataFullQCDCloverTraceless/Nt4_Ns16/DataTraceFmunuLO_Clover_Traceless_Nt4_Ns16_Beta5.9000.txt 
              File2[i]=std::string("DataPloopNt") + NtTempS + "_Ns"+ NsTempS  + "_" + "Beta"+ BetaS[i] + ".txt";
              sprintf(FileName2,"RawData%s_Plaquette_PLoop/Nt%d_Ns%d/%s",QCDTypeS.c_str(), NtTemp, Ns, File2[i].c_str());
              //RawDataFullQCD_Plaquette_PLoop/Nt4_Ns16/DataPloopNt4_Ns16_Beta5.9000.txt
	      Ymean.resize(0); Yerr.resize(0);
	      double mean=0.0, err=0.0;

	 cout<<"\n F3i"<<LOorNLO<<"F3i.Re (Nt,Ns)="<<NtTemp<<","<<Ns<<": Beta="<<Beta[i]<<endl;
	 GM = new TGraph(FileName,"%lf %lf");
	 Mean(GM, &mean, &err, FirstRun, SkipN);
	 Ymean.push_back(mean); Yerr.push_back(err);
	 
	 cout<<" F3i"<<LOorNLO<<"F3i.Im"<<endl;
	 GM = new TGraph(FileName,"%lf %*lf %lf");
	 Mean(GM, &mean, &err, FirstRun, SkipN);
	 Ymean.push_back(mean); Yerr.push_back(err);
	 
	 cout<<" F3i"<<LOorNLO<<"F4i.Re "<<endl;
	 GM = new TGraph(FileName,"%lf %*lf %*lf     %*lf %*lf     %lf");
	 Mean(GM, &mean, &err, FirstRun, SkipN);
	 Ymean.push_back(mean); Yerr.push_back(err);
	 
	 cout<<" F3i"<<LOorNLO<<"F4i.Im"<<endl;
	 GM = new TGraph(FileName,"%lf %*lf %*lf     %*lf %*lf     %*lf %lf ");
	 Mean(GM, &mean, &err, FirstRun, SkipN);
	 Ymean.push_back(mean); Yerr.push_back(err);
	 
	 if(j==0) //LO
	   {
	     cout<<"\n Computing Avg plaq Tr(U1234)/3.0, tadpole factor u0, and bare free energy=-log(|P|) ="<<NtTemp<<","<<Ns<<endl;
	     GM = new TGraph(FileName2,"%lf %lf");
	     Mean(GM, &mean, &err, FirstRun, SkipN);
	     double AveragePlaquette = mean;  double AveragePlaquetteError = err;
	     double TadpoleFactorU0= pow(AveragePlaquette, 1.0/4.0); double TadpoleFactorU0Error=TadpoleFactorU0*err/(4.0*mean);
	     GM = new TGraph(FileName2,"%lf %*lf %*lf %*lf %*lf   %*lf %*lf %*lf %*lf   %*lf %lf %*lf %*lf %*lf");
	     Mean(GM, &mean, &err, FirstRun, SkipN);
	     double BareFreeEnergy = mean; double BareFreeEnergyError= err;
	     if(k==0){sprintf(FileNameOutput,"%s",OutputFile0.c_str());}  //Thermal
	     if(k==1){sprintf(FileNameOutput,"%s",OutputFile00.c_str());} //Vacuum
	     ofstream fout0; fout0.open(FileNameOutput, ios::app);
	     fout0<<Beta[i]<<"\t"<<AveragePlaquette<<"\t"<<AveragePlaquetteError<<"\t"
		 <<TadpoleFactorU0<<"\t"<<TadpoleFactorU0Error<<"\t"
	         <<BareFreeEnergy<<"\t"<<BareFreeEnergyError<<endl;
	     cout<<Beta[i]<<"\t"<<AveragePlaquette<<"\t"<<AveragePlaquetteError<<"\t"
		 <<TadpoleFactorU0<<"\t"<<TadpoleFactorU0Error<<"\t"
	         <<BareFreeEnergy<<"\t"<<BareFreeEnergyError<<endl;
	     fout0.close();
	   }

	 sprintf(FileNameOutput,"%s_%sNt%d_Ns%d.txt",OutputFile1.c_str(), LOorNLO.c_str(),NtTemp, Ns );  //LO  T+V 
	 if(k==0 && j==0){ sprintf(FileNameOutputLO, "%s",FileNameOutput);}
         if(k==0 && j==1){ sprintf(FileNameOutputNLO, "%s",FileNameOutput);}
	 if(k==0 && j==2){ sprintf(FileNameOutputNNLO, "%s",FileNameOutput);}
	 if(k==0 && j==3){ sprintf(FileNameOutputNNNNLO, "%s",FileNameOutput);}
         if(k==1 && j==0){ sprintf(FileNameOutputVacLO, "%s",FileNameOutput);}
	 if(k==1 && j==1){ sprintf(FileNameOutputVacNLO, "%s",FileNameOutput);}
	 if(k==1 && j==2){ sprintf(FileNameOutputVacNNLO, "%s",FileNameOutput);}
	 if(k==1 && j==3){ sprintf(FileNameOutputVacNNNNLO, "%s",FileNameOutput);}

	 fout.open(FileNameOutput, ios::app); double Small=pow(10,-22);
	 if(QCDType<1 && j>1) //NLO, Vacuum
	 {fout<<Beta[i]<<"\t"<<0<<"\t"<<Small*10<<"\t"<<0<<"\t"<<Small*4<<"\t"<<0<<"\t"<<Small*0.2<<"\t"<<0<<"\t"<<Small<<endl;}
	 else {fout<<Beta[i]<<"\t"<<Ymean[0]*Nf<<"\t\t"<<Yerr[0]*Nf<<"\t"<<Ymean[1]*Nf<<"\t\t"<<Yerr[1]*Nf<<"\t"
		   <<Ymean[2]*Nf<<"\t\t"<<Yerr[2]*Nf<<"\t"<<Ymean[3]*Nf<<"\t\t"<<Yerr[3]*Nf<<endl;
	   cout<<Beta[i]<<"\t"<<Ymean[0]*Nf<<"\t\t"<<Yerr[0]*Nf<<"\t"<<Ymean[1]*Nf<<"\t\t"<<Yerr[1]*Nf<<"\t"
	       <<Ymean[2]*Nf<<"\t\t"<<Yerr[2]*Nf<<"\t"<<Ymean[3]*Nf<<"\t\t"<<Yerr[3]*Nf<<endl;
	 }
	 fout.close();
	    }   //end of LO or NLO or NNLO or NNNNLO Fmunu loop 0-1
	}	//end of different beta files loop 0-9
    }           //end of thermal or vacuum loop 0-1
  

  //Plot for Average Plaquette, TadpoleFactor, Bare free energy
  TCanvas *c0 = new TCanvas("c0","Average Plaquette and Bare Free energy",800,400);
  c0->Divide(2,1);c0->cd(1);
  TGraphErrors *gAP  = new TGraphErrors(OutputFile0.c_str(), "%lf  %lf  %lf");
  TGraphErrors *gAP2 = new TGraphErrors(OutputFile00.c_str(),"%lf  %lf  %lf"); //Vacuum
  TGraphErrors *gU0  = new TGraphErrors(OutputFile0.c_str(), "%lf  %*lf %*lf %lf %lf");
  TGraphErrors *gU02 = new TGraphErrors(OutputFile00.c_str(),"%lf  %*lf %*lf %lf %lf"); //Vacuum
  Custom(gAP,2); Custom(gAP2,1); Custom(gU0,4); Custom(gU02,8);
  gAP->Draw("AL p");gAP2->Draw("p same ");gU0->Draw("Lp same "); gU02->Draw("p same ");
  gAP->GetYaxis()->SetRangeUser(0.0,2.0);
  gAP->GetYaxis()->SetTitle("Re<#it{O}> "); gAP->GetXaxis()->SetTitle("#beta");
  gPad->SetMargin(0.2,0.02,0.2,0.1); gAP->GetYaxis()->SetTitleOffset(1.6);
  TLEGEND4(0.2,0.54,0.88,0.92, gAP,"Average Plaquette, #it{n_{#tau}.n_{s}^{3}}","elp", gAP2,"Average Plaquette, #it{n_{s}.n_{s}^{3}}","elp", gU0,"#it{u}_{0}, #it{n_{#tau}.n_{s}^{3}}","elp", gU02,"#it{u}_{0}, #it{n_{s}.n_{s}^{3}}","elp", 0.05);
  c0->cd(2);
  TGraphErrors *gBE  = new TGraphErrors(OutputFile0.c_str(), "%lf  %*lf %*lf %*lf %*lf  %lf %lf");
  TGraphErrors *gBE2 = new TGraphErrors(OutputFile00.c_str(),"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
  Custom(gBE,2); Custom(gBE2,1);  gBE->Draw("AL p");gBE2->Draw("p same ");
  double *MinMax=FindMinMax4(gBE,gBE2,gBE,gBE); double Diff = (MinMax[1]-MinMax[0])/6.0;
  gBE->GetYaxis()->SetRangeUser(MinMax[0]-2*Diff,MinMax[1]+4*Diff);
  gBE->GetYaxis()->SetTitle("Bare free energy (-log(|#it{P}|)) "); gBE->GetXaxis()->SetTitle("#beta");
  gPad->SetMargin(0.2,0.02,0.2,0.1); gBE->GetYaxis()->SetTitleOffset(1.6);
  TLEGEND2(0.5,0.7,0.88,0.88, gBE,"#it{n_{#tau}.n_{s}^{3}}","elp", gBE2,"#it{n_{s}.n_{s}^{3}}","elp", 0.05);
  sprintf(FileNameOutput,"Test/AP_u0_FreeEnergy%s_Nt%d.pdf", QCDTypeS1.c_str(),Nt);
  c0->SaveAs(FileNameOutput);
  
 //(Nt,Ns) //Re-part   LO1, LO2, NLO1, NLO2, NNLO1, NNLO2, NNNNLO1, NNNNLO2
 TGraphErrors *g  = new TGraphErrors(FileNameOutputLO,    "%lf  %lf  %lf");  
 TGraphErrors *g2 = new TGraphErrors(FileNameOutputLO,    "%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *g3 = new TGraphErrors(FileNameOutputNLO,   "%lf  %lf  %lf");
 TGraphErrors *g4 = new TGraphErrors(FileNameOutputNLO,   "%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *g5 = new TGraphErrors(FileNameOutputNNLO,  "%lf  %lf  %lf");
 TGraphErrors *g6 = new TGraphErrors(FileNameOutputNNLO,  "%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *g7 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %lf  %lf");
 TGraphErrors *g8 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
           //Im-part   LO1, LO2, NLO1, NLO2
 TGraphErrors *G  = new TGraphErrors(FileNameOutputLO,"%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G2 = new TGraphErrors(FileNameOutputLO,"%lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G3 = new TGraphErrors(FileNameOutputNLO,"%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G4 = new TGraphErrors(FileNameOutputNLO,"%lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G5 = new TGraphErrors(FileNameOutputNNLO,"%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G6 = new TGraphErrors(FileNameOutputNNLO,"%lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G7 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G8 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf  %lf");
 //(Ns,Ns) //Re-part   LO1, LO2, NLO1, NLO2, NNLO1, NNLO2
 TGraphErrors *gv  = new TGraphErrors(FileNameOutputVacLO,"%lf  %lf  %lf");
 TGraphErrors *gv2 = new TGraphErrors(FileNameOutputVacLO,"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *gv3 = new TGraphErrors(FileNameOutputVacNLO,"%lf  %lf  %lf");
 TGraphErrors *gv4 = new TGraphErrors(FileNameOutputVacNLO,"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *gv5 = new TGraphErrors(FileNameOutputVacNNLO,"%lf  %lf  %lf");
 TGraphErrors *gv6 = new TGraphErrors(FileNameOutputVacNNLO,"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *gv7 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %lf  %lf");
 TGraphErrors *gv8 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
           //Im-part   LO1, LO2, NLO1, NLO2, NNLO1, NNLO2
 TGraphErrors *Gv  = new TGraphErrors(FileNameOutputVacLO,"%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *Gv2 = new TGraphErrors(FileNameOutputVacLO,"%lf  %*lf %*lf %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *Gv3 = new TGraphErrors(FileNameOutputVacNLO,"%lf  %*lf  %*lf %lf %lf");
 TGraphErrors *Gv4 = new TGraphErrors(FileNameOutputVacNLO,"%lf  %*lf %*lf %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *Gv5 = new TGraphErrors(FileNameOutputVacNNLO,"%lf  %*lf  %*lf %lf %lf");
 TGraphErrors *Gv6 = new TGraphErrors(FileNameOutputVacNNLO,"%lf  %*lf %*lf %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *Gv7 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %*lf  %*lf %lf %lf");
 TGraphErrors *Gv8 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %*lf %*lf %*lf %*lf %*lf %*lf  %lf  %lf");

 TGraphErrors *s, *s2, *s3, *s4, *s5, *s6, *s7, *s8; //Real no tad-pole factor 
 TGraphErrors *S, *S2, *S3, *S4, *S5, *S6, *S7, *S8; //Imag no tad-pole factor
 TGraphErrors *stad, *stad2, *stad3, *stad4, *stad5, *stad6, *stad7, *stad8; //Real with tad-pole factor corrected
 TGraphErrors *Stad, *Stad2, *Stad3, *Stad4, *Stad5, *Stad6, *Stad7, *Stad8; //Imag with tad-pole factor corrected

 s=SubtractTGraphErrors(g,gv); s2=SubtractTGraphErrors(g2,gv2); s3=SubtractTGraphErrors(g3,gv3);s4=SubtractTGraphErrors(g4,gv4); s5=SubtractTGraphErrors(g5,gv5); s6=SubtractTGraphErrors(g6,gv6); s7=SubtractTGraphErrors(g7,gv7); s8=SubtractTGraphErrors(g8,gv8);//Real-part LO1 LO2 NLO1 NLO2 NNLO1 NNLO2  NNNNLO1 NNNNLO2
 S=SubtractTGraphErrors(G,Gv); S2=SubtractTGraphErrors(G2,Gv2); S3=SubtractTGraphErrors(G3,Gv3);S4=SubtractTGraphErrors(G4,Gv4); S5=SubtractTGraphErrors(G5,Gv5); S6=SubtractTGraphErrors(G6,Gv6); S7=SubtractTGraphErrors(G7,Gv7);S8=SubtractTGraphErrors(G8,Gv8); //Im-part   LO1 LO2 NLO1 NLO2 NNLO1 NNLO2  NNNNLO1 NNNNLO2
 WriteToFile(s, s2, s3, s4, s5, s6, s7, s8, OutputFile11);//Vacuum subtracted Real-part LO1 LO2 NLO1 NLO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t beta
 WriteToFile(S, S2, S3, S4, S5, S6, S7, S8, OutputFile12);//Vacuum subtracted Im-part LO1 LO2 NLO1 NLO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t beta
 double TadPolePower=8;
 stad=SubtractTGraphErrors(DivideTGraphErrors(g,1,gU0,TadPolePower),DivideTGraphErrors(gv,1,gU02,TadPolePower)); stad2=SubtractTGraphErrors(DivideTGraphErrors(g2,1,gU0,TadPolePower),DivideTGraphErrors(gv2,1,gU02,TadPolePower)); stad3=SubtractTGraphErrors(DivideTGraphErrors(g3,1,gU0,TadPolePower),DivideTGraphErrors(gv3,1,gU02,TadPolePower));stad4=SubtractTGraphErrors(DivideTGraphErrors(g4,1,gU0,TadPolePower),DivideTGraphErrors(gv4,1,gU02,TadPolePower)); stad5=SubtractTGraphErrors(DivideTGraphErrors(g5,1,gU0,TadPolePower),DivideTGraphErrors(gv5,1,gU02,TadPolePower)); stad6=SubtractTGraphErrors(DivideTGraphErrors(g6,1,gU0,TadPolePower),DivideTGraphErrors(gv6,1,gU02,TadPolePower)); stad7=SubtractTGraphErrors(DivideTGraphErrors(g7,1,gU0,TadPolePower),DivideTGraphErrors(gv7,1,gU02,TadPolePower)); stad8=SubtractTGraphErrors(DivideTGraphErrors(g8,1,gU0,TadPolePower),DivideTGraphErrors(gv8,1,gU02,TadPolePower));//Real-part LO1 LO2 NLO1 NLO2 NNLO1 NNLO2 NNNNLO1 NNNNLO2
 Stad=SubtractTGraphErrors(DivideTGraphErrors(G,1,gU0,TadPolePower),DivideTGraphErrors(Gv,1,gU02,TadPolePower)); Stad2=SubtractTGraphErrors(DivideTGraphErrors(G2,1,gU0,TadPolePower),DivideTGraphErrors(Gv2,1,gU02,TadPolePower)); Stad3=SubtractTGraphErrors(DivideTGraphErrors(G3,1,gU0,TadPolePower),DivideTGraphErrors(Gv3,1,gU02,TadPolePower));Stad4=SubtractTGraphErrors(DivideTGraphErrors(G4,1,gU0,TadPolePower),DivideTGraphErrors(Gv4,1,gU02,TadPolePower)); Stad5=SubtractTGraphErrors(DivideTGraphErrors(G5,1,gU0,TadPolePower),DivideTGraphErrors(Gv5,1,gU02,TadPolePower)); Stad6=SubtractTGraphErrors(DivideTGraphErrors(G6,1,gU0,TadPolePower),DivideTGraphErrors(Gv6,1,gU02,TadPolePower)); Stad7=SubtractTGraphErrors(DivideTGraphErrors(G7,1,gU0,TadPolePower),DivideTGraphErrors(Gv7,1,gU02,TadPolePower)); Stad8=SubtractTGraphErrors(DivideTGraphErrors(G8,1,gU0,TadPolePower),DivideTGraphErrors(Gv8,1,gU02,TadPolePower)); //Im-part   LO1 LO2 NLO1 NLO2 NNLO1 NNLO2 NNNNLO1 NNNNLO2  
 WriteToFile(stad, stad2, stad3, stad4, stad5, stad6, stad7, stad8, OutputFile13);//Vacuum subtracted Real-part LO1 LO2 NLO1 NLO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 tad-pole corrected w.r.t beta     
 WriteToFile(Stad, Stad2, Stad3, Stad4, Stad5, Stad6, Stad7, Stad8, OutputFile14);//Vacuum subtracted Im-part LO1 LO2 NLO1 NLO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 tad-pole corrected w.r.t beta
 //convert to temperature dependent form and write to a file
 TGraphErrors *ct, *ct2, *ct3, *ct4, *ct5, *ct6, *ct7, *ct8; //Real
 TGraphErrors *CT, *CT2, *CT3, *CT4, *CT5, *CT6, *CT7, *CT8; //Imag
 ct  = ConvertToT(stad  , Nt, QCDType, Qminus, "LO");  //Real-part LO1
 ct2 = ConvertToT(stad2 , Nt, QCDType, Qminus, "LO");  //Real-part LO2
 ct3 = ConvertToT(stad3 , Nt, QCDType, Qminus, "NLO"); //Real-part NLO1
 ct4 = ConvertToT(stad4 , Nt, QCDType, Qminus, "NLO"); //Real-part NLO2
 ct5 = ConvertToT(stad5 , Nt, QCDType, Qminus, "NNLO"); //Real-part NNLO1
 ct6 = ConvertToT(stad6 , Nt, QCDType, Qminus, "NNLO"); //Real-part NNLO2
 ct7 = ConvertToT(stad7 , Nt, QCDType, Qminus, "NNNNLO"); //Real-part NNNNLO1
 ct8 = ConvertToT(stad8 , Nt, QCDType, Qminus, "NNNNLO"); //Real-part NNNNLO2
 CT  = ConvertToT(Stad  , Nt, QCDType, Qminus, "LO");  //Im-part LO1
 CT2 = ConvertToT(Stad2 , Nt, QCDType, Qminus, "LO");  //Im-part LO2
 CT3 = ConvertToT(Stad3 , Nt, QCDType, Qminus, "NLO"); //Im-part NLO1
 CT4 = ConvertToT(Stad4 , Nt, QCDType, Qminus, "NLO"); //Im-part NLO2
 CT5 = ConvertToT(Stad5 , Nt, QCDType, Qminus, "NNLO"); //Im-part NNLO1
 CT6 = ConvertToT(Stad6 , Nt, QCDType, Qminus, "NNLO"); //Im-part NNLO2
 CT7 = ConvertToT(Stad7 , Nt, QCDType, Qminus, "NNNNLO"); //Im-part NNNNLO1
 CT8 = ConvertToT(Stad8 , Nt, QCDType, Qminus, "NNNNLO"); //Im-part NNNNLO2
 WriteToFile(ct, ct2, ct3, ct4, ct5, ct6, ct7, ct8, OutputFile21);  //Vacuum subtracted Real-part LO1 LO2 NLO1 NLO2 NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t T
 WriteToFile(CT, CT2, CT3, CT4, CT5, CT6, CT7, CT8, OutputFile22);  //Vacuum subtracted   Im-part LO1 LO2 NLO1 NLO2 NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t T 

 NtNs="Nt"+std::to_string(Nt)+"Ns"+std::to_string(Ns);
 TCanvas *c1 = new TCanvas(NtNs.c_str(),NtNs.c_str(),900,800);
 double SmallPad=pow(10,-54);
 c1->Divide(2,2,SmallPad,SmallPad);
 c1->SetFillStyle(4000);
 c1->cd(1);
 Custom(g,2); Custom(g2,41); Custom(g3,28); Custom(g4,8); Custom(g5,4); Custom(g6,6); Custom(g7,1); Custom(g8,46);
 g->SetMarkerStyle(22);  g2->SetMarkerStyle(23);  g3->SetMarkerStyle(25);  g4->SetMarkerStyle(3); g5->SetMarkerStyle(35);  g6->SetMarkerStyle(37); g7->SetMarkerStyle(24); g8->SetMarkerStyle(27);
 g->Draw("AL p");g2->Draw("p same ");g3->Draw("p same "); g4->Draw("p same "); g5->Draw("p same "); g6->Draw("p same"); g7->Draw("p same "); g8->Draw("p same ");
 MinMax=FindMinMax8(g,g2,g3,g4,g5,g6,g7,g8); Diff = (MinMax[1]-MinMax[0])/6.0;
 g->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 g->GetYaxis()->SetTitle("T+V: Real<#it{O}>"); g->GetXaxis()->SetTitle("#beta");
 gPad->SetMargin(0.2,0.02,0.2,0.1); g->GetYaxis()->SetTitleOffset(1.6);

 c1->cd(2);
 Custom(G,2); Custom(G2,41); Custom(G3,28); Custom(G4,8); Custom(G5,4); Custom(G6,6); Custom(G7,1); Custom(G8,46);
 G->SetMarkerStyle(22);  G2->SetMarkerStyle(23);  G3->SetMarkerStyle(25); G4->SetMarkerStyle(3); G5->SetMarkerStyle(35); G6->SetMarkerStyle(37); G7->SetMarkerStyle(24); G8->SetMarkerStyle(27);
 G->Draw("ALp "); G2->Draw(" Lp same"); G3->Draw("Lp same "); G4->Draw("Lp same "); G5->Draw("Lp same "); G6->Draw("Lp same "); G7->Draw(" Lp same"); G8->Draw(" Lp same");
 MinMax=FindMinMax8(G,G2,G3,G4,G5,G6,G7,G8);Diff = (MinMax[1]-MinMax[0])/6.0;
 G->GetYaxis()->SetRangeUser(MinMax[0]-Diff, MinMax[1]+6*Diff);
 G->GetYaxis()->SetTitle("T+V: Imag<#it{O}>");G->GetXaxis()->SetTitle("#beta"); 
 gPad->SetMargin(0.16,0.02,0.2,0.1);G->GetYaxis()->SetTitleOffset(1.5);
 //TLEGEND4(0.25,0.5,0.88,0.88, G,"#it{g^{2}a^{4}(F3iF3i - F4iF4i)}","elp", G2,"#it{g^{2}a^{4}(F3iF4i + F4iF3i)}","elp", G3,"#it{g^{2}a^{5}(F3iDzF3i - F4iDzF4i)}","elp", G4,"#it{g^{2}a^{5}(F3iDzF4i + F4iDzF3i)}","elp", 0.05);
 TLEGEND6(0.25,0.5,0.88,0.88, G,"#it{g^{2}a^{4}(F3iF3i - F4iF4i)}","elp", G2,"#it{g^{2}a^{4}(F3iF4i + F4iF3i)}","elp", G3,"#it{g^{2}a^{5}(F3iDz F3i - F4iDzF4i)}","elp", G4,"#it{g^{2}a^{5}(F3iDzF4i + F4iDzF3i)}","elp", G5,"#it{g^{2}a^{6}(F3iD2zF3i - F4iD2zF4i)}","elp", G6,"#it{g^{2}a^{6}(F3iD2zF3i - F4iD2zF4i)}","elp", 0.05);
 c1->cd(3);
 Custom(gv,2); Custom(gv2,41); Custom(gv3,28); Custom(gv4,8); Custom(gv5,4); Custom(gv6,6); Custom(gv7,1); Custom(gv8,46);
 gv->SetMarkerStyle(22);  gv2->SetMarkerStyle(23);  gv3->SetMarkerStyle(25);  gv4->SetMarkerStyle(3); gv5->SetMarkerStyle(35); gv6->SetMarkerStyle(37); gv7->SetMarkerStyle(24); gv8->SetMarkerStyle(27);
 gv->Draw("AL p");gv2->Draw("Lp same ");gv3->Draw("Lp same "); gv4->Draw("Lp same "); gv5->Draw("Lp same "); gv6->Draw("Lp same "); gv7->Draw("Lp same "); gv8->Draw("Lp same ");
 MinMax=FindMinMax8(gv,gv2,gv3,gv4,gv5,gv6,gv7,gv8); Diff = (MinMax[1]-MinMax[0])/6.0;
 gv->GetYaxis()->SetRangeUser(MinMax[0] -Diff, MinMax[1]+Diff);
 gv->GetYaxis()->SetTitle("Vacuum: Real<#it{O} >"); gv->GetXaxis()->SetTitle("#beta");
 gPad->SetMargin(0.2,0.02,0.12,0.1); gv->GetYaxis()->SetTitleOffset(1.6);

 c1->cd(4);
 Custom(Gv,2); Custom(Gv2,41); Custom(Gv3,28); Custom(Gv4,8); Custom(Gv5,4); Custom(Gv6,6); Custom(Gv7,1); Custom(Gv8,46);
 Gv->SetMarkerStyle(22);  Gv2->SetMarkerStyle(23);  Gv3->SetMarkerStyle(25);  Gv4->SetMarkerStyle(3); Gv5->SetMarkerStyle(35); Gv6->SetMarkerStyle(37); Gv7->SetMarkerStyle(24); Gv8->SetMarkerStyle(27);
 Gv->Draw("ALp ");Gv2->Draw(" Lp same"); Gv3->Draw("Lp same "); Gv4->Draw("Lp same "); Gv5->Draw("Lp same "); Gv6->Draw("Lp same "); Gv7->Draw("Lp same "); Gv8->Draw("Lp same ");
//MinMax=FindMinMax(Gv,Gv2,Gv3,Gv4); Diff = (MinMax[1]-MinMax[0])/6.0;
 Gv->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 Gv->GetYaxis()->SetTitle("Vacuum: Imag<#it{O}>");Gv->GetXaxis()->SetTitle("#beta");
 gPad->SetMargin(0.2,0.02,0.12,0.1);Gv->GetYaxis()->SetTitleOffset(1.8);
 sprintf(FileNameOutput,"Test/BareFmunuFmunuVsBeta%s_Nt%d.pdf", QCDTypeS1.c_str(),Nt); 
 c1->SaveAs(FileNameOutput); 
  
 TCanvas *c2 = new TCanvas("c2","Vacuum subtrated value",900,800);
 c2->Divide(2,2,SmallPad,SmallPad);
 c2->SetFillStyle(4000);
 c2->cd(1);
 Custom(stad,2); Custom(stad2,41); Custom(stad3,28); Custom(stad4,8); Custom(stad5,4); Custom(stad6,6); Custom(stad7,1); Custom(stad8,46);
 stad->SetMarkerStyle(22);  stad2->SetMarkerStyle(23);  stad3->SetMarkerStyle(25);  stad4->SetMarkerStyle(3); stad5->SetMarkerStyle(35); stad6->SetMarkerStyle(37); stad7->SetMarkerStyle(24); stad8->SetMarkerStyle(27);
 stad->Draw("AL p");stad2->Draw("Lp same ");stad3->Draw("Lp same "); stad4->Draw("Lp same "); stad5->Draw("Lp same "); stad6->Draw("Lp same "); stad7->Draw("Lp same "); stad8->Draw("Lp same ");
 MinMax=FindMinMax8(stad,stad2,stad3,stad4,stad5,stad6,stad7,stad8); Diff = (MinMax[1]-MinMax[0])/6.0;
 stad->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 stad->GetYaxis()->SetTitle("(Thermal-Vacuum): Real<#it{O}>"); stad->GetXaxis()->SetTitle("#beta");
 gPad->SetMargin(0.2,0.02,0.2,0.1); stad->GetYaxis()->SetTitleOffset(1.6);

 c2->cd(2);
 Custom(Stad,2); Custom(Stad2,41); Custom(Stad3,28); Custom(Stad4,8); Custom(Stad5,4); Custom(Stad6,6); Custom(Stad7,1); Custom(Stad8,46);
 Stad->SetMarkerStyle(22);  Stad2->SetMarkerStyle(23);  Stad3->SetMarkerStyle(25); Stad4->SetMarkerStyle(3);Stad5->SetMarkerStyle(35); Stad6->SetMarkerStyle(37); Stad7->SetMarkerStyle(24); Stad8->SetMarkerStyle(27);
 Stad->Draw("ALp ");Stad2->Draw(" Lp same"); Stad3->Draw("Lp same "); Stad4->Draw("Lp same "); Stad5->Draw("Lp same "); Stad6->Draw("Lp same "); Stad7->Draw(" Lp same"); Stad8->Draw(" Lp same");
//MinMax=FindMinMax(Stad,Stad2,Stad3,Stad4); Diff = (MinMax[1]-MinMax[0])/6.0;
 Stad->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+6*Diff);
 Stad->GetYaxis()->SetTitle("(Thermal-Vacuum): Imag<#it{O}>");Stad->GetXaxis()->SetTitle("#beta");
 gPad->SetMargin(0.16,0.02,0.2,0.1);Stad->GetYaxis()->SetTitleOffset(1.6);
 TLEGEND8(0.2,0.5,0.88,0.88, Stad,"#it{g^{2}a^{4}(F3iF3i - F4iF4i)}","elp", Stad2,"#it{g^{2}a^{4}(F3iF4i + F4iF3i)}","elp", Stad3,"#it{g^{2}a^{5}(F3iDzF3i - F4iDzF4i)}","elp", Stad4,"#it{g^{2}a^{5}(F3iDzF4i + F4iDzF3i)}","elp",Stad5,"#it{g^{2}a^{6}(F3iD2zF3i - F3iD2zF3i)}", "elp", Stad6,"#it{g^{2}a^{6}(F3iD2zF4i + F4iD2zF3i)}","elp", Stad7,"#it{g^{2}a^{8}(F3iD4zF3i - F4iD4zF4i)}", "elp", Stad8,"#it{g^{2}a^{8}(F3iD4zF4i + F4iD4zF3i)}", "elp", 0.05);

 c2->cd(3);
 Custom(ct,2); Custom(ct2,41); Custom(ct3,28); Custom(ct4,8); Custom(ct5,4); Custom(ct6,6); Custom(ct7,1); Custom(ct8,46);
 ct->SetMarkerStyle(22);  ct2->SetMarkerStyle(23);  ct3->SetMarkerStyle(25);  ct4->SetMarkerStyle(3);ct5->SetMarkerStyle(35); ct6->SetMarkerStyle(37); ct7->SetMarkerStyle(24); ct8->SetMarkerStyle(27);
 ct->Draw("AL p");ct2->Draw("Lp same ");ct3->Draw("Lp same "); ct4->Draw("Lp same ");ct5->Draw("Lp same "); ct6->Draw("Lp same "); ct7->Draw("Lp same "); ct8->Draw("Lp same ");
 MinMax=FindMinMax8(ct,ct2,ct3,ct4,ct5,ct6,ct7,ct8); Diff = (MinMax[1]-MinMax[0])/6.0;
 ct->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 ct->GetYaxis()->SetTitle("(Thermal-Vacuum): Real<#it{O} >"); ct->GetXaxis()->SetTitle("T (MeV)");
 gPad->SetMargin(0.2,0.02,0.12,0.1); ct->GetYaxis()->SetTitleOffset(1.6);
 //TLEGEND4(0.2,0.5,0.88,0.88, ct,"<F^{3i}F^{3i}-F^{4i}F^{4i}>/T^{4}","elp", ct2,"<F^{3i}F^{4i}+F^{4i}F^{3i}>/T^{4}","elp", ct3,"<F^{3i}DzF^{3i}-F^{4i}DzF^{4i}>/(T^{4}q^{-})","elp", ct4,"<F^{3i}DzF^{4i}+F^{4i}DzF^{3i}>/(T^{4}q^{-})","elp", 0.05);

 c2->cd(4);
 Custom(CT,2); Custom(CT2,41); Custom(CT3,28); Custom(CT4,8); Custom(CT5,4); Custom(CT6,6); Custom(CT7,1); Custom(CT8,46);
 CT->SetMarkerStyle(22);  CT2->SetMarkerStyle(23);  CT3->SetMarkerStyle(25);  CT4->SetMarkerStyle(3); CT5->SetMarkerStyle(35); CT6->SetMarkerStyle(37); CT7->SetMarkerStyle(24); CT8->SetMarkerStyle(27);
 CT->Draw("ALp ");CT2->Draw(" Lp same"); CT3->Draw("Lp same "); CT4->Draw("Lp same "); CT5->Draw("Lp same "); CT6->Draw("Lp same "); CT7->Draw("Lp same "); CT8->Draw("Lp same ");
 //MinMax=FindMinMax(CT,CT2,CT3,CT4);
 CT->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 CT->GetYaxis()->SetTitle("(Thermal-Vacuum): Imag<#it{O}>");CT->GetXaxis()->SetTitle("T (MeV)");
 gPad->SetMargin(0.2,0.02,0.12,0.1);CT->GetYaxis()->SetTitleOffset(1.8);
 TLEGEND8(0.2,0.5,0.88,0.88, CT,"#it{<F3iF3i - F4iF4i>/T^{4}}","elp", CT2,"#it{<F3iF4i + F4iF3i>/T^{4}}","elp", CT3,"#it{<F3iDzF3i - F4iDzF4i>/(T^{4}q^{-})}","elp", CT4,"#it{<F3iDzF4i + F4iDzF3i>/(T^{4}q^{-})}","elp", CT5,"#it{<F3iD2zF3i - F4iD2zF4i>/(T^{4}(q^{-})^{2})}","elp", CT6,"#it{<F3iD2zF4i + F4iD2zF3i>/(T^{4}(q^{-})^{2})}","elp",CT7,"#it{<F3iD4zF3i - F4iD4zF4i>/(T^{4}(q^{-})^{4})}","elp",CT8,"#it{<F3iD4zF4i + F4iD4zF3i>/(T^{4}(q^{-})^{4})}","elp", 0.05);
 sprintf(FileNameOutput,"FmunuFmunuVsT%s_Nt%d.pdf", QCDTypeS1.c_str(),Nt); 
 c2->SaveAs(FileNameOutput); 
}

TGraphErrors * ConvertToT(TGraphErrors *g , int Nt, int QCDType, double Qminus, string LOorNLO)
{
  cout<<"Started Convert to temperature  Calculation:: for Nt="<<Nt<<", Qminus="<<Qminus<<", LOorNLO="<<LOorNLO<<"\n";
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  double x[N], y[N], ex[N], ey[N];
  double Beta, T, a0, a0T, g0Square, Rbeta, PreFactor=0;
  for(int i=0; i<N; i++)
    {
      Beta = X[i];
      T  = LatticeSpacing(Nt, QCDType, Beta);
      a0 = 1.0/(Nt*T); g0Square = BareCouplingSquare(QCDType, Beta);
      a0T = 1.0/Nt;
      Rbeta = RbetaFmunuSquare(QCDType, Beta);
      x[i] = T;
      ex[i]=0;//cout<<"X[i] = "<<endl;
      if(LOorNLO=="LO")  { PreFactor=0.5/(g0Square*pow(a0T,4)); }
      if(LOorNLO=="NLO") { PreFactor=0.5*sqrt(2)/(g0Square*pow(a0T,4)*a0*Qminus); }
      if(LOorNLO=="NNLO") { PreFactor=-1.0/(g0Square*pow(a0T,4)*a0*a0*Qminus*Qminus); }
      if(LOorNLO=="NNNNLO") { PreFactor=2.0/(g0Square*pow(a0T,4)*pow(a0*Qminus,4)); }

      cout<<LOorNLO<<"\t,Nt="<<Nt<<",Beta="<<Beta<<",T="<<T<<", Qminus="<<Qminus<<endl;
      cout<<"ConvertToT: PreFactor="<<PreFactor<<",Rbeta="<<Rbeta<<endl;
      y[i]  = Y[i]*PreFactor*Rbeta;
      ey[i] = EY[i]*PreFactor*Rbeta;
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

void QhatRe(int Nt, int QCDType, double *Qminus, string File1Output,  string File3Output)
{
  //File1Output : Vacuum subtracted tad-pole corrected Real-part (LO1) LO2 NLO1 (NLO2) (NNLO1) NNLO2 (NNNNLO1) NNNNLO2 w.r.t beta
  cout<<"Started Real-part of q-hat Calculation::\n";
  double PreFactor;// Qminus=20000, 
  double  MuSquare=20000.0,alphas=0.2, NC=3;
  double TOnePlusTtwo = 2.0;

  PreFactor=(8*sqrt(2)*pow(3.1415926,1.0)*alphas/(NC*TOnePlusTtwo));

  char FileName[50000];
  TGraphErrors *g1 = new TGraphErrors(File1Output.c_str(),"%lf  %lf %lf");
  TGraphErrors *g2 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %*lf %*lf  %lf %lf");
  TGraphErrors *g3 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %*lf %*lf %*lf %*lf  %lf %lf");
  TGraphErrors *g4 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %*lf %*lf %*lf %*lf  %*lf %*lf %*lf %*lf %lf %lf");

  const int N=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2=g2->GetN(); double *X2=g2->GetX(), *Y2=g2->GetY(), *Ex2=g2->GetEX(), *Ey2=g2->GetEY();
  int N3=g3->GetN(); double *X3=g3->GetX(), *Y3=g3->GetY(), *Ex3=g3->GetEX(), *Ey3=g3->GetEY();
  int N4=g4->GetN(); double *X4=g4->GetX(), *Y4=g4->GetY(), *Ex4=g4->GetEX(), *Ey4=g4->GetEY();
  double X[N], Y[N], Ex[N], Ey[N]; double YY[N], Eyy[N];
  double T, a0, Beta, a0T, g0Square, Rbeta, Factor; 
  ofstream f; f.open(File3Output.c_str(),ios::app);
  for(int i=0;i<N;i++)
    {
      Beta=X1[i];
      g0Square=BareCouplingSquare(QCDType, Beta);
      T = LatticeSpacing( Nt, QCDType, Beta);
      a0= 1.0/(Nt*T);
      a0T = 1.0/Nt;
      Rbeta = RbetaFmunuSquare(QCDType, Beta);
      Factor = 1/((g0Square*pow(a0T,4.0)));
      Y[i]=PreFactor*Rbeta*Factor*( (0.5*Y1[i]) + (0.5*sqrt(2.0)*(-Y2[i])/(Qminus[0]*a0)) + (-Y3[i]/pow(Qminus[0]*a0,2)) + (2*Y4[i]/pow(Qminus[0]*a0,4.0)) );
      Ey[i]=PreFactor*Rbeta*Factor*sqrt(pow(0.5*Ey1[i],2.0)  + pow(0.5*sqrt(2.0)*Ey2[i]/(Qminus[0]*a0),2.0) + pow(Ey3[i]/pow(Qminus[0]*a0,2.0),2.0)      +  pow( 2*Ey4[i]/pow(Qminus[0]*a0,4.0),2.0)  );
      Ex[i]=0;
      X[i] = T;
      cout<<"(Beta1,Beta4)=("<<X1[i]<<","<<X4[i]<<") \t (Y1,Y4)=("<<Y1[i]<<","<<Y4[i]<<")"<<endl;
      cout<<"\t T="<<T<<"MeV\t QhatReal="<<Y[i]<<", QhatRealErr="<<Ey[i]<<endl; YY[i]=Y[i]/2.0; Eyy[i]=Ey[i]/2.0;
      f<<X[i]<<"\t"<<Y[i]<<"\t"<<Ey[i]<<endl;
    }
  f.close();cout<<endl;
  TCanvas *cc1 = new TCanvas("cc1","QhatReal",600,500);
  CustomGlobal();
  TGraphErrors *G = new TGraphErrors(N,X,Y,Ex,Ey);
  TGraphErrors *GG = new TGraphErrors(N,X,YY,Ex,Eyy);
  Custom(G,2);Custom(GG,2) ;
  G->Draw("ALP ");//GG->Draw("L same");
  G->GetYaxis()->SetRangeUser(-1.0,4.0);
  G->GetYaxis()->SetTitle("Real(#hat{#it{q}}/#it{T}^{3})");
  G->GetXaxis()->SetTitle("#it{T} (MeV)");
  gPad->SetMargin(0.12,0.02,0.2,0.05);G->GetYaxis()->SetTitleOffset(0.8);
sprintf(FileName,"%s.pdf",File3Output.c_str());
cc1->SaveAs(FileName);
}

void QhatIm(int Nt, int QCDType, double *Qminus, string File1Output, string File3Output)
{
  //File1Output : Vacuum subtracted tad-pole corrected Real-part LO1 (LO2) (NLO1) NLO2 NNLO1 (NNLO2) NNNNLO1 (NNNNLO2) w.r.t beta
  cout<<"Started Imag-part of q-hat Calculation::\n";
  double PreFactor;//, Qminus=20000, 
  double MuSquare=20000.0,alphas=0.2, NC=3;
  double TOnePlusTtwo = 2.0;

  PreFactor=(8*sqrt(2)*pow(3.1415926,1.0)*alphas/(NC*TOnePlusTtwo));

  char FileName[50000];
  TGraphErrors *g1 = new TGraphErrors(File1Output.c_str(),"%lf  %*lf %*lf %lf %lf");
  TGraphErrors *g2 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %lf %lf");
  TGraphErrors *g3 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf %lf");
  TGraphErrors *g4 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %*lf %*lf %*lf %*lf  %*lf %*lf %*lf %*lf %*lf %*lf %lf %lf");
                                                        //        LO1        LO2        NLO1     NLO2       NNLO1     NNLO2    NNNNLO1  NNNNLO2

  const int N=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2=g2->GetN(); double *X2=g2->GetX(), *Y2=g2->GetY(), *Ex2=g2->GetEX(), *Ey2=g2->GetEY();
  int N3=g3->GetN(); double *X3=g3->GetX(), *Y3=g3->GetY(), *Ex3=g3->GetEX(), *Ey3=g3->GetEY();
  int N4=g4->GetN(); double *X4=g4->GetX(), *Y4=g4->GetY(), *Ex4=g4->GetEX(), *Ey4=g4->GetEY();
  double X[N], Y[N], Ex[N], Ey[N]; double YY[N], Eyy[N];
  ofstream f; f.open(File3Output.c_str(),ios::app);
  double T, a0, a0T, Beta, g0Square, Rbeta, Factor;
  for(int i=0;i<N;i++)
    {
      Beta=X1[i];
      g0Square=BareCouplingSquare(QCDType, Beta);
      T = LatticeSpacing(Nt, QCDType, Beta);
      a0= 1.0/(Nt*T);
      a0T=1.0/Nt;
      Rbeta = RbetaFmunuSquare(QCDType, Beta);
      Factor = 1/((g0Square*pow(a0T,4.0)));
      Y[i]=PreFactor*Rbeta*Factor*( (0.5*Y1[i]) + (0.5*sqrt(2.0)*Y2[i]/(Qminus[0]*a0)) + (-Y3[i]/pow(Qminus[0]*a0,2)) + (2*Y4[i]/pow(Qminus[0]*a0,4.0)) );
      Ey[i]=PreFactor*Rbeta*Factor*sqrt(pow(0.5*Ey1[i],2.0)  + pow(0.5*sqrt(2.0)*Ey2[i]/(Qminus[0]*a0),2.0) + pow(Ey3[i]/pow(Qminus[0]*a0,2.0),2.0)      +  pow( 2*Ey4[i]/pow(Qminus[0]*a0,4.0),2.0)  );
      Ex[i]=0;
      X[i] = T;
      cout<<"(Beta2,Beta3)=("<<X2[i]<<","<<X3[i]<<") \t (Y2,Y3)=("<<Y2[i]<<","<<Y3[i]<<")"<<endl;      
      cout<<"T="<<T<<"MeV \t, QhatIm="<<Y[i]<<", QhatImErr="<<Ey[i]<<endl; YY[i]=Y[i]/2.0; Eyy[i]=Ey[i]/2.0;
      f<<X[i]<<"\t"<<Y[i]<<"\t"<<Ey[i]<<endl;
    }
  f.close();
  TCanvas *cc2 = new TCanvas("cc2","QhatIm",600,500);
  CustomGlobal();
  TGraphErrors *G = new TGraphErrors(N,X,Y,Ex,Ey);
  TGraphErrors *GG = new TGraphErrors(N,X,YY,Ex,Eyy);
  Custom(G,2);Custom(GG,1) ;
  G->Draw("ALP ");//GG->Draw("L same");
  G->GetYaxis()->SetRangeUser(-1.0,4.0);
  G->GetYaxis()->SetTitle("Im(#hat{#it{q}}/#it{T}^{3})");
  G->GetXaxis()->SetTitle("#it{T} (MeV)");
  gPad->SetMargin(0.12,0.02,0.2,0.05);G->GetYaxis()->SetTitleOffset(0.8);
sprintf(FileName,"%s.pdf",File3Output.c_str());
cc2->SaveAs(FileName);

}

double BareCouplingSquare(int QCDType, double Beta)
{
  double g0Square=0.0;
  if(QCDType>0) {g0Square=10.0/Beta;}
  if(QCDType<=0){g0Square=6.0/Beta;}
  return g0Square;
}
double RbetaFmunuSquare(int QCDType, double beta)
{ double RBETA=1.0;
  TGraph *g = new TGraph("HISQ_RBeta/R_beta_HISQActionHOTQCD.txt","%lf %lf");
  RBETA = g->Eval(beta);
  if(QCDType<=0)  { RBETA=1.0;}
  return RBETA;
}
double LatticeSpacing(int Nt, int QCDType, double beta)
{
  double *Beta;string *BetaS; int FileCount;
  double *Temperature;

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
      double BETA[10] ={6,        6.215,    6.285,    6.423,    6.664,    6.95,     7.15,      7.373,    7.596,     		7.825};
      string BETAS[10]={"6.0000", "6.2150", "6.2850", "6.4230", "6.6640", "6.9500", "7.1500",  "7.3730", "7.5960",  		"7.8250"};
      double MS[10]  ={0.1138,  0.0862,   0.079,  0.067,   0.0514,  0.0386,  0.032,   0.025,   0.0202,  0.0164};
      double T[10]   ={147,      181,      194,    222,     281,     370,     446,     547,     667,     815};
      FileCount=10; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }
  
  if(QCDType>0 && Nt==8)
    {
      double BETA[10] ={6.285,  6.515,   6.575,   6.664,   6.95,    7.28,    7.5,     7.596,   7.825,   8.2};
      string BETAS[10]={"6.2850", "6.5150", "6.5750",  "6.6640",  "6.9500", "7.2800",  "7.5000", "7.5960", "7.8250",	"8.2000"};
      double MS[10]  ={0.079,  0.0604,  0.0564,  0.0514,  0.0386,  0.0284,  0.0222,  0.0202,  0.0164,  0.01167};
      double T[10]   ={146,    182,     193,     211,     277,     377,     459,     500,     611,     843};
      FileCount=10; Beta= new double[FileCount];  BetaS= new string[FileCount]; Temperature = new double[FileCount];
      for(int i=0;  i<FileCount; i++)
        {Beta[i] = BETA[i];  BetaS[i] =BETAS[i]; Temperature[i] = T[i];
        }
    }

  TGraph *g = new TGraph(FileCount,Beta,Temperature);
  return g->Eval(beta);
}
