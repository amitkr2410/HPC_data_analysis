#include <stdlib.h>
#include "StatTool.h"
#include "Custom.h"

//void Mean(TGraph *g, double *Ymean, double *Yerr, int FirstPoint, int SkipN);
void AnalyzeNtNs(int Nt,int Ns,int QCDType,int FmunuType, int WithOrWithoutTadpoleCorrection, double Qminus,  string File0, string File00, string File1,  string File11, string File12,  string File20, string File21, string File22); 
TGraphErrors * ConvertToT(TGraphErrors *g , int Nt, int QCDType, double Qminus, string LOorNLO);
TGraphErrors * ConvertToTWithoutQminus(TGraphErrors *g , int Nt, int QCDType, string LOorNLO);
void QhatRe(int Nt, int QCDType, double *Qminus, string File1Output, string File3Output);
void QhatIm(int Nt, int QCDType, double *Qminus, string File1Output, string File3Output);
double BareCouplingSquare(int QCDType, double Beta);
double StrongCouplingOneLoop(int QCDType, int Nt, double T);
double RbetaFmunuSquare(int QCDType, double beta);
double LatticeSpacing(int Nt, int QCDType, double beta);
double TadpoleFactorTUMQCD(double beta);
void VacuumResultNt32(double beta, int j, int WithOrWithoutTadpoleCorrection, double *MeanAndErrorFit);

void StatisticalAverageAll_RenormalizationU0(int NT)
{
  //gSystem->Exec("./Clean.sh");
  int Nt=NT;
  int Ns=Nt*4;
  int QCDType=0, FmunuType=2, WithOrWithoutTadpoleCorrection=0;
  double Qminus[5]={100000, 20000, 50000, 80000, 100000}; //in MeV
  //Plaquette, tadpolefactor, and bare free energy
  char Output0[10000], Output00[10000];
  //   T+V LO,         T+V NLO,        V LO            V NLO (w.r.t beta)
  char Output1[10000], Output2[10000], Output3[10000], Output4[10000];
  //   Vacuum Sub_LO    Vacuum Sub_NLO (w.r.t beta)
  char Output11[10000], Output12[10000], Output13[10000], Output14[10000];
  //   Vacuum Sub_LO    Vacuum Sub_NLO (Operators dimensionless w.r.t T) 
  char Output20[10000], Output21[10000], Output22[10000];
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
  sprintf(Output20,"Test/%s%sVacuumSubtractedOperatorVsT_RealPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt);
  sprintf(Output21,"Test/%s%sVacuumSubtractedOperatorReNormalizedVsT_RealPartNt%d.txt",QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt);
  sprintf(Output22,"Test/%s%sVacuumSubtractedOperatorReNormalizedVsT_ImagPartNt%d.txt", QCDTypeS.c_str(),FmunuTypeS.c_str(),Nt);
  AnalyzeNtNs(Nt, Ns, QCDType,FmunuType, WithOrWithoutTadpoleCorrection, Qminus[0], Output0, Output00, Output1,  Output11, Output12,  Output20, Output21, Output22);
  
  //DimensionlessOperatorVsT(Nt, Qminus, Output11, Output12,   Output21, Output22);
  sprintf(Output31,"Test/QhatRealVsT_Nt%d_Type%s%s",Nt,QCDTypeS.c_str(),FmunuTypeS.c_str());
  sprintf(Output32,"Test/QhatImVsT_Nt%d_Type%s%s",Nt, QCDTypeS.c_str(),FmunuTypeS.c_str());
  QhatRe(Nt, QCDType, &Qminus[0], Output11,  Output31);
  QhatIm(Nt, QCDType, &Qminus[0], Output11,  Output32);  

}

void AnalyzeNtNs(int Nt, int Ns, int QCDType, int FmunuType, int WithOrWithoutTadpoleCorrection, double Qminus, string OutputFile0, string OutputFile00, string OutputFile1, string OutputFile11, string OutputFile12, string OutputFile20, string OutputFile21, string OutputFile22)
{
  CustomGlobal();
  int FirstRun=500; int SkipN=10; if(QCDType<=0){SkipN=5;}
  int Norm;
  string  NtNs, QCDTypeS, QCDTypeS1, QCDTypeS2; 

  if(QCDType<=0){QCDTypeS="PureGauge";Norm=1.0;}
  if(QCDType>0){QCDTypeS="FullQCD";Norm=1.0;}
  if(FmunuType==1) { QCDTypeS1=QCDTypeS+"SinglePlaquette";         QCDTypeS2="SinglePlaquette_";}
  if(FmunuType==2) { QCDTypeS1=QCDTypeS+"SinglePlaquetteTraceless"; QCDTypeS2="SinglePlaquette_Traceless_";}
  if(FmunuType==3) { QCDTypeS1=QCDTypeS+"Clover";                  QCDTypeS2="Clover_";}
  if(FmunuType==4) { QCDTypeS1=QCDTypeS+"CloverTraceless";         QCDTypeS2="Clover_Traceless_";}
  int FileCount,NtTemp; double *Beta; string *BetaS; double *Temperature; char NtTempS[100], NsTempS[100];  

  if(QCDType<=0 && Nt==4)
    {
      double BETA[9] = {5.6,      5.7,      5.8,      5.9,      6,        6.2,      6.35,     6.5,      6.6};
      string BETAS[9]= {"5.6000", "5.7000", "5.8000", "5.9000", "6.0000", "6.2000", "6.3500", "6.5000", "6.6000"};
      double T[9] =    {208.61,   270.65,   335.85,    405.79,   481.89,   657.89,   816.19,   1002.76,  1145.77};
      //Old: 191.423 
      //double T[9]   =  { 253.235,  283.295,  316.964,  354.68,   396.931,  497.306,  589.089,  697.975,  781.629};
      //Karsch: 137.753,
      //double T[9]   =  220.99,   270.21,   333.316,  405.588,  485.75,   667.11,   824.04,   1010.12,  1151.89};
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
  double Ymean[3][4], Yerr[3][4];
  double mean[3], err[3]; double MeanAndErrorFit[2]={0,0};
  for(int k=0; k<2;k++)                   //Thermal or Vacuum
    {
      for(int i=0; i<FileCount; i++)      // different beta
	{ 
	  NtTemp=Nt;
	  if(k==1){NtTemp=Ns;}//NtNs= "Nt16_Ns16_"                                                                                          
	  sprintf(NtTempS,"%d",NtTemp);
	  sprintf(NsTempS,"%d",Ns);

	  for(int j=0; j<3;j++)           // Tadpole0 or Tadpole2 or Tadpole4 D4z derivatives terms
	    { if(j==0) {LOorNLO[j]="LO_";     TadpoleOrder="Tadpole0_";}	    
	      if(j==1) {LOorNLO[j]="NNLO_";   TadpoleOrder="Tadpole2_";}
	      if(j==2) {LOorNLO[j]="NNNNLO_"; TadpoleOrder="Tadpole4_";}	     
	      File[j]=std::string("DataTraceFmunu") + TadpoleOrder + QCDTypeS2 + "Nt" + NtTempS + "_Ns"+ NsTempS  +"_" + "Beta"+ BetaS[i] + ".txt";
	      sprintf(FileName[j],"RawData%sRenormalizationU0/Nt%d_Ns%d/%s",QCDTypeS1.c_str(), NtTemp, Ns, File[j].c_str());
	      //RawDataFullQCDCloverTraceless/Nt4_Ns16/DataTraceFmunuTadpole0_Clover_Traceless_Nt4_Ns16_Beta5.9000.txt 
              File2=std::string("DataPloopNt") + NtTempS + "_Ns"+ NsTempS  + "_" + "Beta"+ BetaS[i] + ".txt";
              sprintf(FileName2,"RawData%s_Plaquette_PLoopRenormalizationU0/Nt%d_Ns%d/%s",QCDTypeS.c_str(), NtTemp, Ns, File2.c_str());
              //RawDataFullQCD_Plaquette_PLoop/Nt4_Ns16/DataPloopNt4_Ns16_Beta5.9000.txt
	      for(int m=0;m<4;m++){ Ymean[j][m]=0; Yerr[j][m]=0;}
	      mean[j]=0; err[j]=0;
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
	  
	  if(WithOrWithoutTadpoleCorrection==0){TadpoleFactorU0=1.0;};
	  
	 cout<<"\n Tr(F3i..F3i).Re (Nt,Ns)="<<NtTemp<<","<<Ns<<": Beta="<<Beta[i]<<endl;
	 GM0 = new TGraph(FileName[0],"%lf %lf");
	 GM1 = new TGraph(FileName[1],"%lf %lf");
	 GM2 = new TGraph(FileName[2],"%lf %lf");
	 JackknifeMeanD4z(GM0,GM1,GM2, TadpoleFactorU0, mean, err, FirstRun, SkipN);
	 for(int j=0; j<3;j++) { Ymean[j][0]=mean[j]; Yerr[j][0]=err[j];}
	 
	 cout<<" Tr(F3i..F3i).Im"<<endl;
	 GM0 = new TGraph(FileName[0],"%lf %*lf %lf");
	 GM1 = new TGraph(FileName[1],"%lf %*lf %lf");
	 GM2 = new TGraph(FileName[2],"%lf %*lf %lf");
	 JackknifeMeanD4z(GM0,GM1,GM2, TadpoleFactorU0, mean, err, FirstRun, SkipN);
	 for(int j=0; j<3;j++) { Ymean[j][1]=mean[j]; Yerr[j][1]=err[j];}
	 
	 cout<<" Tr(F3i..F4i).Re "<<endl;
	 GM0 = new TGraph(FileName[0],"%lf %*lf %*lf     %*lf %*lf     %lf");
	 GM1 = new TGraph(FileName[1],"%lf %*lf %*lf     %*lf %*lf     %lf");
	 GM2 = new TGraph(FileName[2],"%lf %*lf %*lf     %*lf %*lf     %lf");
	 JackknifeMeanD4z(GM0,GM1,GM2, TadpoleFactorU0, mean, err, FirstRun, SkipN);
	 for(int j=0; j<3;j++) { Ymean[j][2]=mean[j]; Yerr[j][2]=err[j];}
	 
	 cout<<" Tr(F3i..F4i).Im"<<endl;
	 GM0 = new TGraph(FileName[0],"%lf %*lf %*lf     %*lf %*lf     %*lf %lf ");
	 GM1 = new TGraph(FileName[1],"%lf %*lf %*lf     %*lf %*lf     %*lf %lf ");
	 GM2 = new TGraph(FileName[2],"%lf %*lf %*lf     %*lf %*lf     %*lf %lf ");
	 JackknifeMeanD4z(GM0,GM1,GM2, TadpoleFactorU0, mean, err, FirstRun, SkipN);
         for(int j=0; j<3;j++) { Ymean[j][3]=mean[j]; Yerr[j][3]=err[j];}	 


	 for(int j=0; j<3;j++)
	   { sprintf(FileNameOutput,"%s_%sNt%d_Ns%d.txt",OutputFile1.c_str(), LOorNLO[j].c_str(),NtTemp, Ns );  //LO  T+V
	     if(k==0 && j==0){ sprintf(FileNameOutputLO,    "%s",FileNameOutput);}
	     if(k==0 && j==1){ sprintf(FileNameOutputNNLO,   "%s",FileNameOutput);}
	     if(k==0 && j==2){ sprintf(FileNameOutputNNNNLO, "%s",FileNameOutput);}
	     if(k==1 && j==0){ sprintf(FileNameOutputVacLO,  "%s",FileNameOutput);}
	     if(k==1 && j==1){ sprintf(FileNameOutputVacNNLO,  "%s",FileNameOutput);}
	     if(k==1 && j==2){ sprintf(FileNameOutputVacNNNNLO,  "%s",FileNameOutput);}

	     fout.open(FileNameOutput, ios::app); double Small=pow(10,-22);
	     if( QCDType>0  && NtTemp==32 && k==1 ) //NLO, Vacuum
	       {  VacuumResultNt32(Beta[i], j, WithOrWithoutTadpoleCorrection, MeanAndErrorFit);
		 fout<<Beta[i]<<"\t"<< MeanAndErrorFit[0]<<"\t\t"<< MeanAndErrorFit[1]<<"\t"<<Ymean[j][1]<<"\t\t"<<Yerr[j][1]<<"\t"
		     <<Ymean[j][2]<<"\t\t"<<Yerr[j][2]<<"\t"<<Ymean[j][3]<<"\t\t"<<Yerr[j][3]<<endl;
		 cout<<Beta[i]<<"\t"<<Ymean[j][0]<<"\t\t"<<Yerr[j][0]<<"\t"<<Ymean[j][1]<<"\t\t"<<Yerr[j][1]<<"\t"
		     <<Ymean[j][2]<<"\t\t"<<Yerr[j][2]<<"\t"<<Ymean[j][3]<<"\t\t"<<Yerr[j][3]<<endl;}
	     else {fout<<Beta[i]<<"\t"<<Ymean[j][0]<<"\t\t"<<Yerr[j][0]<<"\t"<<Ymean[j][1]<<"\t\t"<<Yerr[j][1]<<"\t"
		       <<Ymean[j][2]<<"\t\t"<<Yerr[j][2]<<"\t"<<Ymean[j][3]<<"\t\t"<<Yerr[j][3]<<endl;
	           cout<<Beta[i]<<"\t"<<Ymean[j][0]<<"\t\t"<<Yerr[j][0]<<"\t"<<Ymean[j][1]<<"\t\t"<<Yerr[j][1]<<"\t"
		       <<Ymean[j][2]<<"\t\t"<<Yerr[j][2]<<"\t"<<Ymean[j][3]<<"\t\t"<<Yerr[j][3]<<endl;
	           }
	     fout.close();
	   }   //end of LO or  NNLO or NNNNLO Fmunu loop 0-1
	}	//end of different beta files loop 0-9
    }           //end of thermal or vacuum loop 0-1
  
  
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
  c0->SaveAs(FileNameOutput);
  
 //(Nt,Ns) //Re-part   LO1, LO2, NNLO1, NNLO2, NNNNLO1, NNNNLO2
 TGraphErrors *g  = new TGraphErrors(FileNameOutputLO,    "%lf  %lf  %lf");  
 TGraphErrors *g2 = new TGraphErrors(FileNameOutputLO,    "%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *g3 = new TGraphErrors(FileNameOutputNNLO,  "%lf  %lf  %lf");
 TGraphErrors *g4 = new TGraphErrors(FileNameOutputNNLO,  "%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *g5 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %lf  %lf");
 TGraphErrors *g6 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
           //Im-part   LO1, LO2, NNLO1, NNLO2, NNNNLO1, NNNNLO2 
 TGraphErrors *G  = new TGraphErrors(FileNameOutputLO,    "%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G2 = new TGraphErrors(FileNameOutputLO,    "%lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G3 = new TGraphErrors(FileNameOutputNNLO,  "%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G4 = new TGraphErrors(FileNameOutputNNLO,  "%lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G5 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *G6 = new TGraphErrors(FileNameOutputNNNNLO,"%lf  %*lf %*lf %*lf %*lf  %*lf %*lf %lf  %lf");
 //(Ns,Ns) //Re-part   LO1, LO2, NNLO1, NNLO2, NNNNLO1, NNNNLO2 
 TGraphErrors *gv  = new TGraphErrors(FileNameOutputVacLO,    "%lf  %lf  %lf");
 TGraphErrors *gv2 = new TGraphErrors(FileNameOutputVacLO,    "%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *gv3 = new TGraphErrors(FileNameOutputVacNNLO,  "%lf  %lf  %lf");
 TGraphErrors *gv4 = new TGraphErrors(FileNameOutputVacNNLO,  "%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *gv5 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %lf  %lf");
 TGraphErrors *gv6 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %*lf %*lf %*lf %*lf  %lf  %lf");
           //Im-part   LO1, LO2, NNLO1, NNLO2, NNLO1, NNLO2
 TGraphErrors *Gv  = new TGraphErrors(FileNameOutputVacLO,"%lf  %*lf %*lf %lf  %lf");
 TGraphErrors *Gv2 = new TGraphErrors(FileNameOutputVacLO,"%lf  %*lf %*lf %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *Gv3 = new TGraphErrors(FileNameOutputVacNNLO,"%lf  %*lf  %*lf %lf %lf");
 TGraphErrors *Gv4 = new TGraphErrors(FileNameOutputVacNNLO,"%lf  %*lf %*lf %*lf %*lf %*lf %*lf  %lf  %lf");
 TGraphErrors *Gv5 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %*lf  %*lf %lf %lf");
 TGraphErrors *Gv6 = new TGraphErrors(FileNameOutputVacNNNNLO,"%lf  %*lf %*lf %*lf %*lf %*lf %*lf  %lf  %lf");


 TGraphErrors *stad, *stad2, *stad3, *stad4, *stad5, *stad6, *stad7, *stad8; //Real with tad-pole factor corrected
 TGraphErrors *Stad, *Stad2, *Stad3, *Stad4, *Stad5, *Stad6, *Stad7, *Stad8; //Imag with tad-pole factor corrected

 stad=SubtractTGraphErrors(g,gv); stad2=SubtractTGraphErrors(g2,gv2); stad3=SubtractTGraphErrors(g3,gv3);stad4=SubtractTGraphErrors(g4,gv4); stad5=SubtractTGraphErrors(g5,gv5); stad6=SubtractTGraphErrors(g6,gv6); //Real-part LO1 LO2  NNLO1 NNLO2  NNNNLO1 NNNNLO2
 Stad=SubtractTGraphErrors(G,Gv); Stad2=SubtractTGraphErrors(G2,Gv2); Stad3=SubtractTGraphErrors(G3,Gv3);Stad4=SubtractTGraphErrors(G4,Gv4); Stad5=SubtractTGraphErrors(G5,Gv5); Stad6=SubtractTGraphErrors(G6,Gv6); //Im-part   LO1 LO2  NNLO1 NNLO2  NNNNLO1 NNNNLO2
 WriteToFile(stad, stad2, stad3, stad4, stad5, stad6,  OutputFile11);//Vacuum subtracted Real-part LO1 LO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t beta
 WriteToFile(Stad, Stad2, Stad3, Stad4, Stad5, Stad6,  OutputFile12);//Vacuum subtracted Im-part LO1 LO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t beta
 
 //convert to temperature dependent form and write to a file
 TGraphErrors *Gct, *Gct2, *Gct3, *Gct4, *Gct5, *Gct6; //Real, without R_beta and Qminus
 TGraphErrors *ct, *ct2, *ct3, *ct4, *ct5, *ct6, *ct7, *ct8; //Real
 TGraphErrors *CT, *CT2, *CT3, *CT4, *CT5, *CT6, *CT7, *CT8; //Imag

 Gct  = ConvertToTWithoutQminus(stad  , Nt, QCDType,  "LO");  //Real-part LO1, without R_beta and Qminus
 Gct2 = ConvertToTWithoutQminus(stad2 , Nt, QCDType,  "LO");  //Real-part LO2
 Gct3 = ConvertToTWithoutQminus(stad3 , Nt, QCDType,  "NNLO"); //Real-part NNLO1
 Gct4 = ConvertToTWithoutQminus(stad4 , Nt, QCDType,  "NNLO"); //Real-part NNLO2
 Gct5 = ConvertToTWithoutQminus(stad5 , Nt, QCDType,  "NNNNLO"); //Real-part NNNNLO1
 Gct6 = ConvertToTWithoutQminus(stad6 , Nt, QCDType,  "NNNNLO"); //Real-part NNNNLO2

 ct  = ConvertToT(stad  , Nt, QCDType, Qminus, "LO");  //Real-part LO1, with R_beta and Qminus
 ct2 = ConvertToT(stad2 , Nt, QCDType, Qminus, "LO");  //Real-part LO2
 ct3 = ConvertToT(stad3 , Nt, QCDType, Qminus, "NNLO"); //Real-part NNLO1
 ct4 = ConvertToT(stad4 , Nt, QCDType, Qminus, "NNLO"); //Real-part NNLO2
 ct5 = ConvertToT(stad5 , Nt, QCDType, Qminus, "NNNNLO"); //Real-part NNNNLO1
 ct6 = ConvertToT(stad6 , Nt, QCDType, Qminus, "NNNNLO"); //Real-part NNNNLO2
 CT  = ConvertToT(Stad  , Nt, QCDType, Qminus, "LO");  //Im-part LO1
 CT2 = ConvertToT(Stad2 , Nt, QCDType, Qminus, "LO");  //Im-part LO2
 CT3 = ConvertToT(Stad3 , Nt, QCDType, Qminus, "NNLO"); //Im-part NNLO1
 CT4 = ConvertToT(Stad4 , Nt, QCDType, Qminus, "NNLO"); //Im-part NNLO2
 CT5 = ConvertToT(Stad5 , Nt, QCDType, Qminus, "NNNNLO"); //Im-part NNNNLO1
 CT6 = ConvertToT(Stad6 , Nt, QCDType, Qminus, "NNNNLO"); //Im-part NNNNLO2

 WriteToFile(Gct, Gct2, Gct3, Gct4, Gct5, Gct6, OutputFile20);  //Vacuum subtracted Real-part LO1 LO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 without R_beta and Qminus w.r.t T
 WriteToFile(ct, ct2, ct3, ct4, ct5, ct6, OutputFile21);  //Vacuum subtracted Real-part LO1 LO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t T
 WriteToFile(CT, CT2, CT3, CT4, CT5, CT6, OutputFile22);  //Vacuum subtracted   Im-part LO1 LO2  NNLO1 NNLO2 NNNNLO1 NNNNLO2 w.r.t T 

 NtNs="Nt"+std::to_string(Nt)+"Ns"+std::to_string(Ns);
 TCanvas *c1 = new TCanvas(NtNs.c_str(),NtNs.c_str(),900,800);
 double SmallPad=pow(10,-54);
 c1->Divide(2,2,SmallPad,SmallPad);
 c1->SetFillStyle(4000);
 c1->cd(1);
 Custom(g,2); Custom(g2,41); Custom(g3,4); Custom(g4,46); Custom(g5,1); Custom(g6,6);
 g->SetMarkerStyle(22); g2->SetMarkerStyle(23); g3->SetMarkerStyle(24);  g4->SetMarkerStyle(27);  g5->SetMarkerStyle(35);  g6->SetMarkerStyle(37); 
 g->Draw("AL p");g2->Draw("p same ");g3->Draw("p same "); g4->Draw("p same "); g5->Draw("p same "); g6->Draw("p same"); 
 MinMax=FindMinMax6(g,g2,g3,g4,g5,g6); Diff = (MinMax[1]-MinMax[0])/6.0;
 g->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 g->GetYaxis()->SetTitle("T+V: Real<#it{O}>"); g->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.2,0.02,0.2,0.1); g->GetYaxis()->SetTitleOffset(1.6);

 c1->cd(2);
 Custom(G,2); Custom(G2,41); Custom(G3,4); Custom(G4,46); Custom(G5,1); Custom(G6,6);
 G->SetMarkerStyle(22);  G2->SetMarkerStyle(23);  G3->SetMarkerStyle(24); G4->SetMarkerStyle(27); G5->SetMarkerStyle(35); G6->SetMarkerStyle(37); 
 G->Draw("ALp "); G2->Draw(" Lp same"); G3->Draw("Lp same "); G4->Draw("Lp same "); G5->Draw("Lp same "); G6->Draw("Lp same "); 
 MinMax=FindMinMax6(G,G2,G3,G4,G5,G6);Diff = (MinMax[1]-MinMax[0])/6.0;
 G->GetYaxis()->SetRangeUser(MinMax[0]-Diff, MinMax[1]+6*Diff);
 G->GetYaxis()->SetTitle("T+V: Imag<#it{O}>");G->GetXaxis()->SetTitle("#beta_{0}"); 
 gPad->SetMargin(0.16,0.02,0.2,0.1);G->GetYaxis()->SetTitleOffset(1.5);
 //TLEGEND4(0.25,0.5,0.88,0.88, G,"#it{g^{2}a^{4}(F3iF3i - F4iF4i)}","elp", G2,"#it{g^{2}a^{4}(F3iF4i + F4iF3i)}","elp", G3,"#it{g^{2}a^{5}(F3iDzF3i - F4iDzF4i)}","elp", G4,"#it{g^{2}a^{5}(F3iDzF4i + F4iDzF3i)}","elp", 0.05);
 TLEGEND6(0.25,0.5,0.88,0.88, G,"#it{g^{2}a^{4}(F3iF3i - F4iF4i)}","elp", G2,"#it{g^{2}a^{4}(F3iF4i + F4iF3i)}","elp", G3,"#it{g^{2}a^{5}(F3iD2z F3i - F4iD2zF4i)}","elp", G4,"#it{g^{2}a^{5}(F3iD2zF4i + F4iD2zF3i)}","elp", G5,"#it{g^{2}a^{6}(F3iD4zF3i - F4iD4zF4i)}","elp", G6,"#it{g^{2}a^{6}(F3iD4zF4i + F4iD4zF3i)}","elp", 0.05);
 c1->cd(3);
 Custom(gv,2); Custom(gv2,41); Custom(gv3,4); Custom(gv4,46); Custom(gv5,1); Custom(gv6,6); 
 gv->SetMarkerStyle(22);  gv2->SetMarkerStyle(23);  gv3->SetMarkerStyle(24);  gv4->SetMarkerStyle(27); gv5->SetMarkerStyle(35); gv6->SetMarkerStyle(37); 
 gv->Draw("AL p");gv2->Draw("Lp same ");gv3->Draw("Lp same "); gv4->Draw("Lp same "); gv5->Draw("Lp same "); gv6->Draw("Lp same "); 
 MinMax=FindMinMax6(gv,gv2,gv3,gv4,gv5,gv6); Diff = (MinMax[1]-MinMax[0])/6.0;
 gv->GetYaxis()->SetRangeUser(MinMax[0] -Diff, MinMax[1]+Diff);
 gv->GetYaxis()->SetTitle("Vacuum: Real<#it{O} >"); gv->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.2,0.02,0.12,0.1); gv->GetYaxis()->SetTitleOffset(1.6);

 c1->cd(4);
 Custom(Gv,2); Custom(Gv2,41); Custom(Gv3,4); Custom(Gv4,46); Custom(Gv5,1); Custom(Gv6,6);
 Gv->SetMarkerStyle(22);  Gv2->SetMarkerStyle(23);  Gv3->SetMarkerStyle(24);  Gv4->SetMarkerStyle(27); Gv5->SetMarkerStyle(35); Gv6->SetMarkerStyle(37); 
 Gv->Draw("ALp ");Gv2->Draw(" Lp same"); Gv3->Draw("Lp same "); Gv4->Draw("Lp same "); Gv5->Draw("Lp same "); Gv6->Draw("Lp same ");
//MinMax=FindMinMax(Gv,Gv2,Gv3,Gv4); Diff = (MinMax[1]-MinMax[0])/6.0;
 Gv->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 Gv->GetYaxis()->SetTitle("Vacuum: Imag<#it{O}>");Gv->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.2,0.02,0.12,0.1);Gv->GetYaxis()->SetTitleOffset(1.8);
 sprintf(FileNameOutput,"Test/BareFmunuFmunuVsBeta%s_Nt%d.pdf", QCDTypeS1.c_str(),Nt); 
 c1->SaveAs(FileNameOutput); 
 
 TCanvas *c2 = new TCanvas("c2","Vacuum subtrated value",900,800);
 c2->Divide(2,2,SmallPad,SmallPad);
 c2->SetFillStyle(4000);
 c2->cd(1);
 Custom(stad,2); Custom(stad2,41); Custom(stad3,4); Custom(stad4,46); Custom(stad5,1); Custom(stad6,6); 
 stad->SetMarkerStyle(22);  stad2->SetMarkerStyle(23);  stad3->SetMarkerStyle(24);  stad4->SetMarkerStyle(27); stad5->SetMarkerStyle(35); stad6->SetMarkerStyle(37); 
 stad->Draw("AL p");//stad2->Draw("Lp same ");
 stad3->Draw("Lp same "); //stad4->Draw("Lp same ");
 stad5->Draw("Lp same "); //stad6->Draw("Lp same ");
 MinMax=FindMinMax6(stad,stad2,stad3,stad4,stad5,stad6); Diff = (MinMax[1]-MinMax[0])/6.0;
 stad->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 stad->GetYaxis()->SetTitle("(Thermal-Vacuum): Real<#it{O}>"); stad->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.2,0.02,0.2,0.1); stad->GetYaxis()->SetTitleOffset(1.6);
 
 c2->cd(2);
 Custom(Stad,2); Custom(Stad2,41); Custom(Stad3,4); Custom(Stad4,46); Custom(Stad5,1); Custom(Stad6,6); 
 Stad->SetMarkerStyle(22);  Stad2->SetMarkerStyle(23); Stad3->SetMarkerStyle(24); Stad4->SetMarkerStyle(27); Stad5->SetMarkerStyle(35); Stad6->SetMarkerStyle(37); 
 Stad->Draw("ALp ");Stad2->Draw(" Lp same");  Stad3->Draw("Lp same "); Stad4->Draw("Lp same ");  Stad5->Draw("Lp same "); Stad6->Draw("Lp same "); 
//MinMax=FindMinMax(Stad,Stad2,Stad3,Stad4); Diff = (MinMax[1]-MinMax[0])/6.0;
 Stad->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+6*Diff);
 Stad->GetYaxis()->SetTitle("(Thermal-Vacuum): Imag<#it{O}>");Stad->GetXaxis()->SetTitle("#beta_{0}");
 gPad->SetMargin(0.16,0.02,0.2,0.1);Stad->GetYaxis()->SetTitleOffset(1.6);
 TLEGEND6(0.2,0.5,0.88,0.88, Stad,"#it{g^{2}a^{4}(F3iF3i - F4iF4i)}","elp", Stad2,"#it{g^{2}a^{4}(F3iF4i + F4iF3i)}","elp", Stad3,"#it{g^{2}a^{6}(F3iD2zF3i - F3iD2zF3i)}", "elp", Stad4,"#it{g^{2}a^{6}(F3iD2zF4i + F4iD2zF3i)}","elp", Stad5,"#it{g^{2}a^{8}(F3iD4zF3i - F4iD4zF4i)}", "elp", Stad6,"#it{g^{2}a^{8}(F3iD4zF4i + F4iD4zF3i)}", "elp", 0.05);
 
 c2->cd(3);
 Custom(ct,2); Custom(ct2,41); Custom(ct3,4); Custom(ct4,46); Custom(ct5,1); Custom(ct6,6); 
 ct->SetMarkerStyle(22);  ct2->SetMarkerStyle(23);  ct3->SetMarkerStyle(24);  ct4->SetMarkerStyle(27);ct5->SetMarkerStyle(35); ct6->SetMarkerStyle(37); 
 ct->Draw("AL p");//ct2->Draw("Lp same ");
 ct3->Draw("Lp same "); //ct4->Draw("Lp same ");
 ct5->Draw("Lp same "); //ct6->Draw("Lp same "); 
 MinMax=FindMinMax6(ct,ct2,ct3,ct4,ct5,ct6); Diff = (MinMax[1]-MinMax[0])/6.0;
 ct->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 ct->GetYaxis()->SetTitle("(Thermal-Vacuum): Real<#it{O} >"); ct->GetXaxis()->SetTitle("T (MeV)");
 gPad->SetMargin(0.2,0.02,0.12,0.1); ct->GetYaxis()->SetTitleOffset(1.6);
 //TLEGEND4(0.2,0.5,0.88,0.88, ct,"<F^{3i}F^{3i}-F^{4i}F^{4i}>/T^{4}","elp", ct2,"<F^{3i}F^{4i}+F^{4i}F^{3i}>/T^{4}","elp", ct3,"<F^{3i}DzF^{3i}-F^{4i}DzF^{4i}>/(T^{4}q^{-})","elp", ct4,"<F^{3i}DzF^{4i}+F^{4i}DzF^{3i}>/(T^{4}q^{-})","elp", 0.05);
 ct->GetYaxis()->SetRangeUser(MinMax[1]*pow(10,-12),MinMax[1]*8); gPad->SetLogy();

 c2->cd(4);
 Custom(CT,2); Custom(CT2,41); Custom(CT3,4); Custom(CT4,46); Custom(CT5,1); Custom(CT6,6); 
 CT->SetMarkerStyle(22);  CT2->SetMarkerStyle(23);  CT3->SetMarkerStyle(24);  CT4->SetMarkerStyle(27); CT5->SetMarkerStyle(35); CT6->SetMarkerStyle(37); 
 CT->Draw("ALp ");CT2->Draw(" Lp same");  CT3->Draw("Lp same "); CT4->Draw("Lp same "); CT5->Draw("Lp same "); CT6->Draw("Lp same "); 
 //MinMax=FindMinMax(CT,CT2,CT3,CT4);
 CT->GetYaxis()->SetRangeUser(MinMax[0]-Diff,MinMax[1]+Diff);
 CT->GetYaxis()->SetTitle("(Thermal-Vacuum): Imag<#it{O}>");CT->GetXaxis()->SetTitle("T (MeV)");
 gPad->SetMargin(0.2,0.02,0.12,0.1);CT->GetYaxis()->SetTitleOffset(1.8);
 TLEGEND6(0.2,0.5,0.88,0.88, CT,"#it{<F3iF3i - F4iF4i>/T^{4}}","elp", CT2,"#it{<F3iF4i + F4iF3i>/T^{4}}","elp", CT3,"#it{<F3iD2zF3i - F4iD2zF4i>/(T^{4}(q^{-})^{2})}","elp", CT4,"#it{<F3iD2zF4i + F4iD2zF3i>/(T^{4}(q^{-})^{2})}","elp",CT5,"#it{<F3iD4zF3i - F4iD4zF4i>/(T^{4}(q^{-})^{4})}","elp",CT6,"#it{<F3iD4zF4i + F4iD4zF3i>/(T^{4}(q^{-})^{4})}","elp", 0.05);
 sprintf(FileNameOutput,"FmunuFmunuVsT%s_Nt%d.pdf", QCDTypeS1.c_str(),Nt); 
 c2->SaveAs(FileNameOutput); 
 
}

TGraphErrors * ConvertToT(TGraphErrors *g , int Nt, int QCDType, double Qminus, string LOorNLO)
{
  cout<<"Started Convert to temperature  Calculation:: for Nt="<<Nt<<", Qminus="<<Qminus<<", LOorNLO="<<LOorNLO<<"\n";
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  double x[N], y[N], ex[N], ey[N];
  double Beta, T, a0, a0T, g0Square, Rbeta=1.0, PreFactor=0;
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

      y[i]  = Y[i]*PreFactor*Rbeta;
      ey[i] = EY[i]*PreFactor*Rbeta;
      cout<<LOorNLO<<"\t,Nt="<<Nt<<",Beta="<<Beta<<",T="<<T<<", Qminus="<<Qminus<<endl;
      cout<<"ConvertToT: PreFactor="<<PreFactor<<",Rbeta="<<Rbeta<<", Y="<<y[i]<<", ey"<<ey[i]<<endl;
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

TGraphErrors * ConvertToTWithoutQminus(TGraphErrors *g , int Nt, int QCDType, string LOorNLO)
{
  cout<<"Started Convert to temperature without Qminus :: for Nt="<<Nt<<", LOorNLO="<<LOorNLO<<"\n";
  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();
  double x[N], y[N], ex[N], ey[N];
  double Beta, T, a0, a0T, g0Square, PreFactor=0;
  for(int i=0; i<N; i++)
    {
      Beta = X[i];
      T  = LatticeSpacing(Nt, QCDType, Beta);
      a0 = 1.0/(Nt*T); g0Square = BareCouplingSquare(QCDType, Beta);
      a0T = 1.0/Nt;
      x[i] = T;
      ex[i]=0;//cout<<"X[i] = "<<endl;
      if(LOorNLO=="LO")  { PreFactor=1.0/(g0Square*pow(a0T,4)); }
      if(LOorNLO=="NLO") { PreFactor=1.0/(g0Square*pow(a0T,5)); }
      if(LOorNLO=="NNLO") { PreFactor=1.0/(g0Square*pow(a0T,6)); }
      if(LOorNLO=="NNNNLO") { PreFactor=1.0/(g0Square*pow(a0T,8)); }

      y[i]  = Y[i]*PreFactor;
      ey[i] = EY[i]*PreFactor;
      cout<<LOorNLO<<"\t,Nt="<<Nt<<",Beta="<<Beta<<",T="<<T<<endl;
      cout<<"ConvertToT: PreFactor="<<PreFactor<<", Y="<<y[i]<<", ey"<<ey[i]<<endl;
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}
void QhatRe(int Nt, int QCDType, double *Qminus, string File1Output,  string File3Output)
{
  //File1Output : Vacuum subtracted tad-pole corrected Real-part (LO1) LO2 NLO1 (NLO2) (NNLO1) NNLO2 (NNNNLO1) NNNNLO2 w.r.t beta
  cout<<"Started Real-part of q-hat Calculation::\n";
  double PreFactor;// Qminus=20000, 
  double  MuSquare=20000.0,alphas=0.0, NC=3;
  double TOnePlusTtwo = 1.0;

  PreFactor=(8*sqrt(2)*pow(3.1415926,1.0)/(NC*TOnePlusTtwo));

  char FileName[50000];
  TGraphErrors *g1 = new TGraphErrors(File1Output.c_str(),"%lf  %lf %lf");
  TGraphErrors *g2 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %lf %lf");
  TGraphErrors *g3 = new TGraphErrors(File1Output.c_str(),"%lf %*lf %*lf  %*lf %*lf  %*lf %*lf %*lf %*lf  %lf %lf");
                                   //Vacuum subtracted Real-part LO1      LO2         NNLO1     NNLO2     NNNNLO1    NNNNLO2 w.r.t beta
  const int N=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2=g2->GetN(); double *X2=g2->GetX(), *Y2=g2->GetY(), *Ex2=g2->GetEX(), *Ey2=g2->GetEY();
  int N3=g3->GetN(); double *X3=g3->GetX(), *Y3=g3->GetY(), *Ex3=g3->GetEX(), *Ey3=g3->GetEY();
  double X[N], Y[N], Ex[N], Ey[N]; double YY[N], Eyy[N];
  double T, a0, Beta, a0T, g0Square, Rbeta, Factor; 
  ofstream f; f.open(File3Output.c_str(),ios::app);
  for(int i=0;i<N;i++)
    {
      Beta=X1[i];
      g0Square=BareCouplingSquare(QCDType, Beta);
      T   = LatticeSpacing( Nt, QCDType, Beta);
      a0  = 1.0/(Nt*T);
      a0T = 1.0/Nt;
      alphas = StrongCouplingOneLoop(QCDType, Nt, T);
      Rbeta = RbetaFmunuSquare(QCDType, Beta);
      Factor = 1/((g0Square*pow(a0T,4.0)));
      Y[i]=PreFactor*alphas*Rbeta*Factor*( (0.5*Y1[i])  + (-Y2[i]/pow(Qminus[0]*a0,2)) + (2*Y3[i]/pow(Qminus[0]*a0,4.0)) );
      Ey[i]=PreFactor*alphas*Rbeta*Factor*sqrt(pow(0.5*Ey1[i],2.0)  + pow(Ey2[i]/pow(Qminus[0]*a0,2.0),2.0)      +  pow( 2*Ey3[i]/pow(Qminus[0]*a0,4.0),2.0)  );
      Ex[i]=0;
      X[i] = T;
      cout<<"(Beta1,Beta2,Beta3)=("<<X1[i]<<","<<X2[i]<<","<<X3[i]<<") \t (Y1,Y2,Y3)=("<<Y1[i]<<","<<Y2[i]<<","<<Y3[i]<<")"<<endl;
      cout<<"Rbeta="<<Rbeta<<", alphas="<<alphas<<endl;
      cout<<"\t T="<<T<<"MeV\t QhatReal="<<Y[i]<<", QhatRealErr="<<Ey[i]<<"\n"<<endl; YY[i]=Y[i]/2.0; Eyy[i]=Ey[i]/2.0;
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
  G->GetXaxis()->SetRangeUser(100,850);
  G->GetYaxis()->SetTitle("Real(#hat{#it{q}}/#it{T}^{3})");
  G->GetXaxis()->SetTitle("#it{T} (MeV)");
  gPad->SetMargin(0.12,0.02,0.2,0.05);G->GetYaxis()->SetTitleOffset(0.8);
sprintf(FileName,"%s.pdf",File3Output.c_str());
cc1->SaveAs(FileName);
}

void QhatIm(int Nt, int QCDType, double *Qminus, string File1Output, string File3Output)
{
  //File1Output : Vacuum subtracted tad-pole corrected Real-part LO1 (LO2) (NLO1) NLO2 NNLO1 (NNLO2) NNNNLO1 (NNNNLO2) w.r.t beta
  cout<<"Started Imag-part of q-hat Calculation: \n";
  double PreFactor;//, Qminus=20000, 
  double MuSquare=20000.0,alphas=0.0, NC=3;
  double TOnePlusTtwo = 1.0;

  PreFactor=(8*sqrt(2)*pow(3.1415926,1.0)/(NC*TOnePlusTtwo));

  char FileName[50000];
  TGraphErrors *g1 = new TGraphErrors(File1Output.c_str(),"%lf  %*lf %*lf  %lf %lf");
  TGraphErrors *g2 = new TGraphErrors(File1Output.c_str(),"%lf  %*lf %*lf  %*lf %*lf  %*lf %*lf  %lf %lf");
  TGraphErrors *g3 = new TGraphErrors(File1Output.c_str(),"%lf  %*lf %*lf  %*lf %*lf  %*lf %*lf  %*lf %*lf  %*lf %*lf %lf %lf");
                                  //Vacuum subtracted Real-part LO1        LO2        NNLO1      NNLO2      NNNNLO1    NNNNLO2 w.r.t beta
  const int N=g1->GetN(); double *X1=g1->GetX(), *Y1=g1->GetY(), *Ex1=g1->GetEX(), *Ey1=g1->GetEY();
  int N2=g2->GetN(); double *X2=g2->GetX(), *Y2=g2->GetY(), *Ex2=g2->GetEX(), *Ey2=g2->GetEY();
  int N3=g3->GetN(); double *X3=g3->GetX(), *Y3=g3->GetY(), *Ex3=g3->GetEX(), *Ey3=g3->GetEY();
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
      alphas = StrongCouplingOneLoop(QCDType, Nt, T);
      Rbeta = RbetaFmunuSquare(QCDType, Beta);
      Factor = 1/((g0Square*pow(a0T,4.0)));
      Y[i]=PreFactor*alphas*Rbeta*Factor*( (0.5*Y1[i]) + (-Y2[i]/pow(Qminus[0]*a0,2)) + (2*Y3[i]/pow(Qminus[0]*a0,4.0)) );
      Ey[i]=PreFactor*alphas*Rbeta*Factor*sqrt(pow(0.5*Ey1[i],2.0)   + pow(Ey2[i]/pow(Qminus[0]*a0,2.0),2.0)  +  pow( 2*Ey3[i]/pow(Qminus[0]*a0,4.0),2.0)  );
      Ex[i]=0;
      X[i] = T;
      cout<<"(Beta1,Beta2,Beta3)=("<<X1[i]<<","<<X2[i]<<","<<X3[i]<<") \t (Y1,Y2,Y3)=("<<Y1[i]<<","<<Y2[i]<<","<<Y3[i]<<")"<<endl;
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
