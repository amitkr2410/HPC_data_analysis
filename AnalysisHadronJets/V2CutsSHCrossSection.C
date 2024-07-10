//Created by Amit Kumar (kumar.amit@wayne.edu)
// Requires FinalState partons E, Px, Py, Pz as input, also requires pTHatCrossSection (in mb) as input
// Output is a ROOT file which contains following plots
// 1. number of jets vs pT of the jet
// 2. differential jet cross section vs pt of the jet for each pTHat bin. Here differential cross section means= sigmapTHat*dN/(dpT*dEta*TotalEvents)
// 3. final differential jet cross section (i.e. summed over all pTHat bins) vs pt of the jet.
// 4. pTHat cross section vs pTHat bin.
// Note that the plots are saved in two format TH1F and TGraphErrors

#include <iostream>
// FastJet3 library.
#include "/wsu/home/fy/fy41/fy4125/Software/pythia8230/include/Pythia8/Pythia.h"
#include "/wsu/home/fy/fy41/fy4125/Software/pythia8230/include/Pythia8/FJcore.h"
#include <fstream>
#include <cstdio>
#include "iomanip"
#include "TMath.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

// ROOT, for saving file.
#include "TFile.h"

// Histogram coloring.
#include "TStyle.h"
#include "string"

using namespace std;         
using namespace Pythia8;

int main(int argc, char* argv[]){
   int nListJets =1;
   int StartTime = time(NULL);
// Create the ROOT application environment.
   //TApplication theApp("hist", &argc, argv);

// Create file on which histogram(s) can be saved.
// [1] InputDataFileDIR [2] OutputDataFileName  [3]=CMS or Atlas or Alice or STAR, [4]= 200, 2760 or 5020, 7000 [5]=JetRadius, [6]=EtaJetMin, [7]=EtaJetMax    [8] EtaTrackMax [9] pTTrackMin [10] EtaCutFlag =eta or y [11] ChargedFlag 0 or 1 (all or charged)  [12]=JetscapeORPythia, [13] HadronORParton [14] BGS ( 0 or 1 pp, pbpb)  [15] Centrality [16] Species= pp or PbPb pPb    [17] EtaRangeType (2 for two sided, 1 for one sided)

   char InputDataFileDIR[100000], OutputDataFileName[100000], JetscapeORPythia[100], HadronORParton[100];
   //sprintf(InputDataFileDIR,"%s",argv[1]);
   sprintf(OutputDataFileName,"%s","SHCrossSection_5020GeV_CMS_Jetscape_Cent_30-50_ChargedHadron_BGS1"); //,argv[2]);
   sprintf(JetscapeORPythia,"%s","Jetscape");//,argv[12]);
   sprintf(HadronORParton,"%s","Hadron");//,argv[13]);
   
   string ExpName="CMS";//string(argv[3]);
   int Ecm = 5020;//atoi(argv[4]);
   double JetRadius = 0;//atof(argv[1])/10.0;//atof(argv[5])/1.0;                                
   double EtaJetMin = 0.0;//atof(argv[6])/1.0;
   double EtaJetMax= 1.0;//atof(argv[7])/1.0;                                                        
   double EtaTrackMax=1.0;//atof(argv[8])/1.0;
   double pTTrackMin= 0.001;//atof(argv[9])/1.0;
   int EtaCutFlag =1;//atoi(argv[10]); //0 or 1; 0 to use rapididty cut, 1 to use eta cut
   int ChargedFlag = 1;//atoi(argv[11]); // 0 or 1; 0 means all particle, 1 means only charged
   int BGS = 1;//atoi(argv[14]);  //0 pp, 1 pbpb
   string Centrality = "30-50";//string(argv[15]);
   string PPorPbPb="PbPb";//string(argv[16]); //  Species= pp or PbPb pPb
   int EtaRangeType = 2;//atoi(argv[17]); //2 for both sided, 1 for one sided
   double JetpTMin = 5.0; //in GeV
   double AliceTrackpT=0.0;// in GeV jet with atleast single track with pT or above
   double STARTrackpT=0.0; // in GeV jet with atleast single track with pT or above
   double PI=3.1415926;
   //int pTHatMin[1]={500}, pTHatMax[1]={550}, NpTHardBin=1;

   
   int *pTHatMin;
   int *pTHatMax;
   int NpTHardBin;
   
   if(  Ecm == 200 )
     {
       //[23] == 1, 2, 3,
       int Bin1[20] = {  4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90};
       int Bin2[20] = {  5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100};
       NpTHardBin = 20;
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
	 {
	   pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
	 }
     }
   
   if(  Ecm == 2760 )
     {
       //[54]== 1, 2, 3,
       int Bin1[51] = { 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120,                               130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400,                               450, 500,   550, 600, 700, 800, 900, 1000 };
       int Bin2[51] = { 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130,                             140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450,                               500, 550, 600, 700, 800, 900, 1000, 1380 };
       NpTHardBin = 51;
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
         {
           pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
         }
     }

   if(  Ecm == 5020 )
     {  // 66= 1, 2, 3, 4, 5
       int Bin1[41] = {13, 15, 17, 20, 25, 30, 35, 40, 45,  50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,    			150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550};
       int Bin2[41] = {15, 17, 20, 25, 30, 35, 40, 45, 50,  55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,			160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600};

       //int Bin1[50] = {  13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,                       150, 160, 170, 180, 190,			 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900,  1000, 1100, 1200, 1300, 1400};
       //int Bin2[50] = {  15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,                      160, 170, 180, 190, 200,			 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500};

       // int Bin1[62] = { 5, 7, 9,  11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,                                150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600,                                  700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400};
       // int Bin2[62] = { 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,                               160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700,                                  800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510};       

       NpTHardBin = 41; //32 41 62 66
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
         {
           pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
         }
     }  
   

  double *JetpTBin;
  int NpTJetBin;

  if(  Ecm==200 && ( PPorPbPb=="PP" || PPorPbPb=="ppPbPb" || PPorPbPb=="PbPb") && ExpName=="STAR" &&  (JetRadius==0.2  ||  JetRadius==0.3 || JetRadius==0.4 ) && ChargedFlag==1)
    {
      double Bin[10]={6, 7, 8, 10, 12, 14, 16, 20, 25, 30};
      NpTJetBin = 10-1;     JetpTMin = 1.0; STARTrackpT=5.0;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }      
    }

  if(  Ecm == 2760 && (PPorPbPb=="ppPbPb" || PPorPbPb=="PbPb") && ExpName=="CMS" &&  (JetRadius==0.2  ||  JetRadius==0.3 || JetRadius==0.4) && EtaJetMax==2.0 && BGS>=0)
    {
      double Bin[13]  ={70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
      NpTJetBin = 13-1;     JetpTMin = 40;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
	{
	  JetpTBin[i] = Bin[i];    
	}      
    }



  //5.02TeV pp and PbPb

  if(  Ecm == 5020 &&  (PPorPbPb=="ppPbPb" || PPorPbPb=="PbPb") && ExpName=="CMS" && Centrality=="30-50" && BGS>=0)
    {
      double Bin[ 16] ={  5.0, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24,          		       28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6};
      NpTJetBin = 16-1;     JetpTMin = 5.0;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }


  int Error =0; //set to 1 if wants statistical error only
  //double JetpTBinWidth = (JetpTMax - JetpTMin)/NpTJetBin;
  double HardCrossSection[NpTHardBin];
  double HardCrossSectionError[NpTHardBin];
  double dNdpTCount[NpTHardBin][NpTJetBin], dNdpTCountInPlane[NpTHardBin][NpTJetBin], dNdpTCountOutPlane[NpTHardBin][NpTJetBin];   //[ptHatBin] [Regular pt] 
  double dNdpTCount2[NpTHardBin][NpTJetBin], dNdpTCount2InPlane[NpTHardBin][NpTJetBin], dNdpTCount2OutPlane[NpTHardBin][NpTJetBin];
  double dNdpTCountFull[NpTHardBin][NpTJetBin], dNdpTCountFullInPlane[NpTHardBin][NpTJetBin], dNdpTCountFullOutPlane[NpTHardBin][NpTJetBin];
  double pTHardBinJetBinError[NpTHardBin][NpTJetBin], pTHardBinJetBinInPlaneError[NpTHardBin][NpTJetBin], pTHardBinJetBinOutPlaneError[NpTHardBin][NpTJetBin];
  double pTHardBinJetBinError2[NpTHardBin][NpTJetBin], pTHardBinJetBinInPlaneError2[NpTHardBin][NpTJetBin], pTHardBinJetBinOutPlaneError2[NpTHardBin][NpTJetBin];
  double pTHardBinJetBinErrorFull[NpTHardBin][NpTJetBin], pTHardBinJetBinInPlaneErrorFull[NpTHardBin][NpTJetBin], pTHardBinJetBinOutPlaneErrorFull[NpTHardBin][NpTJetBin];
  double TotalCrossSection[NpTJetBin], TotalCrossSectionInPlane[NpTJetBin], TotalCrossSectionOutPlane[NpTJetBin];
  double TotalCrossSection2[NpTJetBin], TotalCrossSection2InPlane[NpTJetBin], TotalCrossSection2OutPlane[NpTJetBin];
  double TotalCrossSectionFull[NpTJetBin], TotalCrossSectionFullInPlane[NpTJetBin], TotalCrossSectionFullOutPlane[NpTJetBin];
  double TotalCrossSectionError[NpTJetBin], TotalCrossSectionInPlaneError[NpTJetBin], TotalCrossSectionOutPlaneError[NpTJetBin];
  double TotalCrossSectionError2[NpTJetBin], TotalCrossSectionInPlaneError2[NpTJetBin], TotalCrossSectionOutPlaneError2[NpTJetBin];
  double TotalCrossSectionErrorFull[NpTJetBin], TotalCrossSectionInPlaneErrorFull[NpTJetBin], TotalCrossSectionOutPlaneErrorFull[NpTJetBin];
  int Events[NpTHardBin], Events2[NpTHardBin], EventsFull[NpTHardBin];
  int NetJetEvents[NpTHardBin], NetJetEvents2[NpTHardBin], NetJetEventsFull[NpTHardBin];
  int HashFlag;
  double JetEta=0;  
  double dEtaJet=0;
  if(EtaRangeType == 2)   //symmetric in both direction, ex 0<|y|<1.0 or 1<|y|<2
    {
      dEtaJet=2.0*(EtaJetMax -EtaJetMin);
    }
  
  if(EtaRangeType == 1)  //asymmetric, one side, ex. -3< eta <-1, or 1 < eta <3
    {
      dEtaJet=fabs(EtaJetMax -EtaJetMin);
    }
// Histograms.
  TH1D *HistTemp = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //5GeV Bin CountVspT
  TH1D *HistTempInPlane = new TH1D("JetSpectrumBinInPlane", "Jet Spectrum pT", NpTJetBin, JetpTBin);
  TH1D *HistTempOutPlane = new TH1D("JetSpectrumBinOutPlane", "Jet Spectrum pT", NpTJetBin, JetpTBin);
  TH1D *HistTemp2 = new TH1D("JetSpectrum2Bin", "Jet Spectrum2 pT", NpTJetBin, JetpTBin);
  TH1D *HistTemp2InPlane =  new TH1D("JetSpectrum2BinInPlane", "Jet Spectrum2 pT", NpTJetBin, JetpTBin);
  TH1D *HistTemp2OutPlane = new TH1D("JetSpectrum2BinOutPlane", "Jet Spectrum2 pT", NpTJetBin, JetpTBin);
  TH1D *HistTempFull =  new TH1D("JetSpectrumFullBin", "Jet Spectrum2 pT", NpTJetBin, JetpTBin);
  TH1D *HistTempFullInPlane =  new TH1D("JetSpectrumFullBinInPlane", "Jet Spectrum2 pT", NpTJetBin, JetpTBin);
  TH1D *HistTempFullOutPlane = new TH1D("JetSpectrumFullBinOutPlane", "Jet Spectrum2 pT", NpTJetBin, JetpTBin);

  ofstream EventR, EventA, foutput, foutputInPlane, foutputOutPlane; char VarFileName[100000];
  ofstream foutput2, foutput2InPlane, foutput2OutPlane, foutputFull, foutputFullInPlane, foutputFullOutPlane;
  
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/EventRecord%s.txt",OutputDataFileName);
  EventR.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/CodeCheck%s.txt",OutputDataFileName);
  EventA.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/%s.root",OutputDataFileName); // name of the output root file               
  TFile* outFile = new TFile( VarFileName, "RECREATE");
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/High0-30V2Cuts%s.txt",OutputDataFileName);
  foutput.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/High0-30InPlaneV2Cuts%s.txt",OutputDataFileName);
  foutputInPlane.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/High0-30OutPlaneV2Cuts%s.txt",OutputDataFileName);
  foutputOutPlane.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/Low0-30V2Cuts%s.txt",OutputDataFileName);
  foutput2.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/Low0-30InPlaneV2Cuts%s.txt",OutputDataFileName);
  foutput2InPlane.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/Low0-30OutPlaneV2Cuts%s.txt",OutputDataFileName);
  foutput2OutPlane.open(VarFileName,ios::out);

  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/Full0-100V2Cuts%s.txt",OutputDataFileName);
  foutputFull.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/Full0-100InPlaneV2Cuts%s.txt",OutputDataFileName);
  foutputFullInPlane.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/Files/Full0-100OutPlaneV2Cuts%s.txt",OutputDataFileName);
  foutputFullOutPlane.open(VarFileName,ios::out);

  cout<<"Final result will be store in  File="<<VarFileName<<endl;
  EventA<<"OutputFile Name is "<<VarFileName<<endl;
  EventR<<"pTHatBin \t "<<"TotalEventsRegistered"<<endl;

  std::vector<double> RecoilParticleE, RecoilParticlePx, RecoilParticlePy, RecoilParticlePz, RecoilParticleEta, RecoilParticlePhi;
  std::vector <fjcore::PseudoJet> fjInputs;
  Pythia pythia; string CentS="30-40", InputDataFileDIRS, SoftFile;    
  double Q2VectorCuts[2]={111.0, 77.5}; double DeltaPsi=0, PhiJet=0, PSI2New=0;
  // For loop to open different pTHat bin files
  for (int k = 0; k<NpTHardBin; ++k)
    {

      char HistName[1000];
      HistTemp->Reset(); 
      sprintf(HistName,"CountVspTSpectrumBin%i_%iHigh0-30v2",pTHatMin[k],pTHatMax[k]);
      HistTemp->SetName(HistName);

      HistTempInPlane->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iHigh0-30v2InPlane",pTHatMin[k],pTHatMax[k]);
      HistTempInPlane->SetName(HistName);
      
      HistTempOutPlane->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iHigh0-30v2OutPlane",pTHatMin[k],pTHatMax[k]);
      HistTempOutPlane->SetName(HistName);

      HistTemp2->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iLow0-30v2",pTHatMin[k],pTHatMax[k]);
      HistTemp2->SetName(HistName);

      HistTemp2InPlane->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iLow0-30v2InPlane",pTHatMin[k],pTHatMax[k]);
      HistTemp2InPlane->SetName(HistName);

      HistTemp2OutPlane->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iLow0-30v2OutPlane",pTHatMin[k],pTHatMax[k]);
      HistTemp2OutPlane->SetName(HistName);

      HistTempFull->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iFull0-100v2",pTHatMin[k],pTHatMax[k]);
      HistTempFull->SetName(HistName);

      HistTempFullInPlane->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iFull0-100v2InPlane",pTHatMin[k],pTHatMax[k]);
      HistTempFullInPlane->SetName(HistName);

      HistTempFullOutPlane->Reset();
      sprintf(HistName,"CountVspTSpectrumBin%i_%iFull0-100v2OutPlane",pTHatMin[k],pTHatMax[k]);
      HistTempFullOutPlane->SetName(HistName);

      RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
      RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
      fjInputs.resize(0);
      HashFlag=0;
      Events[k] =0; Events2[k] =0; EventsFull[k] =0;
      NetJetEvents[k] = 0; NetJetEvents2[k] = 0; NetJetEventsFull[k] = 0;

      for(int CK=0; CK<2; CK++)
	{
	  char HadronFile[100000], HardCrossSectionFile[100000], pTBinString[10000];
	  CentS="30-40";
	  if(CK==1){CentS="40-50";}
	  InputDataFileDIRS = "/wsu/home/groups/maj-shen/AAPaperData/MATTER_LBT_RunningAlphaS_Q2qhat/5020_PbPb_" + string(CentS) + "_0.30_2.0_1";
	  sprintf(HadronFile,"%s/%s%sListBin%i_%i.out",InputDataFileDIRS.c_str(),JetscapeORPythia,HadronORParton, pTHatMin[k],pTHatMax[k]);
      if(PPorPbPb == "pppPb" || PPorPbPb == "pPb")
	{
	  sprintf(HardCrossSectionFile,"%s/SigmaHardBin%i_%i.out", InputDataFileDIR, pTHatMin[k], pTHatMax[k]);
	}
      else
	{
	  sprintf(HardCrossSectionFile,"/wsu/home/groups/maj-shen/AAPaperData/pTHatCrossSection_%iGeV/SigmaHardBin%i_%i.out",  Ecm, pTHatMin[k],pTHatMax[k]);
	}

      SoftFile = "/wsu/home/groups/maj-shen/AAPaperData/MATTER_LBT_RunningAlphaS_Q2qhat/SoftV2_Q2VectorAnalysis_5TeV_" + string(CentS) + ".txt";
      sprintf(pTBinString,"Current pTHatBin is %i (%i,%i) GeV",k,pTHatMin[k],pTHatMax[k]);
      EventA<<pTBinString<<endl;
      EventA<<"Hadron File ="<<HadronFile<<endl;
      EventA<<"HardCrossSectionFile="<<HardCrossSectionFile<<endl;
      EventA<<"SoftFile="<<SoftFile<<endl;
      ifstream myfile, myfile2, myfile3;
      myfile.open(HadronFile,ios::in);
      myfile2.open(HardCrossSectionFile,ios::in);
      myfile3.open(SoftFile.c_str(),ios::in);
     int  NewEvent=0, SN=0,PID;
     string x, y,z,t;
     double Px, Py, Pz, E, Eta, Phi, pStat, PSI2; double SNv2, Mult, v2soft, Psi2, q2Vector=0;
     
     // Reset
     // Parameters for the jet finders. 

// Read file
     HashFlag=0; int NetCount=0;
     string EventLabel, MyString("#");
     while ( myfile >> EventLabel  ) 
       {
	 if( MyString.compare(EventLabel)==0 && fjInputs.size()>0 )
	   {
	     if(q2Vector > Q2VectorCuts[0] ) {NetJetEvents[k] = NetJetEvents[k] +1;	     }
	     if(q2Vector < Q2VectorCuts[1] ) {NetJetEvents2[k] = NetJetEvents2[k] +1;          }
	     NetJetEventsFull[k] = NetJetEventsFull[k] +1; 
		 // Run Fastjet algorithm and sort jets in pT order.

		     if(BGS==1)		     
		       {//cout << "\n -------- List of Negative/recoil particles to be subtracted ------------------\n";
			 //cout << "i       E      pT     Eta      Phi"<<endl;
			 for (int j = 0; j < int(RecoilParticleE.size()); ++j)
			   {
			     //cout << setw(4) << j << fixed << setprecision(3) << setw(9) << RecoilParticleE[j];
			     //cout << setw(9)  << TMath::Sqrt( pow(RecoilParticlePx[j],2.0) + pow(RecoilParticlePy[j],2.0) );
			     //cout << setw(9) << RecoilParticleEta[j] << setw(9)  << RecoilParticlePhi[j] << endl;
			   }
			 //cout << "---------------------------------" << endl;
		       }		   		 		 		 
		     
		     for (int i = 0; i < fjInputs.size(); ++i)
		       {		     			 
			 PhiJet = fjInputs[i].phi()*180/PI;
			 if(PSI2 < 0) {PSI2New = 360 + (PSI2*180.0/PI);}
			 DeltaPsi = fabs(PhiJet - PSI2New);
			 if(DeltaPsi > 180) {DeltaPsi = 360 - DeltaPsi;}
			 if (DeltaPsi > 90.0) {DeltaPsi = 180 - DeltaPsi;}
			 cout<<"(Phi,PSI2New,DeltaPsi) = "<<PhiJet<<", "<<PSI2New <<", "<<DeltaPsi<<endl;
			 //High-v2 events
			 if(q2Vector > Q2VectorCuts[0] )
			   { HistTemp->Fill( fjInputs[i].perp() );				 
			     if( DeltaPsi < 30 ){HistTempInPlane->Fill(fjInputs[i].perp() );}
			     if( DeltaPsi > 60 ){HistTempOutPlane->Fill(fjInputs[i].perp() );}
			   }
			 //Low-v2 events
			 if(q2Vector < Q2VectorCuts[1] )
			   { HistTemp2->Fill( fjInputs[i].perp() );
			     if( DeltaPsi < 30 ){HistTemp2InPlane->Fill(fjInputs[i].perp() );}
			     if( DeltaPsi > 60 ){HistTemp2OutPlane->Fill(fjInputs[i].perp() );}
			   }
			 //All v2 events
			 HistTempFull->Fill(fjInputs[i].perp());
			 if( DeltaPsi < 30 ){HistTempFullInPlane->Fill(fjInputs[i].perp() );}
			 if( DeltaPsi > 60 ){HistTempFullOutPlane->Fill(fjInputs[i].perp() );}
		       }			                            			                          
		     
		     
		     
		     RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
		     RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
		     fjInputs.resize(0);
		     NewEvent = 0;
		     // cout<<"Found a Jet \t "<<pTBinString<<"\t NetJetevents is \t"<<NetJetEvents[k]<<endl;
	   }
	 
	 if(MyString.compare(EventLabel)!=0)
	   {
	     myfile >> PID >> pStat >> E >> Px >> Py >> Pz >> Eta >> Phi;
	     //cout<<" "<<PID<<" "<< pStat<< " " << E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi<<endl;
	     double PT = TMath::Sqrt( pow(Px ,2.0) + pow(Py , 2.0) );
	     //pStat=0;//Uncomment  THIS WHEN RUNNING PbPb heavy ION DATA
	     int Charge = pythia.particleData.charge( PID );
	     if(ChargedFlag==0)
	       {  
		 if( fabs(Eta) < EtaTrackMax   && pStat!=-1 && PT > pTTrackMin  && PID!=12 && PID!=14 && PID!=16 && PID!=18)
		   {
		     fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
		     fjInputs[int(fjInputs.size())-1].set_user_index(PID);
		   }
		 if( fabs(Eta) < EtaTrackMax  && pStat==-1 && PT > pTTrackMin)
		   {
		     RecoilParticleE.push_back(E);   RecoilParticlePx.push_back(Px);   RecoilParticlePy.push_back(Py);
		     RecoilParticlePz.push_back(Pz); RecoilParticleEta.push_back(Eta); RecoilParticlePhi.push_back(Phi);
		   }
	       }
	     
	     if(ChargedFlag==1)
	       {
                 if( fabs(Eta) < EtaTrackMax   && Charge!=0 && pStat!=-1 && PT > pTTrackMin  && PID!=12 && PID!=14 && PID!=16 && PID!=18){
		   fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
		   fjInputs[int(fjInputs.size())-1].set_user_index(PID);
		 }
		 if( fabs(Eta) < EtaTrackMax   && Charge!=0 && pStat==-1 && PT > pTTrackMin)
                   {
                     RecoilParticleE.push_back(E);   RecoilParticlePx.push_back(Px);   RecoilParticlePy.push_back(Py);
                     RecoilParticlePz.push_back(Pz); RecoilParticleEta.push_back(Eta); RecoilParticlePhi.push_back(Phi);
                   }
               }

	     if(HashFlag==1 )
	       {
                 HashFlag=0;
                 if(q2Vector > Q2VectorCuts[0] ){Events[k]= Events[k] +1;}
		 if(q2Vector < Q2VectorCuts[1] ){Events2[k]= Events2[k] +1;}
		 EventsFull[k]= EventsFull[k] +1;
		 NetCount =NetCount+1;
               }
	   }

	 if(MyString.compare(EventLabel)==0)
           { myfile >> PSI2; cout<<EventLabel<<" ";	    
	     for(int i=0;i<7;i++){myfile >> EventLabel; cout<<EventLabel<<" ";
	     }
	     PSI2 = int(PSI2*10000)/10000.0;
	     myfile3>> SNv2 >> Mult >> v2soft >> Psi2;
	     Psi2 = int(Psi2*10000)/10000.0;
	     if(PSI2!=Psi2){EventA<<SNv2<<"\t"<<PSI2<<"\t"<<Psi2<<" Error -- ERROR ---, NOTE Psi2 angle DOESNOT MATCH \n";}
	     q2Vector = v2soft*sqrt(Mult);
	     HashFlag=1;
	     cout<<endl;
           }
	 
       }
     myfile2 >> HardCrossSection[k]>>HardCrossSectionError[k];
     myfile.close();
     myfile2.close();
     myfile3.close();
	}
     
    
     //scalling and adding temp histograms     
     //myfile2 >> HardCrossSection[k]>>HardCrossSectionError[k];
     cout<<std::fixed<<std::setprecision(10);
     cout<<"Sigma = "<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<"\t"<<Events[k]<<"\t"<<Events2[k]<<endl; 
     EventA<<"HardCrossSection="<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl;
     if(Error == 1)
       {
	 HardCrossSectionError[k]=0.0;
       }

     cout<<"For pT-hatBin k="<<k<<endl;
     for(int j=0;j<NpTJetBin;j++)
       {	 
         dNdpTCount[k][j]= HistTemp->GetBinContent(j+1);
	 dNdpTCountInPlane[k][j]= HistTempInPlane->GetBinContent(j+1);
	 dNdpTCountOutPlane[k][j]= HistTempOutPlane->GetBinContent(j+1);
	 dNdpTCount2[k][j]= HistTemp2->GetBinContent(j+1);
	 dNdpTCount2InPlane[k][j]= HistTemp2InPlane->GetBinContent(j+1);
         dNdpTCount2OutPlane[k][j]= HistTemp2OutPlane->GetBinContent(j+1);
	 dNdpTCountFull[k][j]= HistTempFull->GetBinContent(j+1);
         dNdpTCountFullInPlane[k][j]= HistTempFullInPlane->GetBinContent(j+1);
         dNdpTCountFullOutPlane[k][j]= HistTempFullOutPlane->GetBinContent(j+1);
	 
	 if(dNdpTCount[k][j]> 5.0)
	   {
	     pTHardBinJetBinError[k][j] = (dNdpTCount[k][j]*HardCrossSection[k]/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCount[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     pTHardBinJetBinInPlaneError[k][j] = (dNdpTCountInPlane[k][j]*HardCrossSection[k]/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCountInPlane[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     pTHardBinJetBinOutPlaneError[k][j] = (dNdpTCountOutPlane[k][j]*HardCrossSection[k]/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCountOutPlane[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*JetpTBinWidth*2.0*(EtaJetMax-EtaJetMin))<<endl;
	   }
	 else
	   {
	     pTHardBinJetBinError[k][j] = 0.0;
	     pTHardBinJetBinInPlaneError[k][j] = 0.0;
	     pTHardBinJetBinOutPlaneError[k][j] = 0.0;
	     dNdpTCount[k][j]=0.0;
	     dNdpTCountInPlane[k][j] =0;
	     dNdpTCountOutPlane[k][j] = 0;
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
	   }
	 if(dNdpTCount2[k][j]> 5.0)
           {
	     pTHardBinJetBinError2[k][j] = (dNdpTCount2[k][j]*HardCrossSection[k]/(Events2[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCount2[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     pTHardBinJetBinInPlaneError2[k][j] = (dNdpTCount2InPlane[k][j]*HardCrossSection[k]/(Events2[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCount2InPlane[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     pTHardBinJetBinOutPlaneError2[k][j] = (dNdpTCount2OutPlane[k][j]*HardCrossSection[k]/(Events2[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCount2OutPlane[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	   }
	 else { 
	   pTHardBinJetBinError2[k][j] = 0.0;          dNdpTCount2[k][j]=0.0;
	   pTHardBinJetBinInPlaneError2[k][j] = 0.0;   dNdpTCount2InPlane[k][j]=0.0;
	   pTHardBinJetBinOutPlaneError2[k][j] = 0.0;  dNdpTCount2OutPlane[k][j]=0.0;
	 }

	 if(dNdpTCountFull[k][j]> 5.0)
           {
	     pTHardBinJetBinErrorFull[k][j] = (dNdpTCountFull[k][j]*HardCrossSection[k]/(EventsFull[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCountFull[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
             pTHardBinJetBinInPlaneErrorFull[k][j] = (dNdpTCountFullInPlane[k][j]*HardCrossSection[k]/(EventsFull[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCountFullInPlane[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     pTHardBinJetBinOutPlaneErrorFull[k][j] = (dNdpTCountFullOutPlane[k][j]*HardCrossSection[k]/(EventsFull[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCountFullOutPlane[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	   }
	 else {
	       pTHardBinJetBinErrorFull[k][j] = 0.0;          dNdpTCountFull[k][j]=0.0;
	       pTHardBinJetBinInPlaneErrorFull[k][j] = 0.0;   dNdpTCountFullInPlane[k][j]=0.0;
	       pTHardBinJetBinOutPlaneErrorFull[k][j] = 0.0;  dNdpTCountFullOutPlane[k][j]=0.0;
       }

	 cout<<dEtaJet <<"\t"<<HistTemp->GetBinContent(j+1)<<"\t"<<JetpTBin[j+1]-JetpTBin[j]<<endl;
	 cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet))<<endl;
	 EventR<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet))<<endl;
       }     
       HistTemp->Write(); HistTempInPlane->Write();  HistTempOutPlane->Write();  
       HistTemp2->Write(); HistTemp2InPlane->Write();  HistTemp2OutPlane->Write();
       HistTempFull->Write(); HistTempFullInPlane->Write();  HistTempFullOutPlane->Write();
     //HistTempScaled->Scale(HardCrossSection[k]/(2.0*(EtaJetMax-EtaJetMin)*Events[k]*JetpTBinWidth));
     //HistTemp->Draw();
     //HistTempScaled->Write();
     //HistFinal->Add(HistTempScaled);
     EventR<<"For bin k="<<k<<" pTHat=("<<pTHatMin[k]<<"-"<<pTHatMax[k]<<")GeV, Events= \t"<<Events[k]<<", \t"<<NetJetEvents[k]<<endl;
     EventR<<"For bin k="<<k<<" pTHat=("<<pTHatMin[k]<<"-"<<pTHatMax[k]<<")GeV, Events= \t"<<Events2[k]<<", \t"<<NetJetEvents2[k]<<endl;

// End file for k-loop
   }
  
  EventR.close();
  cout<<"Wrote the Final Histogram in ROOT file "<<endl;

  //Saving the hardpTBin crosssection vs pTBin
  double pTHatError[NpTHardBin], pTHatAvg[NpTHardBin];
  for(int k=0; k<NpTHardBin; k++)
    {
      pTHatError[k]= (pTHatMax[k] - pTHatMin[k])/2.0;
      pTHatAvg[k]= ((pTHatMin[k] + pTHatMax[k])/2.0);
    }
  TGraphErrors *GHardCrossSection = new TGraphErrors(NpTHardBin,pTHatAvg,HardCrossSection,pTHatError,HardCrossSectionError);
  GHardCrossSection->SetName("GHardCrossSection");
  GHardCrossSection->Write();

//Computing CrossSection with Error and saving as TGraphError
  double CrossSection[NpTJetBin], CrossSectionError[NpTJetBin], CrossSectionInPlane[NpTJetBin], CrossSectionInPlaneError[NpTJetBin], CrossSectionOutPlane[NpTJetBin], CrossSectionOutPlaneError[NpTJetBin], pT[NpTJetBin],pTError[NpTJetBin];
  double  CrossSection2[NpTJetBin], CrossSectionError2[NpTJetBin], CrossSection2InPlane[NpTJetBin], CrossSectionInPlaneError2[NpTJetBin], CrossSection2OutPlane[NpTJetBin], CrossSectionOutPlaneError2[NpTJetBin];
  double  CrossSectionFull[NpTJetBin], CrossSectionErrorFull[NpTJetBin], CrossSectionFullInPlane[NpTJetBin], CrossSectionInPlaneErrorFull[NpTJetBin], CrossSectionFullOutPlane[NpTJetBin], CrossSectionOutPlaneErrorFull[NpTJetBin];
  TGraphErrors *GE, *GE2;
  for(int j=0; j<NpTJetBin;j++)
    {
      CrossSection[j] =0.0; CrossSectionError[j] = 0.0; CrossSectionInPlane[j] =0.0; CrossSectionInPlaneError[j] = 0.0; CrossSectionOutPlane[j] =0.0; CrossSectionOutPlaneError[j] = 0.0;
      CrossSection2[j] =0.0; CrossSectionError2[j] = 0.0;CrossSection2InPlane[j] =0.0; CrossSectionInPlaneError2[j] = 0.0; CrossSection2OutPlane[j] =0.0; CrossSectionOutPlaneError2[j] = 0.0;
      CrossSectionFull[j] =0.0; CrossSectionErrorFull[j] = 0.0;CrossSectionFullInPlane[j] =0.0; CrossSectionInPlaneErrorFull[j] = 0.0; CrossSectionFullOutPlane[j] =0.0; CrossSectionOutPlaneErrorFull[j] = 0.0;      
      TotalCrossSection[j] =0.0; TotalCrossSectionError[j] = 0.0; TotalCrossSectionInPlane[j] =0.0; TotalCrossSectionInPlaneError[j] = 0.0; TotalCrossSectionOutPlane[j] =0.0; TotalCrossSectionOutPlaneError[j] = 0.0;
      TotalCrossSection2[j] =0.0; TotalCrossSectionError2[j] = 0.0; TotalCrossSection2InPlane[j] =0.0; TotalCrossSectionInPlaneError2[j] = 0.0; TotalCrossSection2OutPlane[j] =0.0; TotalCrossSectionOutPlaneError2[j] = 0.0;
      TotalCrossSectionFull[j] =0.0; TotalCrossSectionFullInPlane[j] =0.0; TotalCrossSectionFullOutPlane[j] =0.0;
      TotalCrossSectionErrorFull[j] = 0.0; TotalCrossSectionInPlaneErrorFull[j] = 0.0; TotalCrossSectionOutPlaneErrorFull[j] = 0.0;
      pT[j] = 0.0;
      pTError[j] = 0.0;
    }
  cout<<"Initialized array's to zero"<<endl;
  for(int k=0;k<NpTHardBin;k++)
    {cout<<"For ptHardBin = "<<k+1<<"\t CrossSection is Below "<<endl;
      for(int j=0; j<NpTJetBin;j++)
	{
	  CrossSection[j] = (dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
	  CrossSectionError[j] = pTHardBinJetBinError[k][j];
	  TotalCrossSection[j] = TotalCrossSection[j] + CrossSection[j];
	  TotalCrossSectionError[j] = TotalCrossSectionError[j] + TMath::Power( pTHardBinJetBinError[k][j], 2.0);

	  CrossSectionInPlane[j] = (dNdpTCountInPlane[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
          CrossSectionInPlaneError[j] = pTHardBinJetBinInPlaneError[k][j];
          TotalCrossSectionInPlane[j] = TotalCrossSectionInPlane[j] + CrossSectionInPlane[j];
	  TotalCrossSectionInPlaneError[j] = TotalCrossSectionInPlaneError[j] + TMath::Power( pTHardBinJetBinInPlaneError[k][j], 2.0);

	  CrossSectionOutPlane[j] = (dNdpTCountOutPlane[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
	  CrossSectionOutPlaneError[j] = pTHardBinJetBinOutPlaneError[k][j];
	  TotalCrossSectionOutPlane[j] = TotalCrossSectionOutPlane[j] + CrossSectionOutPlane[j];
	  TotalCrossSectionOutPlaneError[j] = TotalCrossSectionOutPlaneError[j] + TMath::Power( pTHardBinJetBinOutPlaneError[k][j], 2.0);

	  CrossSection2[j] = (dNdpTCount2[k][j]*HardCrossSection[k])/(Events2[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
          CrossSectionError2[j] = pTHardBinJetBinError2[k][j];
          TotalCrossSection2[j] = TotalCrossSection2[j] + CrossSection2[j];
          TotalCrossSectionError2[j] = TotalCrossSectionError2[j] + TMath::Power( pTHardBinJetBinError2[k][j], 2.0);

	  CrossSection2InPlane[j] = (dNdpTCount2InPlane[k][j]*HardCrossSection[k])/(Events2[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
          CrossSectionInPlaneError2[j] = pTHardBinJetBinInPlaneError2[k][j];
          TotalCrossSection2InPlane[j] = TotalCrossSection2InPlane[j] + CrossSection2InPlane[j];
	  TotalCrossSectionInPlaneError2[j] = TotalCrossSectionInPlaneError2[j] + TMath::Power( pTHardBinJetBinInPlaneError2[k][j], 2.0);

	  CrossSection2OutPlane[j] = (dNdpTCount2OutPlane[k][j]*HardCrossSection[k])/(Events2[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
	  CrossSectionOutPlaneError2[j] = pTHardBinJetBinOutPlaneError2[k][j];
          TotalCrossSection2OutPlane[j] = TotalCrossSection2OutPlane[j] + CrossSection2OutPlane[j];
          TotalCrossSectionOutPlaneError2[j] = TotalCrossSectionOutPlaneError2[j] + TMath::Power( pTHardBinJetBinOutPlaneError2[k][j], 2.0);

	  CrossSectionFull[j] = (dNdpTCountFull[k][j]*HardCrossSection[k])/(EventsFull[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
          CrossSectionErrorFull[j] = pTHardBinJetBinErrorFull[k][j];
          TotalCrossSectionFull[j] = TotalCrossSectionFull[j] + CrossSectionFull[j];
          TotalCrossSectionErrorFull[j] = TotalCrossSectionErrorFull[j] + TMath::Power( pTHardBinJetBinErrorFull[k][j], 2.0);
	  
	  CrossSectionFullInPlane[j] = (dNdpTCountFullInPlane[k][j]*HardCrossSection[k])/(EventsFull[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
          CrossSectionInPlaneErrorFull[j] = pTHardBinJetBinInPlaneErrorFull[k][j];
          TotalCrossSectionFullInPlane[j] = TotalCrossSectionFullInPlane[j] + CrossSectionFullInPlane[j];
	  TotalCrossSectionInPlaneErrorFull[j] = TotalCrossSectionInPlaneErrorFull[j] + TMath::Power( pTHardBinJetBinInPlaneErrorFull[k][j], 2.0);

	  CrossSectionFullOutPlane[j] = (dNdpTCountFullOutPlane[k][j]*HardCrossSection[k])/(EventsFull[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
          CrossSectionOutPlaneErrorFull[j] = pTHardBinJetBinOutPlaneErrorFull[k][j];
          TotalCrossSectionFullOutPlane[j] = TotalCrossSectionFullOutPlane[j] + CrossSectionFullOutPlane[j];
          TotalCrossSectionOutPlaneErrorFull[j] = TotalCrossSectionOutPlaneErrorFull[j] + TMath::Power( pTHardBinJetBinOutPlaneErrorFull[k][j], 2.0);

	  
	  pT[j] = JetpTBin[j] + ((JetpTBin[j+1]-JetpTBin[j])/2.0);
	  pTError[j] = (JetpTBin[j+1]-JetpTBin[j])/2.0;
	  // cout<<pT[j]<<"\t"<<CrossSection[j]<<"\t"<<CrossSectionError[j]<<endl;
	}
      
      GE  = new TGraphErrors(NpTJetBin,pT,CrossSection,pTError,CrossSectionError);
      GE2 = new TGraphErrors(NpTJetBin,pT,CrossSection2,pTError,CrossSectionError2);
      char MyGraphName[100];
      sprintf(MyGraphName,"JetCrossSectionGraphErrorBin%i_%i",pTHatMin[k],pTHatMax[k]);
      GE->SetName(MyGraphName);
      sprintf(MyGraphName,"JetCrossSectionGraphError2Bin%i_%i",pTHatMin[k],pTHatMax[k]);
      GE2->SetName(MyGraphName);
      GE->Write();      GE2->Write();
    }

  //Final Plot TotalJetCrossSection with error
  cout<<"Final results for total crossSection"<<endl;
  cout<<"JetpT \t"<<"TotalCrossSection \t"<<"Error"<<endl;
  foutput<<"#JetpT \t"<<"TotalCrossSection \t"<<"JetpTBinWidth/2.0 \t"<<"TCSError \t"<<"TCSErr%"<<endl;
  for(int j=0;j<NpTJetBin;j++)
    {
      TotalCrossSectionError[j] =TMath::Sqrt(TotalCrossSectionError[j]);
      TotalCrossSectionInPlaneError[j] =TMath::Sqrt(TotalCrossSectionInPlaneError[j]);
      TotalCrossSectionOutPlaneError[j] =TMath::Sqrt(TotalCrossSectionOutPlaneError[j]);
      TotalCrossSectionError2[j] =TMath::Sqrt(TotalCrossSectionError2[j]);
      TotalCrossSectionInPlaneError2[j] =TMath::Sqrt(TotalCrossSectionInPlaneError2[j]);
      TotalCrossSectionOutPlaneError2[j] =TMath::Sqrt(TotalCrossSectionOutPlaneError2[j]);
      TotalCrossSectionErrorFull[j] =TMath::Sqrt(TotalCrossSectionErrorFull[j]);
      TotalCrossSectionInPlaneErrorFull[j] =TMath::Sqrt(TotalCrossSectionInPlaneErrorFull[j]);
      TotalCrossSectionOutPlaneErrorFull[j] =TMath::Sqrt(TotalCrossSectionOutPlaneErrorFull[j]);

      cout<<pT[j]<<"\t"<<TotalCrossSection[j]<<" +/- \t"<<TotalCrossSectionError[j]<<endl;
      cout<<pT[j]<<"\t"<<TotalCrossSection2[j]<<" +/- \t"<<TotalCrossSectionError2[j]<<endl;
      foutput<<pT[j]<<"\t"<<TotalCrossSection[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionError[j]<<"\t"<<TotalCrossSectionError[j]*100/TotalCrossSection[j]<<endl;
      foutputInPlane<<pT[j]<<"\t"<<TotalCrossSectionInPlane[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionInPlaneError[j]<<"\t"<<TotalCrossSectionInPlaneError[j]*100/TotalCrossSectionInPlane[j]<<endl;
      foutputOutPlane<<pT[j]<<"\t"<<TotalCrossSectionOutPlane[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionOutPlaneError[j]<<"\t"<<TotalCrossSectionOutPlaneError[j]*100/TotalCrossSectionOutPlane[j]<<endl;
      foutput2<<pT[j]<<"\t"<<TotalCrossSection2[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionError2[j]<<"\t"<<TotalCrossSectionError2[j]*100/TotalCrossSection2[j]<<endl;
      foutput2InPlane<<pT[j]<<"\t"<<TotalCrossSection2InPlane[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionInPlaneError2[j]<<"\t"<<TotalCrossSectionInPlaneError2[j]*100/TotalCrossSection2InPlane[j]<<endl;
      foutput2OutPlane<<pT[j]<<"\t"<<TotalCrossSection2OutPlane[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionOutPlaneError2[j]<<"\t"<<TotalCrossSectionOutPlaneError2[j]*100/TotalCrossSection2OutPlane[j]<<endl;
      
      foutputFull<<pT[j]<<"\t"<<TotalCrossSectionFull[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionErrorFull[j]<<"\t"<<TotalCrossSectionErrorFull[j]*100/TotalCrossSectionFull[j]<<endl;
      foutputFullInPlane<<pT[j]<<"\t"<<TotalCrossSectionFullInPlane[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionInPlaneErrorFull[j]<<"\t"<<TotalCrossSectionInPlaneErrorFull[j]*100/TotalCrossSectionFullInPlane[j]<<endl;
      foutputFullOutPlane<<pT[j]<<"\t"<<TotalCrossSectionFullOutPlane[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionOutPlaneErrorFull[j]<<"\t"<<TotalCrossSectionOutPlaneErrorFull[j]*100/TotalCrossSectionFullOutPlane[j]<<endl;
    }
  cout<<"Create Final TGraphError plot "<<endl;
  //Final Plot and save also
  TGraphErrors *GETotal = new TGraphErrors(NpTJetBin,pT,TotalCrossSection,pTError,TotalCrossSectionError);
  GETotal->SetName("TotalJetCrossSectionGraphError");
  TGraphErrors *GETotalInPlane = new TGraphErrors(NpTJetBin,pT,TotalCrossSectionInPlane,pTError,TotalCrossSectionInPlaneError);
  GETotalInPlane->SetName("TotalJetCrossSectionGraphErrorInPlane");
  TGraphErrors *GETotalOutPlane = new TGraphErrors(NpTJetBin,pT,TotalCrossSectionOutPlane,pTError,TotalCrossSectionOutPlaneError);
  GETotalOutPlane->SetName("TotalJetCrossSectionGraphErrorOutPlane");

  TGraphErrors *GETotal2 = new TGraphErrors(NpTJetBin,pT,TotalCrossSection2,pTError,TotalCrossSectionError2);
  GETotal2->SetName("TotalJetCrossSectionGraphError2");
  TGraphErrors *GETotal2InPlane = new TGraphErrors(NpTJetBin,pT,TotalCrossSection2InPlane,pTError,TotalCrossSectionInPlaneError2);
  GETotal2InPlane->SetName("TotalJetCrossSectionGraphError2InPlane");
  TGraphErrors *GETotal2OutPlane = new TGraphErrors(NpTJetBin,pT,TotalCrossSection2OutPlane,pTError,TotalCrossSectionOutPlaneError2);
  GETotal2OutPlane->SetName("TotalJetCrossSectionGraphError2OutPlane");
  
  TGraphErrors *GETotalFull = new TGraphErrors(NpTJetBin,pT,TotalCrossSectionFull,pTError,TotalCrossSectionErrorFull);
  GETotalFull->SetName("TotalJetCrossSectionGraphErrorFull");
  TGraphErrors *GETotalFullInPlane = new TGraphErrors(NpTJetBin,pT,TotalCrossSectionFullInPlane,pTError,TotalCrossSectionInPlaneErrorFull);
  GETotalFullInPlane->SetName("TotalJetCrossSectionGraphErrorFullInPlane");
  TGraphErrors *GETotalFullOutPlane = new TGraphErrors(NpTJetBin,pT,TotalCrossSectionFullOutPlane,pTError,TotalCrossSectionOutPlaneErrorFull);
  GETotalFullOutPlane->SetName("TotalJetCrossSectionGraphErrorFullOutPlane");
  GETotal->Write();  GETotalInPlane->Write(); GETotalOutPlane->Write();
  GETotal2->Write(); GETotal2InPlane->Write(); GETotal2OutPlane->Write();
  GETotalFull->Write(); GETotalFullInPlane->Write(); GETotalFullOutPlane->Write();
  foutput.close(); foutputInPlane.close();foutputOutPlane.close();
  foutput2.close(); foutput2InPlane.close();foutput2OutPlane.close();
  foutputFull.close(); foutputFullInPlane.close();foutputFullOutPlane.close();


  // Deleting Dynamic pointers
  delete HistTemp, HistTempInPlane, HistTempOutPlane, HistTemp2, HistTemp2InPlane, HistTemp2OutPlane; 
  delete HistTempFull, HistTempFullInPlane, HistTempFullOutPlane;
  delete GE, GE2; delete GETotalFull, GETotalFullInPlane, GETotalFullOutPlane;
  delete GETotal, GETotalInPlane, GETotalOutPlane, GETotal2, GETotal2InPlane, GETotal2OutPlane;
  delete GHardCrossSection;
  outFile->Close();
  //Done
  int EndTime = time(NULL);
  int Hour = (EndTime-StartTime)/3600;
  int Minute = ((EndTime-StartTime)/60)-Hour*60;
  int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
  cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
  EventR<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
  return 0;
}
