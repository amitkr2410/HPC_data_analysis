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
//#include "Pythia8Plugins/FastJet3.h"
//#include "fastjet/PseudoJet.hh"
//#include "fastjet/ClusterSequence.hh"

#include <fstream>
#include <cstdio>
#include "iomanip"
#include "TMath.h"

// ROOT, for histogramming.
#include"TString.h"
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
using namespace fjcore;
int main(int argc, char **argv){
   int nListJets =1;
   int StartTime = time(NULL);
   //TApplication theApp("hist", &argc, argv);

// Create file on which histogram(s) can be saved.
   char outFileName[50000];

   int Ecm = atoi(argv[4]);
   int pTHatMin[1]={25};
   int pTHatMax[1]={30};
   int NpTHardBin=1; 
   
   //int *pTHatMin;
   //int *pTHatMax;
   //int NpTHardBin;
   /*
   if(  Ecm == 200 )
     {
       int Bin1[23] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90};
       int Bin2[23] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, -1};
       NpTHardBin = 23;
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
	 {
	   pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
	 }
     }
   
   if(  Ecm == 2760 )
     {
       int Bin1[54] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120,                               130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400,                               450, 500,   550, 600, 700, 800, 900, 1000 };
       int Bin2[54] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130,                             140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450,                               500, 550, 600, 700, 800, 900, 1000, -1};
       NpTHardBin = 54;
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
         {
           pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
         }
     }

   if(  Ecm == 5020 )
     {
       int Bin1[66] = { 1, 2, 3, 4, 5, 7, 9,  11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,                                150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600,                                  700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400};
       int Bin2[66] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,                               160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700,                                  800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510};       
       NpTHardBin = 66;
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
         {
           pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
         }
     }

   if(  Ecm == 7000 )
     {
       int Bin1[69] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180,                            190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1100, 1200, 1300,                             1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000 };
       int Bin2[69] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,                        150, 160, 170, 180, 190, 200,  210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550,                         600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600,                           2800, 3000, -1};
       NpTHardBin = 69;       
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
         {
           pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
         }       
     }
   */
     //[1]=JetRadius,  [2]=Jetscape or Pythia, [3]=CMS or Atlas or Alice or Star, [4]= 200 or 2760 or 7000, [5]=EtaJetMin, [6]=EtaJetMax, [7]EtaCutFlag, [8]=ChargedFlag=0 or 1 (all or charged), [9]Hadron or Parton, [10]FileNameOutPut  [11] EtaTrackMax [12] pTTrackMin   [13]HadronDirectoryRun

  double JetRadius = atof(argv[1])/1.0;
  int EtaCutFlag =atoi(argv[7]); //0 or 1; 0 to use rapididty cut, 1 to use eta cut
  double EtaJetMin = atof(argv[5])/1.0;
  double EtaJetMax=atof(argv[6])/1.0;
  int ChargedFlag = atoi(argv[8]); // 0 or 1; 0 means all particle, 1 means only charged
  double EtaTrackMax=atof(argv[11])/1.0;
  double pTTrackMin= atof(argv[12])/1.0;
  // int Ecm = atoi(argv[4]); //Center-of-mass energy of PP collision
  
  string ExpName=std::string(argv[3]);
  double JetpTMin = 10; //in GeV //double JetpTMax = 500; //in GeV

  double *JetpTBin;
  int NpTJetBin;
  
  sprintf(outFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/%s.root",argv[10]);
  //sprintf(outFileName,"RootFiles/JetCrossSection_R_0p4_AntikT_pp_JETSCAPE_CMS_2760GeV_EtaJetMin0EtaMax2_FS_allHadron.root",argv  [1],argv[3],argv[2]);  
  TFile* outFile = new TFile(outFileName, "RECREATE");

  if(  Ecm == 200 && ExpName=="Star" &&  JetRadius==0.6  && (EtaJetMax==0.5 || EtaJetMax==1.0) )
    {      
      double Bin[11]  = {10, 12, 14, 16, 19, 22, 26.5, 32, 37, 44, 52};
      NpTJetBin = 11-1;     JetpTMin = 1;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 2760 && ExpName=="CMS" &&  (JetRadius==0.2  ||  JetRadius==0.3 || JetRadius==0.4) && EtaJetMax==2.0)
    {
      double Bin[13]  ={70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
      NpTJetBin = 13-1;     JetpTMin = 10;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
	{
	  JetpTBin[i] = Bin[i];    
	}      
    }

  if(  Ecm == 2760 && ExpName=="CMS" &&  JetRadius==0.7 && EtaJetMax==0.5)
    {       
      double Bin[20] = {74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592} ;
      NpTJetBin = 20-1;     JetpTMin = 10; 
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 2760 && ExpName=="Alice" &&  (JetRadius==0.2  || JetRadius==0.4) && EtaJetMax==0.5)
    {
      double Bin[10] = {20, 25, 30, 35, 40, 50, 60, 70, 85, 105};
      NpTJetBin = 10-1;     JetpTMin = 5;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 5020 && ExpName=="Atlas" &&  JetRadius==0.4 && (EtaJetMax==0.3 || EtaJetMax==2.8) )
    {
      double Bin[15] = {40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501,631, 800, 1000};
      NpTJetBin = 15-1;     JetpTMin = 5;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 5020 && ExpName=="Atlas" &&  JetRadius==0.4 && (EtaJetMin==1.6 && EtaJetMax==2.1) )
    {
      double Bin[13] = {40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501, 631};
      NpTJetBin = 13-1;     JetpTMin = 5;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 7000 && ExpName=="CMS" &&  JetRadius==0.7 && EtaJetMax==2.5)
    {
      double Bin[22] = {56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592};
      NpTJetBin = 22-1;       JetpTMin = 10;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 7000 && ExpName=="CMS" &&  (JetRadius==0.7 || JetRadius==0.5) && EtaJetMax==0.5)
    {
      double Bin[34] = {56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905,  967, 1032, 1101, 1172, 1248, 1327};
      NpTJetBin = 34-1;       JetpTMin = 10;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
	{
          JetpTBin[i] = Bin[i];
        }
    }
  

  if(  Ecm == 7000 && ExpName=="Atlas" &&  JetRadius==0.6 && EtaJetMax==0.3 && EtaCutFlag==0)
    {      
      double Bin[10] = {60, 80, 110, 160, 210, 260, 310, 400, 500, 600};
      NpTJetBin = 10-1;       JetpTMin = 10;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 7000 && ExpName=="Atlas" &&  ( (JetRadius==0.4 && EtaJetMax==0.3) || (JetRadius==0.4 && EtaJetMax==2.1) || (JetRadius==0.6 && EtaJetMax==2.1)  ) && EtaCutFlag==0)
    {
      double Bin[9] = {60, 80, 110, 160, 210, 260, 310, 400, 500};
      NpTJetBin = 9-1;       JetpTMin = 10;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 7000 && ExpName=="Alice" &&  ((JetRadius==0.2 && EtaJetMax==0.7 && EtaTrackMax==0.9 && pTTrackMin==0.15) || (JetRadius==0.4 && EtaJetMax==0.5 && EtaTrackMax==0.9 && pTTrackMin==0.15 )) && EtaCutFlag==1 && ChargedFlag==1)
    {
      double Bin[12] = {20, 24, 28, 32, 38, 44, 50, 58, 68, 74, 88, 98};
      NpTJetBin = 12-1;      JetpTMin = 5;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 7000 && ExpName=="Alice" &&  ((JetRadius==0.2 && EtaJetMax==0.3 && EtaTrackMax==0.9 && pTTrackMin==0.15) || (JetRadius==0.4 && EtaJetMax==0.3 && EtaTrackMax==0.9 && pTTrackMin==0.15 )) && EtaCutFlag==1 && ChargedFlag==1)
    {
      double Bin[8] = {20, 24, 28, 32, 38,  50, 66, 100};
      NpTJetBin = 8-1;      JetpTMin = 5;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }


  if(  Ecm == 7000 && ExpName=="Atlas" &&  JetRadius==0.4 && EtaJetMax==0.5  && EtaCutFlag==0 && ChargedFlag==1)
    {
      double Bin[34] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 45, 50, 60, 70, 80, 90, 100};
      NpTJetBin = 34-1;          JetpTMin = 1;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }


  int BGS=0; //for recoil substraction set it to 1
  int Error =0; //set to 1 if wants statistical error only
  //double JetpTBinWidth = (JetpTMax - JetpTMin)/NpTJetBin;
  double HardCrossSection[NpTHardBin];
  double HardCrossSectionError[NpTHardBin];
  double dNdpTCount[NpTHardBin][NpTJetBin];  //[ptHatBin] [Regular pt]
  double pTHardBinJetBinError[NpTHardBin][NpTJetBin];
  double TotalCrossSection[NpTJetBin];
  double TotalCrossSectionError[NpTJetBin];
  int Events[NpTHardBin];
  int NetJetEvents[NpTHardBin];
// Histograms.
  TH1D *HistTemp = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //5GeV Bin CountVspT
  TH1D *HistTempScaled = new TH1D("JetSpectrumScaledBin", "Jet Spectrum Scaled pT", NpTJetBin, JetpTBin); // Scaled VspT
  TH1D *HistTempPhoton; 
  TH1D *HistTempPhotonScaled;
  TGraphErrors *GE;

  ofstream EventR, foutput;

  char EventFile[50000];sprintf(EventFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/EventRecord%s.txt",argv[10]);
  EventR.open(EventFile,ios::out);

  sprintf(EventFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/%s.txt",argv[10]);
  foutput.open(EventFile,ios::out);
  EventR<<"pTHatBin \t "<<"TotalEventsRegistered"<<endl;

  std::vector<double> RecoilParticleE, RecoilParticlePx, RecoilParticlePy, RecoilParticlePz, RecoilParticleEta, RecoilParticlePhi;
  std::vector <PseudoJet> fjInputs;
  Pythia pythia;    

  // For loop to open different pTHat bin files
  for (int k = 0; k<NpTHardBin; ++k)
    {
      char PartonFile[30000], HardCrossSectionFile[30000], pTBinString[10000];
      sprintf(PartonFile,"%s/%s%sListBin%i_%i.out",argv[13], argv[2],argv[9], pTHatMin[k],pTHatMax[k]);
      sprintf(HardCrossSectionFile,"%s/SigmaHardBin%i_%i.out",argv[13],  pTHatMin[k],pTHatMax[k]);
      //sprintf(PartonFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/%s",argv[13]);
      //sprintf(HardCrossSectionFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/FakePartonColored200GeV/SigmaHardFakePartonPz100.out");

      sprintf(pTBinString,"Current pTHatBin is %i (%i,%i) GeV",k,pTHatMin[k],pTHatMax[k]);      
      cout<<pTBinString<<"\n"<<HardCrossSectionFile<<endl;
      cout<<PartonFile<<"\n"<<argv[13]<<endl;
      ifstream myfile, myfile2;
      myfile.open(PartonFile,ios::in);
      myfile2.open(HardCrossSectionFile,ios::in);

     int  NewEvent=0, SN=0,PID;
     string x, y,z,t;
     double Px, Py, Pz, E, Eta, Phi, pStat;
     
     // Reset
     char HistName[100];
     HistTemp->Reset();
     sprintf(HistName,"CountVspTSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
     HistTemp->SetName(HistName);
     HistTempScaled->Reset();
     sprintf(HistName,"WeightedJetSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
     HistTempScaled->SetName(HistName);
     RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
     RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
     fjInputs.resize(0);

     // Parameters for the jet finders. 
     JetDefinition jetDef(antikt_algorithm, JetRadius);

// Read file
    
     Events[k] =0;
     NetJetEvents[k] = 0;
     string EventLabel, MyString("#");
     myfile >>EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>>EventLabel>> EventLabel>> EventLabel>>EventLabel;

     while ( myfile >> EventLabel ) 
       {
	 //cout<<EventLabel;
          	 
	 if( MyString.compare(EventLabel)==0 && fjInputs.size()>0 )
	   {
	     Events[k] = Events[k] +1;	     
	     //cout<<pTBinString<<"\t event is \t"<<Events[k]<<endl;
	     //cout<<"(E,Px,Py,Pz,Eta) = "<<E<<"\t"<<Px<<"\t"<<Py<<"\t"<<Pz<<"\t"<<Eta<<endl;
	     //cout<<EventLabel << " "<<PID<<" "<< pStat<< " " << E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi<<endl;

		 // Run Fastjet algorithm and sort jets in pT order.
	     vector <fjcore::PseudoJet> UnSortedJets, SortedJets;
		 fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
		 UnSortedJets = clustSeq.inclusive_jets(JetpTMin);
		 SortedJets    = sorted_by_pt(UnSortedJets);
		 
		 // List first few FastJet jets and some info about them.
		 if (nListJets) 
		   {
		     cout << "\n ---FastJet jets before recoil substraction, R = " << JetRadius << "\t anti-kt with pTHatBin = "<<k<<"\t JetpTMin="<<JetpTMin
			  << "  --------------------------------------------------\n\n "
			  << "  i         pT        y     eta      phi  " << endl;
		     for (int i = 0; i < int(SortedJets.size()); ++i) 
		       {
			 vector<fjcore::PseudoJet> constituents = SortedJets[i].constituents();			 
			 cout << setw(4) << i << fixed << setprecision(3) << setw(11)
			      << SortedJets[i].perp() << setw(9)  << SortedJets[i].rap()
			      << setw(9) << SortedJets[i].eta() << setw(9)  << SortedJets[i].phi_std() << endl;
		       }		       
		     cout << "\n --------  End FastJet Listing  ------------------"
			  << "---------------------------------" << endl;
		     if(BGS==1)		     
		       {cout << "\n -------- List of Negative/recoil particles to be subtracted ------------------\n"
			     << "i        E        pT         Eta            Phi"<<endl;
			 for (int j = 0; j < int(RecoilParticleE.size()); ++j)
			   {
			     cout << setw(4) << j << fixed << setprecision(3) << setw(11)
				  << RecoilParticleE[j] << setw(9)  << TMath::Sqrt( pow(RecoilParticlePx[j],2.0) + pow(RecoilParticlePy[j],2.0) )
				  << setw(9) << RecoilParticleEta[j] << setw(9)  << RecoilParticlePhi[j] << endl;
			   }
			 cout << "---------------------------------" << endl;
		       }
		   }
		 		 
		 int pFast = SortedJets.size();
		 if(BGS==1)
                   {
                     for (int j = 0; j < RecoilParticleEta.size(); ++j)
                       {
                         for (int i = 0; i < pFast; ++i)
                           {
                             double RecoilRadius = sqrt(pow(SortedJets[i].eta() - RecoilParticleEta[j]  ,2.0 ) + pow(SortedJets[i].phi() - RecoilParticlePhi[j],2.0));
                             if(RecoilRadius < JetRadius)
                               {
                                 double p1 = SortedJets[i].px() - RecoilParticlePx[j];
                                 double p2 = SortedJets[i].py() - RecoilParticlePy[j];
                                 double p3 = SortedJets[i].pz() - RecoilParticlePz[j];
                                 double p0 = SortedJets[i].e() - RecoilParticleE[j];
                                 SortedJets[i].reset_momentum(p1,p2,p3,p0);
                                 i=pFast;
                               }
                           }
                       }
		     if (nListJets)
                       {
                         cout << "\n --  FastJet jets after Subtracting recoil, R = " << JetRadius << "\t anti-kt for pTHatBin = "<<k
                              << "  --------------------------------------------------\n\n "
                              << "  i         pT        y     eta      phi  " << endl;
                         for (int i = 0; i < int(SortedJets.size()); ++i)
                           {
                             cout << setw(4) << i << fixed << setprecision(3) << setw(11)
                                  << SortedJets[i].perp() << setw(9)  << SortedJets[i].rap()
                                  << setw(9) << SortedJets[i].eta() << setw(9)  << SortedJets[i].phi_std() << endl;
                           }
                         cout << "\n --------  End FastJet Listing  ------------------"
                              << "---------------------------------" << endl;
                       }
		   }

		 for (int i = 0; i < pFast; ++i)
		   {
		     if( EtaCutFlag==1 && EtaJetMin <= fabs(SortedJets[i].eta()) && fabs(SortedJets[i].eta()) < EtaJetMax )
		       {
			 HistTemp->Fill( SortedJets[i].perp() );       		    
			 HistTempScaled->Fill( SortedJets[i].perp() );
		       }
		     if( EtaCutFlag==0 && EtaJetMin <= fabs(SortedJets[i].rap()) && fabs(SortedJets[i].rap()) < EtaJetMax )
                       {
                         HistTemp->Fill( SortedJets[i].perp() );
                         HistTempScaled->Fill( SortedJets[i].perp() );
                       }
		   }		 
		 NetJetEvents[k]= NetJetEvents[k] +1;
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
	     pStat=0;//Uncomment  THIS WHEN RUNNING PbPb heavy ION DATA
	     int Charge = pythia.particleData.charge( PID );
	     if(ChargedFlag==0)
	       {  
		 if( fabs(Eta) < EtaTrackMax   && pStat>=0 && PT > pTTrackMin  && PID!=12 && PID!=14 && PID!=16 && PID!=18){fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));}
	       }
	     
	     if(ChargedFlag==1)
	       {
                 if( fabs(Eta) < EtaTrackMax   && Charge!=0 && pStat>=0 && PT > pTTrackMin  && PID!=12 && PID!=14 && PID!=16 && PID!=18){fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));}
               }

	     if( fabs(Eta) < EtaTrackMax  && pStat==-1)
	       {
		 RecoilParticleE.push_back(E);   RecoilParticlePx.push_back(Px);   RecoilParticlePy.push_back(Py);
		 RecoilParticlePz.push_back(Pz); RecoilParticleEta.push_back(Eta); RecoilParticlePhi.push_back(Phi);
	       }
	   }

	 if(MyString.compare(EventLabel)==0)
           {	    
	     for(int i=0;i<8;i++){myfile >> EventLabel; //cout<<EventLabel<<" ";
	     }
	     //cout<<endl;
           }
	 
       }
     
    
     //scalling and adding temp histograms     
     myfile2 >> HardCrossSection[k]>>HardCrossSectionError[k];
     cout<<std::fixed<<std::setprecision(10);
     cout<<"Sigma = "<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl;     
     if(Error == 1)
       {
	 HardCrossSectionError[k]=0.0;
       }

     for(int j=0;j<NpTJetBin;j++)
       {	 
         dNdpTCount[k][j]= HistTemp->GetBinContent(j+1);
	 if(dNdpTCount[k][j]> 0.0)
	   {
	     pTHardBinJetBinError[k][j] = (dNdpTCount[k][j]*HardCrossSection[k]/(NetJetEvents[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*(EtaJetMax-EtaJetMin)))*TMath::Sqrt( (1/dNdpTCount[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(NetJetEvents[k]*JetpTBinWidth*2.0*(EtaJetMax-EtaJetMin))<<endl;
	   }
	 else
	   {
	     pTHardBinJetBinError[k][j] = 0.0;
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
	   }
       }     
     HistTemp->Write();     
     //HistTempScaled->Scale(HardCrossSection[k]/(2.0*(EtaJetMax-EtaJetMin)*NetJetEvents[k]*JetpTBinWidth));
     //HistTemp->Draw();
     //HistTempScaled->Write();
     //HistFinal->Add(HistTempScaled);
     EventR<<pTHatMin[k]<<"-"<<pTHatMax[k]<<"\t"<<NetJetEvents[k]<<endl;
     myfile.close();
     myfile2.close();

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
  double CrossSection[NpTJetBin],CrossSectionError[NpTJetBin],pT[NpTJetBin],pTError[NpTJetBin];  
  for(int j=0; j<NpTJetBin;j++)
    {
      CrossSection[j] =0.0;
      CrossSectionError[j] = 0.0;
      TotalCrossSection[j] =0.0;
      TotalCrossSectionError[j] = 0.0;
      pT[j] = 0.0;
      pTError[j] = 0.0;
    }
  cout<<"Initialized array's to zero"<<endl;
  for(int k=0;k<NpTHardBin;k++)
    {cout<<"For ptHardBin = "<<k+1<<"\t CrossSection is Below "<<endl;
      for(int j=0; j<NpTJetBin;j++)
	{
	  CrossSection[j] = (dNdpTCount[k][j]*HardCrossSection[k])/(NetJetEvents[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*(EtaJetMax-EtaJetMin));
	  CrossSectionError[j] = pTHardBinJetBinError[k][j];
	  TotalCrossSection[j] = TotalCrossSection[j] + CrossSection[j];
	  TotalCrossSectionError[j] = TotalCrossSectionError[j] + TMath::Power( pTHardBinJetBinError[k][j], 2.0);
	  pT[j] = JetpTBin[j] + ((JetpTBin[j+1]-JetpTBin[j])/2.0);
	  pTError[j] = (JetpTBin[j+1]-JetpTBin[j])/2.0;
	  // cout<<pT[j]<<"\t"<<CrossSection[j]<<"\t"<<CrossSectionError[j]<<endl;
	}
      
      GE = new TGraphErrors(NpTJetBin,pT,CrossSection,pTError,CrossSectionError);
      char MyGraphName[100];
      sprintf(MyGraphName,"JetCrossSectionGraphErrorBin%i_%i",pTHatMin[k],pTHatMax[k]);
      GE->SetName(MyGraphName);
      GE->Write();      
    }

  //Final Plot TotalJetCrossSection with error
  cout<<"Final results for total crossSection"<<endl;
  cout<<"JetpT \t"<<"TotalCrossSection \t"<<"Error"<<endl;
  foutput<<"#JetpT \t"<<"DifferentialCrossSection \t"<<"JetpTBinWidth/2.0 \t"<<"DCSError \t"<<"DCSErr%"<<endl;
  for(int j=0;j<NpTJetBin;j++)
    {
      TotalCrossSectionError[j] =TMath::Sqrt(TotalCrossSectionError[j]);
      cout<<pT[j]<<"\t"<<TotalCrossSection[j]<<" +/- \t"<<TotalCrossSectionError[j]<<endl;
      foutput<<pT[j]<<"\t"<<TotalCrossSection[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionError[j]<<"\t"<<TotalCrossSectionError[j]*100/TotalCrossSection[j]<<endl;
    }
  cout<<"Create Final TGraphError plot "<<endl;
  //Final Plot and save also
  TGraphErrors *GETotal = new TGraphErrors(NpTJetBin,pT,TotalCrossSection,pTError,TotalCrossSectionError);
  GETotal->SetName("TotalJetCrossSectionGraphError");
  GETotal->Write();
  foutput.close();

  // Deleting outFile.
  delete HistTemp;
  delete HistTempScaled;
  //delete HistFinal;
  delete GE;
  delete GETotal;
  delete GHardCrossSection;  
  //Done
  int EndTime = time(NULL);
  int Hour = (EndTime-StartTime)/3600;
  int Minute = ((EndTime-StartTime)/60)-Hour*60;
  int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
  cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
        
  return 0;
}
