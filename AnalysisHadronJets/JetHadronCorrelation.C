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
   sprintf(InputDataFileDIR,"%s",argv[1]);
   sprintf(OutputDataFileName,"%s",argv[2]);
   sprintf(JetscapeORPythia,"%s",argv[12]);
   sprintf(HadronORParton,"%s",argv[13]);
   
   string ExpName=string(argv[3]);
   int Ecm = atoi(argv[4]);
   double JetRadius = atof(argv[5])/1.0;                                                             
   double EtaJetMin = atof(argv[6])/1.0;
   double EtaJetMax=atof(argv[7])/1.0;                                                        
   double EtaTrackMax=atof(argv[8])/1.0;
   double pTTrackMin= atof(argv[9])/1.0;
   int EtaCutFlag =atoi(argv[10]); //0 or 1; 0 to use rapididty cut, 1 to use eta cut
   int ChargedFlag = atoi(argv[11]); // 0 or 1; 0 means all particle, 1 means only charged
   int BGS = atoi(argv[14]);  //0 pp, 1 pbpb
   string Centrality = string(argv[15]);
   string PPorPbPb=string(argv[16]); //  Species= pp or PbPb pPb  
   int EtaRangeType = atoi(argv[17]); //2 for both sided, 1 for one sided
   double JetpTMin = 10; //in GeV
   double AliceTrackpT=0.0;// in GeV jet with atleast single track with pT or above
   double STARTrackpT=0.0; // in GeV jet with atleast single track with pT or above
   double InelasticCS = 70.0; //
   double PI=3.1415926;
   //int pTHatMin[2]={200,210}, pTHatMax[2]={210,220}, NpTHardBin=2;
   //int pTHatMin[2]={500,550}, pTHatMax[2]={550,600}, NpTHardBin=2;
    
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
       //int Bin1[32] = {  50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,    			150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550};
       //int Bin2[32] = {  55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,			160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600};

       //int Bin1[50] = {  13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,                       150, 160, 170, 180, 190,			 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900,  1000, 1100, 1200, 1300, 1400};
       //int Bin2[50] = {  15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,                      160, 170, 180, 190, 200,			 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500};

       int Bin1[62] = { 5, 7, 9,  11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,                                150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600,                                  700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400};
        int Bin2[62] = { 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,                               160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700,                                  800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510};       

       NpTHardBin = 62; //32 41 62 66
       InelasticCS = 70.0;
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


  //5.02TeV pp

  if(  Ecm==5020 && ( PPorPbPb=="pp" || PPorPbPb=="ppPbPb" || PPorPbPb=="PbPb") && ExpName=="Alice" && Centrality=="0-10" &&  (JetRadius==0.2 || JetRadius==0.4 || JetRadius==0.5 ) && BGS>=0)
    {
      double Bin[11]  ={10, 15, 20, 25, 30, 40, 50, 60,  80, 110, 140};
      NpTJetBin = 11-1;     JetpTMin = 5;   AliceTrackpT=5.0;
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
  double dNdpTCount[NpTHardBin][NpTJetBin];  //[ptHatBin] [Regular pt]
  double pTHardBinJetBinError[NpTHardBin][NpTJetBin];
  double TotalCrossSection[NpTJetBin];
  double TotalCrossSectionError[NpTJetBin];
  int Events[NpTHardBin];
  int NetJetEvents[NpTHardBin];
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
  TH1D *HistTemp = new TH1D("TriggerNormalizedJetSpectrumBin", "Triggered Jet Spectrum pT", NpTJetBin, JetpTBin); //5GeV Bin CountVspT
  TH1D *HistTempSignal = new TH1D("Trigger1NormalizedJetSpectrumBin", "Triggered Jet Spectrum pT", NpTJetBin, JetpTBin); 
  TH1D *HistTempRef = new TH1D("Trigger2NormalizedJetSpectrumBin", "Triggered Jet Spectrum pT", NpTJetBin, JetpTBin);
  TH1D *HistTotalSignal = new TH1D("TotalTrigger1NormalizedJetSpectrumBin", "Triggered Jet Spectrum pT", NpTJetBin, JetpTBin);
  TH1D *HistTotalRef = new TH1D("TotalTrigger2NormalizedJetSpectrumBin", "Triggered Jet Spectrum pT", NpTJetBin, JetpTBin);      HistTotalSignal->Reset();      HistTotalRef->Reset();
  HistTemp->Sumw2(); HistTempSignal->Sumw2();HistTempRef->Sumw2();HistTotalSignal->Sumw2();HistTotalRef->Sumw2();
  ofstream EventR, EventA, foutput; char VarFileName[100000];
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFiles/EventRecord%s.txt",OutputDataFileName);
  EventR.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFiles/CodeCheck%s.txt",OutputDataFileName);
  EventA.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFiles/%s.root",OutputDataFileName); // name of the output root file               
  TFile* outFile = new TFile( VarFileName, "RECREATE");
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFiles/%s.txt",OutputDataFileName);
  foutput.open(VarFileName,ios::out);
  cout<<"Final result will be store in  File="<<VarFileName<<endl;
  EventA<<"OutputFile Name is "<<VarFileName<<endl;
  EventR<<"pTHatBin \t "<<"TotalEventsRegistered"<<endl;
  int NRefHadron=0,NSignalHadron=0; double PhiJet=0,PhiHadron=0,DeltaPhi=0;
  std::vector<double> RecoilParticleE, RecoilParticlePx, RecoilParticlePy, RecoilParticlePz, RecoilParticleEta, RecoilParticlePhi;
  std::vector <fjcore::PseudoJet> fjInputs, TriggerSignalList,TriggerRefList;
  Pythia pythia;    

  // For loop to open different pTHat bin files
  for (int k = 0; k<NpTHardBin; ++k)
    {
      char HadronFile[100000], HardCrossSectionFile[100000], pTBinString[10000];
      sprintf(HadronFile,"%s/%s%sListBin%i_%i.out",InputDataFileDIR,JetscapeORPythia,HadronORParton, pTHatMin[k],pTHatMax[k]);
      if(PPorPbPb == "pppPb" || PPorPbPb == "pPb")
	{
	  sprintf(HardCrossSectionFile,"%s/SigmaHardBin%i_%i.out", InputDataFileDIR, pTHatMin[k], pTHatMax[k]);
	}
      else
	{
	  sprintf(HardCrossSectionFile,"/wsu/home/groups/maj-shen/AAPaperData/pTHatCrossSection_%iGeV/SigmaHardBin%i_%i.out",  Ecm, pTHatMin[k],pTHatMax[k]);
	}

      
      sprintf(pTBinString,"Current pTHatBin is %i (%i,%i) GeV",k,pTHatMin[k],pTHatMax[k]);
      EventA<<pTBinString<<endl;
      EventA<<"Hadron File ="<<HadronFile<<endl;
      EventA<<"HardCrossSectionFile="<<HardCrossSectionFile<<endl;
      ifstream myfile, myfile2;
      myfile.open(HadronFile,ios::in);
      myfile2.open(HardCrossSectionFile,ios::in);

     int  NewEvent=0, SN=0,PID;
     string x, y,z,t;
     double Px, Py, Pz, E, Eta, Phi, pStat;
     
     // Reset
     char HistName[1000];
     HistTemp->Reset(); 
     sprintf(HistName,"CountVspTSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]); HistTemp->SetName(HistName);
     HistTempSignal->Reset();
     sprintf(HistName,"SignalVspTSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]); HistTempSignal->SetName(HistName);
     HistTempRef->Reset();
     sprintf(HistName,"RefVspTSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]); HistTempRef->SetName(HistName);

     RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
     RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
     fjInputs.resize(0); TriggerSignalList.resize(0); TriggerRefList.resize(0);

     // Parameters for the jet finders. 
     fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius);

// Read file
     HashFlag=0;
     Events[k] =0;
     NetJetEvents[k] = 0;
     string EventLabel, MyString("#");
     while ( myfile >> EventLabel && Events[k] < 500 ) 
       {
	 if( MyString.compare(EventLabel)==0 && fjInputs.size()>0 )
	   {
	     NetJetEvents[k] = NetJetEvents[k] +1;	     
		 // Run Fastjet algorithm and sort jets in pT order.
		 vector <fjcore::PseudoJet> UnSortedJets, SortedJets;
		 fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
		 UnSortedJets = clustSeq.inclusive_jets(JetpTMin);
		 SortedJets    = sorted_by_pt(UnSortedJets);

		 // List first few FastJet jets and some info about them.
		 if (nListJets) 
		   {
		     cout<< "\n ---FastJet jets before recoil substraction, R = ";
		     cout<< JetRadius << "\t anti-kt with pTHatBin = "<<k<<"=("<<pTHatMin[k]<<","<<pTHatMax[k]<<")";
		     cout<<"\t JetpTMin="<<JetpTMin<<"\n"<<endl;
		     cout<< "  i         pT        E        y        eta          phi_std         phi" << endl;
		     for (int i = 0; i < int(SortedJets.size()); ++i) 
		       {
			 vector<fjcore::PseudoJet> constituents = SortedJets[i].constituents();			 
			 cout << setw(4) << i << fixed << setprecision(3) << setw(9)
			      << SortedJets[i].perp() << setw(9)  << SortedJets[i].e()<< setw(9)  << SortedJets[i].rap()
			      << setw(9) << SortedJets[i].eta()<<setw(9)<< SortedJets[i].phi_std()<<setw(9) << SortedJets[i].phi()<< endl;
		       }		       
		     cout << "\n --------  End FastJet Listing  ------------------"<<endl;

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
                         cout << "\n --  FastJet jets after Subtracting recoil, R = ";
			 cout << JetRadius << "\t anti-kt for pTHatBin = "<<k<<"=("<<pTHatMin[k]<<","<<pTHatMax[k]<<")......\n";
                         cout << "  i         pT           y        eta         phi_std      phi  " << endl;
                         for (int i = 0; i < int(SortedJets.size()); ++i)
                           {
                             cout << setw(4) << i << fixed << setprecision(3) << setw(11);
			     cout << SortedJets[i].perp() << setw(9)  << SortedJets[i].rap();
			     cout << setw(9) << SortedJets[i].eta() << setw(9)  << SortedJets[i].phi_std() << setw(9)  << SortedJets[i].phi()<< endl;
                           }
                         cout << "\n --------  End FastJet Listing  ------------------";
			 cout << "---------------------------------" << endl;
                       }
		   }

		 for (int i = 0; i < pFast; ++i)
		   {
		     if( EtaCutFlag==1){JetEta = SortedJets[i].eta(); }
		     if( EtaCutFlag==0){JetEta = SortedJets[i].rap(); }    		     			 
		     NRefHadron=0; NSignalHadron=0;
			 if( EtaJetMin <= fabs(JetEta) && fabs(JetEta) < EtaJetMax)       
			   {
			     PhiJet = SortedJets[i].phi();
			     for(int j=0; j<TriggerSignalList.size(); j++)
			       {
                                 PhiHadron = TriggerSignalList[j].phi();
				 DeltaPhi = PhiHadron - PhiJet;
				 if (fabs(DeltaPhi > PI)){ DeltaPhi=2*PI - fabs(DeltaPhi); }
				 if(fabs( fabs(DeltaPhi) - PI )<0.6 ) {NSignalHadron = NSignalHadron+1;} 	   
			       }
			     for(int j=0; j<TriggerRefList.size(); j++)
			       {
				 PhiHadron = TriggerRefList[j].phi();
				 DeltaPhi = PhiHadron - PhiJet;
				 if (fabs(DeltaPhi > PI)){ DeltaPhi=2*PI - fabs(DeltaPhi); }
                                 if(fabs( fabs(DeltaPhi) - PI )<0.6 ) {NRefHadron = NRefHadron+1;}
			       }

			     cout<<"NSignalHadron and Ref="<<NSignalHadron<<","<<NRefHadron<<endl;
			     cout<<"Total Signal, Ref="<<TriggerSignalList.size()<<","<<TriggerRefList.size()<<endl;
			     if(TriggerSignalList.size() > 0) 
			       { HistTemp->Fill( SortedJets[i].perp(), 1.0*NSignalHadron/TriggerSignalList.size());
				 HistTempSignal->Fill(SortedJets[i].perp(),1.0*NSignalHadron/TriggerSignalList.size());
			       }
			     if(TriggerRefList.size() > 0) 
			       { HistTemp->Fill( SortedJets[i].perp(), -1.0*NRefHadron/TriggerRefList.size()); 
				 HistTempRef->Fill( SortedJets[i].perp(), 1.0*NRefHadron/TriggerRefList.size());
			       }
			   }			 
		   }	 
		 RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
		 RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
		 fjInputs.resize(0);TriggerSignalList.resize(0); TriggerRefList.resize(0);
		 NewEvent = 0;
		 // cout<<"Found a Jet \t "<<pTBinString<<"\t NetJetevents is \t"<<NetJetEvents[k]<<endl;
	   }
	 
	 if(MyString.compare(EventLabel)!=0)
	   {
	     myfile >> PID >> pStat >> E >> Px >> Py >> Pz >> Eta >> Phi;
	     cout<<" "<<PID<<" "<< pStat<< " " << E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi<<endl;
	     double PT = TMath::Sqrt( pow(Px ,2.0) + pow(Py , 2.0) );
	     //pStat=0;//Uncomment  THIS WHEN RUNNING PbPb heavy ION DATA
	     int Charge = pythia.particleData.charge( PID );

                 if( fabs(Eta) < EtaTrackMax   && Charge!=0 && pStat!=-1 && PT > pTTrackMin  && PID!=12 && PID!=14 && PID!=16 && PID!=18){
		   fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
		   fjInputs[int(fjInputs.size())-1].set_user_index(PID);
		   if(20.0 <= PT && PT <= 50.0 ) {TriggerSignalList.push_back(fjcore::PseudoJet(Px,Py,Pz,E));}
		   if(5.0 <= PT && PT <= 7.0 ) {TriggerRefList.push_back(fjcore::PseudoJet(Px,Py,Pz,E));}
		 }
		 if( fabs(Eta) < EtaTrackMax   && Charge!=0 && pStat==-1 && PT > pTTrackMin)
                   {
                     RecoilParticleE.push_back(E);   RecoilParticlePx.push_back(Px);   RecoilParticlePy.push_back(Py);
                     RecoilParticlePz.push_back(Pz); RecoilParticleEta.push_back(Eta); RecoilParticlePhi.push_back(Phi);
                   }

	     if(HashFlag==1 )
	       {
                 HashFlag=0;
                 Events[k]= Events[k] +1;
               }
	   }

	 if(MyString.compare(EventLabel)==0)
           {	    
	     for(int i=0;i<8;i++){myfile >> EventLabel; cout<<EventLabel<<" ";
	     }
	     HashFlag=1;
	     cout<<endl;
           }
	 
       }
     
    
     //scalling and adding temp histograms     
     myfile2 >> HardCrossSection[k]>>HardCrossSectionError[k];
     cout<<std::fixed<<std::setprecision(10);
     cout<<"Sigma = "<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl; 
     EventA<<"HardCrossSection="<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl;
     if(Error == 1)
       {
	 HardCrossSectionError[k]=0.0;
       }

     cout<<"For pT-hatBin k="<<k<<endl;
     for(int j=0;j<NpTJetBin;j++)
       {	 
         dNdpTCount[k][j]= HistTemp->GetBinContent(j+1);
	 if(dNdpTCount[k][j]> 0.0)
	   {
	     pTHardBinJetBinError[k][j] = (dNdpTCount[k][j]*HardCrossSection[k]/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet)))*TMath::Sqrt( (1/dNdpTCount[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*JetpTBinWidth*2.0*(EtaJetMax-EtaJetMin))<<endl;
	   }
	 else
	   {
	     pTHardBinJetBinError[k][j] = 0.0;
	     dNdpTCount[k][j]=0.0;
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
	   }
	 cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet))<<endl;
	 EventR<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet))<<endl;
       }     
     HistTemp->Write();      
     HistTempSignal->Scale(HardCrossSection[k]/(2.0*(EtaJetMax-EtaJetMin)*Events[k]*InelasticCS));
     HistTempRef->Scale(HardCrossSection[k]/(2.0*(EtaJetMax-EtaJetMin)*Events[k]*InelasticCS));
     HistTotalSignal->Add(HistTempSignal);
     HistTotalRef->Add(HistTempRef);
     HistTempSignal->Write(); HistTempRef->Write();
     EventR<<"For bin k="<<k<<" pTHat=("<<pTHatMin[k]<<"-"<<pTHatMax[k]<<")GeV, Events= \t"<<Events[k]<<", \t"<<NetJetEvents[k]<<endl;
     myfile.close();
     myfile2.close();

// End file for k-loop
   }
  
  EventR.close();
  HistTotalSignal->Scale(1.0,"width");
  HistTotalRef->Scale(1.0,"width");
  HistTotalSignal->Write(); HistTotalRef->Write();
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
  TGraphErrors *GE;
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
	  CrossSection[j] = (dNdpTCount[k][j]*HardCrossSection[k])/(InelasticCS*Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
	  CrossSectionError[j] = pTHardBinJetBinError[k][j]/InelasticCS;
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
  foutput<<"#JetpT \t"<<"TotalCrossSection \t"<<"JetpTBinWidth/2.0 \t"<<"TCSError \t"<<"TCSErr% \t Signal \t Ref \t Signal-Ref"<<endl;
  double ans=0;
  for(int j=0;j<NpTJetBin;j++)
    {
      TotalCrossSectionError[j] =TMath::Sqrt(TotalCrossSectionError[j]);
      ans = HistTotalSignal->GetBinContent(j+1) - HistTotalRef->GetBinContent(j+1);
      cout<<pT[j]<<"\t"<<TotalCrossSection[j]<<" +/- \t"<<TotalCrossSectionError[j]<<endl;
      foutput<<pT[j]<<"\t"<<TotalCrossSection[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionError[j]<<"\t"<<TotalCrossSectionError[j]*100/TotalCrossSection[j]<<"\t"<<HistTotalSignal->GetBinContent(j+1)<<"\t"<<HistTotalRef->GetBinContent(j+1)<<"\t"<<ans<<endl;
    }
  cout<<"Create Final TGraphError plot "<<endl;
  //Final Plot and save also
  TGraphErrors *GETotal = new TGraphErrors(NpTJetBin,pT,TotalCrossSection,pTError,TotalCrossSectionError);
  GETotal->SetName("TotalJetCrossSectionGraphError");
  GETotal->Write();
  foutput.close();

  // Deleting Dynamic pointers
  delete HistTemp;
  delete GE;
  delete GETotal;
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
