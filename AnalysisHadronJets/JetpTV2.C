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
// [1] InputDataFileDIR [2] OutputDataFileName  [3]=CMS or Atlas or Alice or STAR, [4]= 200, 2760 or 5020, 7000 [5]=JetRadius, [6]=EtaJetMin, [7]=EtaJetMax    [8] EtaTrackMax [9] pTTrackMin [10] EtaCutFlag =eta or y [11] ChargedFlag 0 or 1 (all or charged)  [12]=JetscapeORPythia, [13] HadronORParton [14] BGS ( 0 or 1 pp, pbpb)  [15] Centrality [16] Species= pp or PbPb pPb [17]  SoftV2Flag=0 or 1 (off, on)
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
   int SoftV2Flag = atoi(argv[17]);  //0= without soft v2, 1=with soft v2 

   double JetpTMin = 10; //in GeV
   double dEtaJet=2.0*(EtaJetMax -EtaJetMin);
   double AliceTrackpT=0.0;// in GeV jet with atleast single track with pT or above
   int *GoodJetStatus;//1 is fill, 0 not fill
   //int pTHatMin[1]={500}, pTHatMax[1]={550}, NpTHardBin=1;

   
   int *pTHatMin;
   int *pTHatMax;
   int NpTHardBin;
    
   if(  Ecm == 200 )
     {
       int Bin1[23] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90};
       int Bin2[23] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100};
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
       int Bin2[54] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130,                             140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450,                               500, 550, 600, 700, 800, 900, 1000, 1380 };
       NpTHardBin = 54;
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
         {
           pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
         }
     }

   if(  Ecm == 5020 )
     {  // 66= 1, 2, 3, 4, 5
       int Bin1[66] = { 1, 2, 3, 4, 5, 7, 9,  11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,                            160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800,                             900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400};
       int Bin2[66] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160,                           170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900,                             1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510};
       NpTHardBin = 66;
       pTHatMin = new int[NpTHardBin];
       pTHatMax = new int[NpTHardBin];
       for(int i=0; i<NpTHardBin; i++)
         {
           pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
         }
     }  
   

  double *JetpTBin;
  int NpTJetBin;


  double Bin[11]  ={50, 60, 71, 79, 89, 100, 126, 158, 200, 251, 300};
  NpTJetBin = 11-1;     JetpTMin = 20;
  JetpTBin = new double[NpTJetBin+1];
  for(int i=0; i<NpTJetBin+1; i++)
    {
      JetpTBin[i] = Bin[i];    
    }      


  int Error =0; //set to 1 if wants statistical error only
  //double JetpTBinWidth = (JetpTMax - JetpTMin)/NpTJetBin;
  double HardCrossSection[NpTHardBin];
  double HardCrossSectionError[NpTHardBin];
  double dNdpTCount[NpTHardBin][NpTJetBin];  //[ptHatBin] [Regular pt]
  double dNCosinedpTCount[NpTHardBin][NpTJetBin];
  double SumSquareCosineWeight[NpTHardBin][NpTJetBin];
  double SumSquareSoftV2[NpTHardBin];
  double pTHardBinJetBinError[NpTHardBin][NpTJetBin];
  double TotalCrossSection[NpTJetBin];
  double TotalCrossSectionError[NpTJetBin];
  int Events[NpTHardBin];
  int NetJetEvents[NpTHardBin];
  int HashFlag;
// Histograms.
  TH1D *HistTemp = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //5GeV Bin CountVspT
  TH1D *HistTempV2 = new TH1D("JetV2SpectrumBin", "Jet v2 Spectrum pT", NpTJetBin, JetpTBin); //5GeV Bin CountVspT
  TH1D *HistTempV2SingleEvent = new TH1D("JetV2SpectrumSingleEvent", "Jet V2 Spectrum pT for single event", NpTJetBin, JetpTBin); 
  ofstream EventR, EventA, foutput, foutputJetV2; char VarFileName[100000];
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFilesV2/EventRecord%s.txt",OutputDataFileName);
  EventR.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFilesV2/CodeCheck%s.txt",OutputDataFileName);
  EventA.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFilesV2/%s.root",OutputDataFileName); // name of the output root file               
  TFile* outFile = new TFile( VarFileName, "RECREATE");
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFilesV2/InclusiveJCS%s.txt",OutputDataFileName);
  foutput.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/JetFilesV2/%s.txt",OutputDataFileName);
  foutputJetV2.open(VarFileName,ios::out);
  cout<<"Final result will be store in  File="<<VarFileName<<endl;
  EventA<<"OutputFile Name is "<<VarFileName<<endl;
  EventR<<"pTHatBin \t "<<"TotalEventsRegistered"<<endl;

  std::vector<double> RecoilParticleE, RecoilParticlePx, RecoilParticlePy, RecoilParticlePz, RecoilParticleEta, RecoilParticlePhi;
  std::vector <fjcore::PseudoJet> fjInputs;
  Pythia pythia;    

  // For loop to open different pTHat bin files
  for (int k = 0; k<NpTHardBin; ++k)
    {
      char HadronFile[100000], HardCrossSectionFile[100000], pTBinString[10000], SoftV2Coefficient[100000];
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
      ifstream myfile, myfile2, myfile3;
      myfile.open(HadronFile,ios::in);
      myfile2.open(HardCrossSectionFile,ios::in);

     int  NewEvent=0, SN=0,PID;
     string x, y,z,t; double EventPlaneAngle=0;
     double Px, Py, Pz, E, Eta, Phi, pStat;
     double soft_v2, soft_v2_error, psi2, psi2_error, soft_v3, soft_v3_error, psi3, psi3_error;
     int EventsPerFile=0;
     // Reset
     char HistName[1000];
     for(int j=0;j<NpTJetBin;j++)
       {
	 SumSquareCosineWeight[k][j] = 0;
       }
     SumSquareSoftV2[k] =0;

     HistTemp->Reset(); 
     sprintf(HistName,"CountVspTSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
     HistTemp->SetName(HistName);
     HistTempV2->Reset();
     sprintf(HistName,"CosineCountVspTSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
     HistTempV2->SetName(HistName);
     HistTempV2SingleEvent->Reset();
     RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
     RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
     fjInputs.resize(0);

     // Parameters for the jet finders. 
     fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius);

// Read file
     HashFlag=0;
     Events[k] =0;
     NetJetEvents[k] = 0;
     string EventLabel, MyString("#"), SoftV2Label, DummyS;

     if(SoftV2Flag==1)
       {
	 sprintf(SoftV2Coefficient,"/wsu/home/groups/maj-shen/AAPaperData/MATTER_LBT_RunningAlphaS_Q2qhat/SoftAnisotropy_Coefficients_5TeV_%s.txt",Centrality.c_str());
	 myfile3.open(SoftV2Coefficient,ios::in);
	 myfile3>> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label;
	 EventA<<"Soft v2 File ="<<SoftV2Coefficient<<endl;
       }

     while ( myfile >> EventLabel   ) 
       {
	 if( MyString.compare(EventLabel)==0 && fjInputs.size()>0 )
	   {
	     NetJetEvents[k] = NetJetEvents[k] +1;	     
		 // Run Fastjet algorithm and sort jets in pT order.
		 vector <fjcore::PseudoJet> UnSortedJets, SortedJets;
		 fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
		 UnSortedJets = clustSeq.inclusive_jets(JetpTMin);
		 SortedJets    = sorted_by_pt(UnSortedJets);
		 int *GoodJetStatus; GoodJetStatus = new int[int(SortedJets.size())];

		 for (int i = 0; i < int(SortedJets.size()); i++)
		   {GoodJetStatus[i]=1;}
		 
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
		     /*
		     if(BGS==1)		     
		       {cout << "\n -------- List of Negative/recoil particles to be subtracted ------------------\n";
			cout << "i       E      pT     Eta      Phi"<<endl;
			 for (int j = 0; j < int(RecoilParticleE.size()); ++j)
			   {
			     cout << setw(4) << j << fixed << setprecision(3) << setw(9) << RecoilParticleE[j];
			     cout << setw(9)  << TMath::Sqrt( pow(RecoilParticlePx[j],2.0) + pow(RecoilParticlePy[j],2.0) );
			     cout << setw(9) << RecoilParticleEta[j] << setw(9)  << RecoilParticlePhi[j] << endl;
			   }
			 cout << "---------------------------------" << endl;
		       }
		     */
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
		     ///if(GoodJetStatus[i]==1)
		       {
			 double COSINE = soft_v2*std::cos(2.0*(SortedJets[i].phi_std() - EventPlaneAngle));
			 if( EtaCutFlag==1 && EtaJetMin <= fabs(SortedJets[i].eta()) && fabs(SortedJets[i].eta()) < EtaJetMax )
		       {
			 HistTemp->Fill( SortedJets[i].perp() );
			 HistTempV2->Fill( SortedJets[i].perp(), COSINE);
			 HistTempV2SingleEvent->Fill( SortedJets[i].perp(), COSINE);
		       }
			 if( EtaCutFlag==0 && EtaJetMin <= fabs(SortedJets[i].rap()) && fabs(SortedJets[i].rap()) < EtaJetMax )
                       {
                         HistTemp->Fill( SortedJets[i].perp() );
			 HistTempV2->Fill( SortedJets[i].perp(), COSINE);
			 HistTempV2SingleEvent->Fill( SortedJets[i].perp(), COSINE);
                       }
		       }
		   }		 
		 RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
		 RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
		 fjInputs.resize(0);
		 delete [] GoodJetStatus;
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
                 Events[k]= Events[k] +1;
               }
	   }

	 if(MyString.compare(EventLabel)==0)
           { myfile >> EventPlaneAngle;	 cout<<EventPlaneAngle<<" ";   
	     for(int i=0;i<7;i++){myfile >> EventLabel; cout<<EventLabel<<" ";}
	     HashFlag=1;
	     soft_v2 = 1.0;
	     if(SoftV2Flag==1)
	       {
		 myfile3>>SoftV2Label >> soft_v2 >> soft_v2_error >> psi2 >> psi2_error >> DummyS >> DummyS >> DummyS >> DummyS;
		 cout<<SoftV2Label <<" "<< soft_v2 <<" "<< soft_v2_error <<" "<< psi2 <<" "<< psi2_error <<" "<<DummyS  <<" "<< DummyS <<" "<< DummyS <<" "<< DummyS <<" "<<Events[k]+1<<endl;
		 if( int((EventPlaneAngle-psi2)*10000) !=0 )
		   {cout<<"\n Error soft coefficient are mismatched:("<<EventPlaneAngle<<","<<psi2<<"), pTHat= "<<pTHatMin[k]<<", Events="<<Events[k]+1<<endl;
		     EventA<<" \n Error soft coefficent are mismatched:("<<EventPlaneAngle<<","<<psi2<<"), pTHat= "<<pTHatMin[k]<<", Events="<<Events[k]+1<<endl;//cin>>DummyS; 
		     soft_v2 = 1.0;
		   }
	       }
	     for(int j=0;j<NpTJetBin;j++)
	       {
		 SumSquareCosineWeight[k][j] = SumSquareCosineWeight[k][j] + pow(HistTempV2SingleEvent->GetBinContent(j+1),2);
	       }
	     SumSquareSoftV2[k] = SumSquareSoftV2[k] + pow(soft_v2, 2);
	     HistTempV2SingleEvent->Reset();	     
	     cout<<endl;
           }
	 
       }
     
    
     //scalling and adding temp histograms     
     myfile2 >> HardCrossSection[k]>>HardCrossSectionError[k];
     cout<<std::fixed<<std::setprecision(10);
     cout<<"Sigma = "<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl; 
     EventA<<"HardCrossSection="<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl;
     //Last event of the pTHatBin
     for(int j=0;j<NpTJetBin;j++)
       {
	 SumSquareCosineWeight[k][j] = SumSquareCosineWeight[k][j] + pow(HistTemp->GetBinContent(j+1),2);
       }

     if(Error == 1) { HardCrossSectionError[k]=0.0;   }

     for(int j=0;j<NpTJetBin;j++)
       {	 
         dNdpTCount[k][j]= HistTemp->GetBinContent(j+1);
	 dNCosinedpTCount[k][j]= HistTempV2->GetBinContent(j+1);	
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
       }     
     HistTemp->Write();
     HistTempV2->Write(); 
     
     EventR<<"For bin k="<<k<<" pTHat=("<<pTHatMin[k]<<"-"<<pTHatMax[k]<<")GeV, Events= \t"<<Events[k]<<", \t"<<NetJetEvents[k]<<endl;
     myfile.close();
     myfile2.close();
     myfile3.close();
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
	  CrossSection[j] = (dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*(dEtaJet));
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
  foutput<<"#JetpT \t"<<"TotalCrossSection \t"<<"JetpTBinWidth/2.0 \t"<<"TCSError \t"<<"TCSErr%"<<endl;
  for(int j=0;j<NpTJetBin;j++)
    {
      TotalCrossSectionError[j] =TMath::Sqrt(TotalCrossSectionError[j]);
      cout<<pT[j]<<"\t"<<TotalCrossSection[j]<<" +/- \t"<<TotalCrossSectionError[j]<<endl;
      foutput<<pT[j]<<"\t"<<TotalCrossSection[j]<<"\t"<<pTError[j]<<"\t"<<TotalCrossSectionError[j]<<"\t"<<TotalCrossSectionError[j]*100/TotalCrossSection[j]<<endl;
    }

  //Now computing jet V2
  double JetV2[NpTJetBin], JetV2Error[NpTJetBin];
  double ErrorAverageCosineTerm[NpTHardBin][NpTJetBin], AverageSquareV2Soft[NpTHardBin];
  for(int j=0; j<NpTJetBin;j++)
    {
      JetV2[j]=0; JetV2Error[j]=0;
      for(int k=0;k<NpTHardBin;k++)
	{
	  ErrorAverageCosineTerm[k][j] = (SumSquareCosineWeight[k][j] - Events[k]*pow(dNCosinedpTCount[k][j]/Events[k],2) )/(Events[k] -1);
	  ErrorAverageCosineTerm[k][j] = sqrt(ErrorAverageCosineTerm[k][j]/(1.0*Events[k]));
	  AverageSquareV2Soft[k] = sqrt(SumSquareSoftV2[k]/(1.0*Events[k]));

	  JetV2[j] = JetV2[j] + ((dNCosinedpTCount[k][j]*HardCrossSection[k]/(AverageSquareV2Soft[k]*Events[k]*(JetpTBin[j+1]-JetpTBin[j])*1.0*dEtaJet))/TotalCrossSection[j]);
	  JetV2Error[j] = JetV2Error[j] + pow(ErrorAverageCosineTerm[k][j]*HardCrossSection[k]/(AverageSquareV2Soft[k]*TotalCrossSection[j]),2);
	}
      JetV2Error[j] = JetV2Error[j] + pow(JetV2[j]*TotalCrossSectionError[j]/TotalCrossSection[j],2);
      JetV2Error[j] = sqrt( JetV2Error[j]);
    }

  //write to a file
  foutputJetV2<<"#pT \t  JetV2 \t  pTError \t JetV2Error \t Err% " <<endl;
  for(int j=0; j<NpTJetBin;j++)
    {
      cout<< pT[j]<<"\t"<< JetV2[j] <<"\t"<< pTError[j] << "\t" <<JetV2Error[j] <<"\t"<< JetV2Error[j]*100.0/JetV2[j]<<endl;
      foutputJetV2<< pT[j]<<"\t"<< JetV2[j] <<"\t"<< pTError[j] << "\t" <<JetV2Error[j] <<"\t"<< JetV2Error[j]*100.0/JetV2[j]<<endl;
    }

  //save the graph
  TGraphErrors *GEJetV2  = new TGraphErrors(NpTJetBin, pT, JetV2,  pTError, JetV2Error);
  GEJetV2->SetName("JetV2");
  GEJetV2->Write();

  cout<<"Create Final TGraphError plot "<<endl;
  //Final Plot and save also
  TGraphErrors *GETotal = new TGraphErrors(NpTJetBin,pT,TotalCrossSection,pTError,TotalCrossSectionError);
  GETotal->SetName("TotalJetCrossSectionGraphError");
  GETotal->Write();
  foutput.close();
  foutputJetV2.close();

  // Deleting Dynamic pointers
  delete HistTemp, HistTempV2, HistTempV2SingleEvent;
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
