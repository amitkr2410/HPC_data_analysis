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
// [1] InputDataFileDIR [2] OutputDataFileName  [3]=CMS or Atlas or Alice or STAR, [4]= 200, 2760 or 5020, 7000 [5]=JetRadius, [6]=EtaJetMin, [7]=EtaJetMax    [8] EtaTrackMax [9] pTTrackMin [10] EtaCutFlag =eta or y [11] ChargedFlag 0 or 1 (all or charged)  [12]=JetscapeORPythia, [13] HadronORParton [14] BGS ( 0 or 1 pp, pbpb)  

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
   double JetpTMin = 1; //in GeV

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
       int Bin2[54] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130,                             140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450,                               500, 550, 600, 700, 800, 900, 1000, -1 };
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
   

  double *JetpTBin;
  int NpTJetBin;

  if(  Ecm == 2760 && ExpName=="CMS" &&  (JetRadius==0.2  ||  JetRadius==0.3 || JetRadius==0.4) && EtaJetMax==2.0 && BGS>=0)
    {
      double Bin[13]  ={70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
      NpTJetBin = 13-1;     JetpTMin = 40;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
	{
	  JetpTBin[i] = Bin[i];    
	}      
    }

  if(  Ecm == 2760 && ExpName=="CMS" &&  JetRadius==0.7 && EtaJetMax==0.5 && BGS>=0)
    {       
      double Bin[20] = {74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592} ;
      NpTJetBin = 20-1;     JetpTMin = 40; 
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  //5.02TeV pp
  if(  Ecm == 5020 && ExpName=="Atlas" &&  JetRadius==0.4 && (EtaJetMax==0.3 || EtaJetMax==2.8) && BGS<0)
    {
      double Bin[15] = {40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501,631, 800, 1000};
      NpTJetBin = 15-1;     JetpTMin = 5;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }

  if(  Ecm == 5020 && ExpName=="Atlas" &&  JetRadius==0.4 && (EtaJetMin==1.6 && EtaJetMax==2.1) && BGS<0)
    {
      double Bin[13] = {40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501, 631};
      NpTJetBin = 13-1;     JetpTMin = 5;
      JetpTBin = new double[NpTJetBin+1];
      for(int i=0; i<NpTJetBin+1; i++)
        {
          JetpTBin[i] = Bin[i];
        }
    }
  //JetRAA  pp and PbPb
  if(  Ecm == 5020 && ExpName=="Atlas" &&  JetRadius==0.4 && EtaJetMax==2.8 && BGS>=0)
    {
      double Bin[16] = {100, 112,   125,  141, 158, 177, 199, 223, 251, 281, 316, 354, 398, 501, 630, 999};
      NpTJetBin = 16-1;     JetpTMin = 5;
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
// Histograms.
  TH1D *HistTemp = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //5GeV Bin CountVspT

  ofstream EventR, EventA, foutput; char VarFileName[100000];
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/JetFiles/EventRecord%s.txt",OutputDataFileName);
  EventR.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/JetFiles/CodeCheck%s.txt",OutputDataFileName);
  EventA.open(VarFileName,ios::out);
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/JetFiles/%s.root",OutputDataFileName); // name of the output root file               
  TFile* outFile = new TFile( VarFileName, "RECREATE");
  sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/JetFiles/%s.txt",OutputDataFileName);
  foutput.open(VarFileName,ios::out);
  cout<<"Final result will be store in  File="<<VarFileName<<endl;
  EventA<<"OutputFile Name is "<<VarFileName<<endl;
  EventR<<"pTHatBin \t "<<"TotalEventsRegistered"<<endl;

  std::vector<double> RecoilParticleE, RecoilParticlePx, RecoilParticlePy, RecoilParticlePz, RecoilParticleEta, RecoilParticlePhi;
  std::vector <fjcore::PseudoJet> fjInputs;
  Pythia pythia;    

  // For loop to open different pTHat bin files
  for (int k = 0; k<NpTHardBin; ++k)
    {
      char HadronFile[100000], HardCrossSectionFile[100000], pTBinString[10000];
      sprintf(HadronFile,"%s/%s%sListBin%i_%i.out",InputDataFileDIR,JetscapeORPythia,HadronORParton, pTHatMin[k],pTHatMax[k]);
      sprintf(HardCrossSectionFile,"%s/SigmaHardBin%i_%i.out",InputDataFileDIR,  pTHatMin[k],pTHatMax[k]);
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
     sprintf(HistName,"CountVspTSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
     HistTemp->SetName(HistName);
     RecoilParticleE.resize(0);RecoilParticlePx.resize(0);RecoilParticlePy.resize(0);RecoilParticlePz.resize(0);
     RecoilParticleEta.resize(0);RecoilParticlePhi.resize(0);
     fjInputs.resize(0);

     // Parameters for the jet finders. 
     fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius);

// Read file
     HashFlag=0;
     Events[k] =0;
     NetJetEvents[k] = 0;
     string EventLabel, MyString("#");
     while ( myfile >> EventLabel ) 
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
		     if( EtaCutFlag==1 && EtaJetMin <= fabs(SortedJets[i].eta()) && fabs(SortedJets[i].eta()) < EtaJetMax )
		       {
			 HistTemp->Fill( SortedJets[i].perp() );
		       }
		     if( EtaCutFlag==0 && EtaJetMin <= fabs(SortedJets[i].rap()) && fabs(SortedJets[i].rap()) < EtaJetMax )
                       {
                         HistTemp->Fill( SortedJets[i].perp() );
                       }
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
		   {fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));}
		 if( fabs(Eta) < EtaTrackMax  && pStat==-1 && PT > pTTrackMin)
		   {
		     RecoilParticleE.push_back(E);   RecoilParticlePx.push_back(Px);   RecoilParticlePy.push_back(Py);
		     RecoilParticlePz.push_back(Pz); RecoilParticleEta.push_back(Eta); RecoilParticlePhi.push_back(Phi);
		   }
	       }
	     
	     if(ChargedFlag==1)
	       {
                 if( fabs(Eta) < EtaTrackMax   && Charge!=0 && pStat!=-1 && PT > pTTrackMin  && PID!=12 && PID!=14 && PID!=16 && PID!=18){fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));}
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
           {	    
	     for(int i=0;i<8;i++){myfile >> EventLabel; //cout<<EventLabel<<" ";
	     }
	     HashFlag=1;
	     //cout<<endl;
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

     for(int j=0;j<NpTJetBin;j++)
       {	 
         dNdpTCount[k][j]= HistTemp->GetBinContent(j+1);
	 if(dNdpTCount[k][j]> 0.0)
	   {
	     pTHardBinJetBinError[k][j] = (dNdpTCount[k][j]*HardCrossSection[k]/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*(EtaJetMax-EtaJetMin)))*TMath::Sqrt( (1/dNdpTCount[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*JetpTBinWidth*2.0*(EtaJetMax-EtaJetMin))<<endl;
	   }
	 else
	   {
	     pTHardBinJetBinError[k][j] = 0.0;
	     //cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTemp->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
	   }
       }     
     HistTemp->Write();     
     //HistTempScaled->Scale(HardCrossSection[k]/(2.0*(EtaJetMax-EtaJetMin)*Events[k]*JetpTBinWidth));
     //HistTemp->Draw();
     //HistTempScaled->Write();
     //HistFinal->Add(HistTempScaled);
     EventR<<pTHatMin[k]<<"-"<<pTHatMax[k]<<"\t"<<Events[k]<<"\t"<<NetJetEvents[k]<<endl;
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
	  CrossSection[j] = (dNdpTCount[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*(EtaJetMax-EtaJetMin));
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
      
  return 0;
}
