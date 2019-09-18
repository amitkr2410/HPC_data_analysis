//Created by Amit Kumar (kumar.amit@wayne.edu)
// Requires FinalState partons E, Px, Py, Pz as input, also requires pTHatCrossSection (in mb) as input
// Output is a ROOT file which contains following plots

#include <iostream>
// FastJet3 library.
#include "/wsu/home/fy/fy41/fy4125/Software/pythia8230/include/Pythia8/Pythia.h"
#include "/wsu/home/fy/fy41/fy4125/Software/pythia8230/include/Pythia8/FJcore.h"
#include <fstream>
#include <cstdio>
#include "iomanip"
#include "TMath.h"

// ROOT, for histogramming.
#include"TString.h"
#include "TH1.h"
#include "TH2.h"
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
int main(int argc, char* argv[]){
   int StartTime = time(NULL);
// Create the ROOT application environment.
   TApplication theApp("hist", &argc, argv);
   double PI=3.1415926;

// Create file on which histogram(s) can be saved.
   char outFileName[50000];
   int Ecm = atoi(argv[3]);

     //[1]=Jetscape or Pythia, [2]=CMS or Atlas or Alice or Star, [3]= 200 or 2760 or 7000, [4]EtaCutFlag, [5]=ChargedFlag=0 or 1 (all or charged), [6]=Hadron or Parton, [7]=EtaTrackMax, [8]=pTTrackMin, [9]=pTTrigMin [10]=pTTrigMax [11]=pTAssMin  [12]=pTAssMax [13]=FileNameOutPut  [14]=HadronDirectoryRun  [15] NumberOfFiles

  string ExpName=std::string(argv[2]);
  int EtaCutFlag =atoi(argv[4]); //0 or 1; 0 to use rapididty cut, 1 to use eta cut
  int ChargedFlag = atoi(argv[5]); // 0 or 1; 0 means all particle, 1 means only charged  
  double EtaTrackMax=atof(argv[7])/1.0;
  double pTTrackMin= atof(argv[8])/1.0;
  double pTTrigMin = atof(argv[9])/1.0;
  double pTTrigMax = atof(argv[10])/1.0;
  double pTAssMin = atof(argv[11])/1.0;
  double pTAssMax = atof(argv[12])/1.0;
  int NumberOfFiles = atoi(argv[15]);
  // int Ecm = atoi(argv[3]); //Center-of-mass energy of PP collision

  double DeltaPhiMin=-PI/2.0; double DeltaPhiMax=3.0*PI/2.0;
  double DeltaEtaMin=-2.0; double DeltaEtaMax=2.0;
  double *DeltaPhiBin;  int NDeltaPhiBin=50;
  double *DeltaEtaBin;  int NDeltaEtaBin=50;
  cout<<argv[12]<<endl;
  cout<<argv[13]<<endl;
  cout<<argv[14]<<endl;
  sprintf(outFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/%s.root",argv[13]);
  TFile* outFile = new TFile(outFileName, "RECREATE");
  cout<<"ROOT FIle = "<<outFileName<<endl;

  int BGS=0; //for recoil substraction set it to 1
  int Error =0; //set to 1 if wants statistical error only
  double dNdDeltaEtaDeltaPhiCount[NDeltaPhiBin][NDeltaEtaBin];  //[ptHatBin] [Regular pt]
  double pTHardDeltaPhiDeltaEtaBinError[NDeltaPhiBin][NDeltaEtaBin];
  double Normalized2PSpectra[NDeltaPhiBin][NDeltaEtaBin];
  double Normalized2PSpectraError[NDeltaPhiBin][NDeltaEtaBin];

// Histograms.
  TH2D *HistSignal = new TH2D("CountVsDeltaEtaDeltaPhiSpectrum", "2Particle Spectrum pT", NDeltaEtaBin, DeltaEtaMin, DeltaEtaMax, NDeltaPhiBin, DeltaPhiMin, DeltaPhiMax); //0.25 Bin size 
  TH2D *HistBG = new TH2D("CountVsDeltaEtaDeltaPhiBG", "2Particle background", NDeltaEtaBin, DeltaEtaMin, DeltaEtaMax, NDeltaPhiBin, DeltaPhiMin, DeltaPhiMax); // 0.25 Bin size
  
  
  ofstream EventR, foutput;
  char EventFile[50000];sprintf(EventFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/EventRecord%s.txt",argv[13]);
  EventR.open(EventFile,ios::out);
  char DataFile[50000];
  sprintf(DataFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/%s.txt",argv[13]);
  foutput.open(DataFile,ios::out);
  EventR<<"pTHatBin \t "<<"TotalEventsRegistered"<<endl;
  
  std::vector<double> EtaVector, PhiVector, pTVector;
  double MaxpT=0,MaxpTIndex=0;
  int Events=0;
  int NetEvents=0;
  Pythia pythia;    

  // For loop to open different pTHat bin files
  for (int k = 0; k<NumberOfFiles; ++k)
    {
      char PartonHadronFile[30000];
      sprintf(PartonHadronFile,"/wsu/home/fy/fy41/fy4125/RUN/%s_%d.dat",argv[14],k);

      cout<<"PartonHadronFile Name = "<<PartonHadronFile<<"\n"<<argv[14]<<endl;
      ifstream myfile, myfile2;
      myfile.open(PartonHadronFile,ios::in);
     int  NewEvent=0, SN=0,PID=0;
     string x, y,z,t;
     double Px, Py, Pz, E,pStat;
     MaxpT=0,MaxpTIndex=0;
     // Reset
     char HistName[100];
     EtaVector.resize(0);PhiVector.resize(0);pTVector.resize(0);
// Read file    
     NetEvents =0;
     string EventLabel, MyString("#");
     myfile >>EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>>EventLabel>> EventLabel;

     while ( myfile >> EventLabel ) 
       {
	 //cout<<EventLabel;          	 
	 if( MyString.compare(EventLabel)==0 && EtaVector.size()>0 )
	   { NetEvents = NetEvents +1;	Events = Events +1;     
	     cout<<"For File k="<<k<<"\t event is \t"<<NetEvents<<endl;
	     //pStat=0;//Uncomment  THIS WHEN RUNNING PbPb heavy ION DATA
	     //int Charge = pythia.particleData.charge( PID );
	     for(int i=0; i<pTVector.size(); i++)
               { 
		 for(int j=i+1; j<pTVector.size(); j++)
                   {	
		     double DeltaPhi = PhiVector[i] - PhiVector[j];
		     if(DeltaPhi< -0.5*PI)
		       {
			 DeltaPhi = DeltaPhi + 2.0*PI;
		       }

		     if(DeltaPhi> 1.5*PI)
                       {
                         DeltaPhi = DeltaPhi - 2.0*PI;
                       }

		     HistSignal->Fill( EtaVector[i] - EtaVector[j],  DeltaPhi, 2.0/(pTVector.size()*(pTVector.size()-1.0)));
		     for(int k=0; k<1; k++)
		       {
			 double RandPhi1 = 0.001*(rand() % 1000)*2*PI - (0.5*PI);
			 double RandEta2 = EtaTrackMax - (0.001*(rand() % 1000)*2.0*EtaTrackMax);
			 double RandEta3 = EtaTrackMax - (0.001*(rand() % 1000)*2.0*EtaTrackMax);
			 HistBG->Fill( RandEta2 - RandEta3 , RandPhi1, 2.0/(pTVector.size()*(pTVector.size()-1.0)));
			 HistBG->Fill( RandEta3 - RandEta2 , RandPhi1, 2.0/(pTVector.size()*(pTVector.size()-1.0)));
		       }
		     //cout<<"(Phi, Eta1, Eta2) = ("<<Rand1<<","<<Rand2<<","<<Rand3<<")"<<endl;
		     //cout<<PhiVector[i] - PhiVector[j]<<", "<<EtaVector[i] - EtaVector[j]<<endl;
                   }
               }
	     PhiVector.resize(0); EtaVector.resize(0);pTVector.resize(0); MaxpT=0; MaxpTIndex=0;
	     for(int i=0; i<6; i++){myfile >> EventLabel;}
	   }
         else
	   {
	     myfile >> PID >> pStat >> E >> Px >> Py >> Pz;
	     double pT = sqrt(Px*Px + Py*Py); double ModP = sqrt(Px*Px + Py*Py + Pz*Pz);
	     double Phi = atan2(Py,Px); double Eta = 0.5*log((ModP+Pz)/(ModP-Pz));
	     if(pTTrigMin< pT && pT < pTTrigMax &&  fabs(Eta) < EtaTrackMax ) 
	       {
		 PhiVector.push_back(Phi); EtaVector.push_back(Eta); pTVector.push_back(pT);
	       }
	   }
	 
       }
     cout<<"File="<<k<<", TotalEvents="<<Events<<endl;
     EventR<<k<<","<<NetEvents<<","<<Events<<endl;
     myfile.close();
    }    //End file for k-loop
  //scalling and adding temp histograms     
  HistSignal->Scale(NDeltaEtaBin*NDeltaPhiBin/((DeltaEtaMax-DeltaEtaMin)*(DeltaPhiMax-DeltaPhiMin)*Events));
  HistBG->Scale(NDeltaEtaBin*NDeltaPhiBin/((DeltaEtaMax-DeltaEtaMin)*(DeltaPhiMax-DeltaPhiMin)*Events));
  HistSignal->Write();
  HistBG->Write();
  
  EventR.close();
  foutput.close();
  cout<<"Wrote the Final Histogram in ROOT file "<<endl;
  
  
  // Deleting outFile.
  delete HistSignal;
  delete HistBG;

  //Done
  int EndTime = time(NULL);
  int Hour = (EndTime-StartTime)/3600;
  int Minute = ((EndTime-StartTime)/60)-Hour*60;
  int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
  cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
        
  return 0;
}
