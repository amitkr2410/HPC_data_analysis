// This program reads "test_out.dat" produced from the JETSCAPE and Generates spectrums for jets and charged hadron yield.
// JETSCAPE by default gives pTHatCrossSection in milibarn and particle's momentum in GeV.
// With slight change, this program also allows reading multiple pTHatBins data and produce weighted spectra. 
// Output of this program is a ROOT file which contains following plots
// 1. Number of jets vs pT of the jet. (Graph name is CountVspTJetSpectrumBinpTHatMin_pTHatMax)
// 2. Number of charged hadron vs pT of the hadron. (Graph name is CountVspTSingleHadronSpectrumBinpTHatMin_pTHatMax)
// 3. Weighted differential jet cross section vs pT of the jet. Here Weighted differential crosssection means = sigmapTHat*dN/(dpTJet*dEta*TotalEvents)
//    (Graph name is DifferentialJetCrossSectionBinpTHatMin_pTHatMax)
//
// 4. Weighted charged hadron differential yield vs pT of hadron. Here Weighted differential Yield = sigmapTHat*dN/(TotalInelasticCrossSection*2*PI*dpTHadron*dEta*TotalEvents)  
//    (Graph name is DifferentialSingleHadronYieldBinpTHatMin_pTHatMax)
//
// 5. Total differential jet cross section (i.e. summed over all pTHat bins) vs pT of the jet. (if you have multiple pTHatBins data)
//    (Graph name is TotalDifferentialJetCrossSection)
//
// 6. Total differential charged hadron yield (i.e. summed over all pTHat bins) vs pT of the hadron.
//    (Graph name is TotalDifferentialSingleHadronYield)
//
// 7. pTHat crosssection i.e. hard scattering crosssection vs pTHat bin.
//    (Graph name is HardCrossSection) 
// 
// Note, plots are saved in two ROOT format TH1D and TGraphErrors.
// For jet spectrum, we use anti-kT algorithm with jet radius, R=0.3 and |eta_jet|<2.0. Inside jet cone, we include all particle except neutrinos (CMS def).
// For charged hadron spectrum, we use |eta_hadron|<1.0. Only charged hadrons are included in the spectrum.

// Authorship: written by Amit Kumar

//C++ header
#include "string"
#include <iostream>
#include <fstream>
#include "iomanip"

#include <memory>
#include <chrono>
#include <thread>
#include "/wsu/home/fy/fy41/fy4125/Software/pythia8230/include/Pythia8/Pythia.h"
#include "/wsu/home/fy/fy41/fy4125/Software/pythia8230/include/Pythia8/FJcore.h"

//ROOT headers
#include "TString.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TVector.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"

using namespace std;
using namespace Pythia8;
using namespace fjcore;
int main(int argc, char* argv[]){
  double PI=3.1415926;
    int nListJets =1;  // for output control
    int StartTime = time(NULL);
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);
    
    ////List here all pTHatBins for which you have generated "test_out.dat"
    //int pTHatMin[48] = {1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90,  100,  110, 120, 130, 140, 150, 160, 170, 180, 190,200,                       210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600};
    //int pTHatMax[48] = {2, 3, 4, 5, 10,15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110,  120, 130, 140, 150, 160, 170, 180, 190, 200,210,                       220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 1380};
    //int pTHatMin[69] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 		     170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000,   1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000 };
    // int pTHatMax[69] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170,                      180, 190, 200,  210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000, 3500};

    //By default we want to look at only one pTHatBin spectrum
    int pTHatMin[1]={100};
    int pTHatMax[1]={110};
    int NpTHardBin = sizeof(pTHatMin)/sizeof(pTHatMin[0]); //sizeof(pTHatMin[]); or sizeof(pTHatMax[]); // # of bin we analyze
    double DetectorEtaCut= 3.0;
    int NumberOfFiles=atoi(argv[3]);
    int Events[NpTHardBin];
    double HardCrossSection[NpTHardBin], HardCrossSectionError[NpTHardBin];
    //Variables for jet spectrum
    //(For CMS at 2760 GeV with jet radius R=2,0.3, or 0.4, |eta_jet|<2.0)
    //double JetpTBin[20] = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300}; //in GeV 
    double JetpTBin[26] = {2, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
    int NpTJetBin = 26-1;     double JetpTMin = 0.5; //in GeV
    
    double JetEtaCut = 2.0;
    double JetRadius = 0.7;
    double dNdpTCountJet[NpTHardBin][NpTJetBin];  //[ptHatBin] [Regular pt]  // count for jet in each bin
    double pTHardBinJetBinError[NpTHardBin][NpTJetBin];        // error
    double TotalDifferentialJetCrossSection[NpTJetBin];        // weighted cross section 
    double TotalDifferentialJetCrossSectionError[NpTJetBin];   // weighted error
    
    //Variable for single hadron spectrum
    double SingleHadronpTBin[17] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.2, 10.8, 14.4, 21.6, 28.8, 38.4, 48.0, 67.2, 86.4, 112.2};
    int NpTSingleHadronBin = 17-1;
    
    double SingleHadronEtaCut =3.0;
    double dNdpTCountSingleHadron[NpTHardBin][NpTSingleHadronBin];  //[ptHatBin] [Regular pt]
    double pTHardBinSingleHadronBinError[NpTHardBin][NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYield[NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYieldError[NpTSingleHadronBin];
   
    // for jet substructure
    double JetpTCut = 200;
    //Variables for jet shape
    double JetShaperBin[7] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
    int NrJetShapeBin = 6;
    //Variables for  jet fragmentation function
    double JetFFzBin[11] = {0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1};
    int NzJetFFBin = 10;
    //Variables for jet mass
    double JetMassBin[13] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};
    int NJetMassBin = 12;
    
    //Two-particle correlation
    double pTTrigMin = 1.0;
    double pTTrigMax = 10.0;
    double pTAssMin = 1.0;
    double pTAssMax = 10.0;
    double EtaTrackMax=2.0;
    double DeltaPhiMin=-PI/2.0; double DeltaPhiMax=3.0*PI/2.0;
    double DeltaEtaMin=-2.0; double DeltaEtaMax=2.0;
    double *DeltaPhiBin;  int NDeltaPhiBin=50;
    double *DeltaEtaBin;  int NDeltaEtaBin=50;
    std::vector<double> EtaVector, PhiVector, pTVector;

    // Histograms.1. Number of jets vs pT of the jet. 2. Number of charged hadron vs pT of the hadron
    TH1D *HistTempJet = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //CountVspT for jets
    TH1D *HistTempSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
    TH1D *HistTempJetShape = new TH1D("JetShapeBin", "Jet Shape", NrJetShapeBin, JetShaperBin); //CountVsr for Jet Shape
    TH1D *HistTempJetFF = new TH1D("JetFragmentaionFunctionBin", "Jet Fragmentaion Function", NzJetFFBin, JetFFzBin);
    TH1D *HistTempJetMass = new TH1D("JetMassBin", "Jet Mass", NJetMassBin, JetMassBin);
    TH2D *HistTemp2PCSignal = new TH2D("CountVsDeltaEtaDeltaPhiSpectrum", "2Particle Spectrum pT", NDeltaEtaBin, DeltaEtaMin, DeltaEtaMax, NDeltaPhiBin, DeltaPhiMin, DeltaPhiMax);
    ofstream EventR, EventA;
    char EventFile[50000];sprintf(EventFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/EventRecord%s.txt",argv[2]);
    EventR.open(EventFile,ios::out);
    char DataFile[50000];
    sprintf(DataFile,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/CodeCheck%s.txt",argv[2]);
    EventA.open(DataFile,ios::out);
    char outFileName[10000];
    sprintf(outFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/%s.root",argv[2]); // name of the output root file  
    TFile* outFile = new TFile( outFileName, "RECREATE");

    std::vector <fjcore::PseudoJet> fjInputs;
    std::vector <int> chargeList;
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius); // define jet finding algorithm
    
    Pythia8::Pythia pythia;//("",false);
    int FakeEventCount=0;

    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files, in this default example, only one file
    for (int k = 0; k<NpTHardBin; ++k) //NpTHardBin
    {
        char HadronFile[300], pTBinString[100];        
	//sprintf(HadronFile,"test_out.dat");  // name of the input file        
        sprintf(pTBinString,"Current pTHatBin is %i (%i,%i) GeV",k,pTHatMin[k],pTHatMax[k]);
	EventA<<pTBinString<<endl;
        int  SN=0,PID=0,pStat=0;
        double Px, Py, Pz, E, Eta, Phi;
        Events[k]=0;
	int EventsPerFile=0;
        int TriggeredJetNumber=0;
        double SumE=0, SumPx=0, SumPy=0, SumPz=0;

        // Reset for each pTHardBin
        char HistName[10000];        
        HistTempJet->Reset();
        sprintf(HistName,"CountVspTJetSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJet->SetName(HistName);
        
        HistTempSingleHadron->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadron->SetName(HistName);
        
        HistTempJetShape->Reset();
        sprintf(HistName,"CountVsrJetShapeBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJetShape->SetName(HistName);
        HistTempJetShape->Sumw2();
        
        HistTempJetFF->Reset();
        sprintf(HistName,"CountVszJetFFBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJetFF->SetName(HistName);
        HistTempJetFF->Sumw2();
        
        HistTempJetMass->Reset();
        sprintf(HistName,"CountVsJetMassBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJetMass->SetName(HistName);
        HistTempJetMass->Sumw2();
        
        fjInputs.resize(0);
        chargeList.resize(0);
	EtaVector.resize(0); PhiVector.resize(0); pTVector.resize(0);
	for(int FileIndex=0; FileIndex<NumberOfFiles; FileIndex++)
	  {
	    char HadronFile[30000], HardCrossSectionFile[100000], pTBinString[10000];
	    sprintf(HadronFile,"/wsu/home/groups/maj-shen/%s/JetscapePartonListBin%i_%i.out",argv[1], pTHatMin[k],pTHatMax[k]);  // name of the input file  
	    sprintf(HardCrossSectionFile,"/wsu/home/groups/maj-shen/%s/SigmaHardBin%i_%i.out",argv[1],  pTHatMin[k],pTHatMax[k]);
	    sprintf(pTBinString,"Current Bin is %i pTHatBin(%i,%i) GeV and file is ",k,pTHatMin[k],pTHatMax[k], FileIndex);        
	    EventA<<pTBinString<<endl;
	    EventsPerFile =0;
	    fjInputs.resize(0);chargeList.resize(0);EtaVector.resize(0); PhiVector.resize(0); pTVector.resize(0);
	    // Read file
	    string EventLabel;// MyString("#");
	    string MyString("Event");
	    ifstream myfile;  myfile.open(HadronFile,ios::in);
	    ifstream myfile2; myfile2.open(HardCrossSectionFile,ios::in);
	    //myfile >>EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel;
	    while (myfile >>EventLabel )
	      {
		// Add this particle for Jet spectrum
		if( MyString.compare(EventLabel)==0 && fjInputs.size()>0 )  // construct jet for the previous event
		  {
		    // Run Fastjet algorithm and sort jets in pT order.
		    vector <fjcore::PseudoJet> UnSortedJets, SortedJets;
		    fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
		    UnSortedJets = clustSeq.inclusive_jets(JetpTMin);
		    SortedJets    = sorted_by_pt(UnSortedJets);                    
		    // List first few FastJet jets and some info about them.
		    // if (nListJets && EventsPerFile)
		      {
			cout << "\n --------  FastJet jets, R = " << JetRadius << "anti-kt, pTHatBin="<<k<<", fileIndex="<<FileIndex
			     << "  i         pT        y     eta      phi  " << endl;
			for (int i = 0; i < int(SortedJets.size()); i++)
			  {
			    vector<fjcore::PseudoJet> constituents = SortedJets[i].constituents();
			    cout << setw(4) << i << fixed << setprecision(3) << setw(11)
				 << SortedJets[i].perp() << setw(9)  << SortedJets[i].rap()
				 << setw(9) << SortedJets[i].eta() << setw(9)  << SortedJets[i].phi_std() << endl;
			  }
			cout << "\n --------  End FastJet Listing  ------------------"
			     << "---------------------------------" << endl;
		      }
		    
		    int pFast = SortedJets.size();
		    cout<<"pFast = "<<pFast<<endl;
		    for (int i = 0; i < pFast; i++)  // fill jets from this event into observable bins
		      {
			if(-JetEtaCut < SortedJets[i].eta() && SortedJets[i].eta()< JetEtaCut )
			  {
			    HistTempJet->Fill( SortedJets[i].perp() ); // fill the pT spectrum
			    if(SortedJets[i].perp() > 120.0){FakeEventCount++;}
			    if( SortedJets[i].perp() >= JetpTCut )
			      { // now for substructure				  
				TriggeredJetNumber ++;
				vector<fjcore::PseudoJet> constituents = SortedJets[i].constituents();
				cout<<"For i="<<i<<" jet pT="<<SortedJets[i].perp()<<endl;
				//Jet Shape---
				for( int j = 0; j < constituents.size(); j++ )
				  {				      
				    double delta_eta = constituents[j].eta() - SortedJets[i].eta();
				    double delta_phi = SortedJets[i].delta_phi_to(constituents[j]);
				    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
				    HistTempJetShape->Fill( delta_r, constituents[j].perp() ); // fill the pT vs. r, second variable is the weight to fill in
				    if(SortedJets[i].perp() > 120.0){cout<<"(E,px,py,pz,pt)="<<constituents[j].e() <<","<< constituents[j].px() <<","<<constituents[j].py() <<","<<constituents[j].pz() <<","<<constituents[j].perp() <<"\n"; }
				  }
				
				//Jet Fragmentation Function---
				for( int j = 0; j < fjInputs.size(); j++ )
				  {
				    double delta_eta = fjInputs[j].eta() - SortedJets[i].eta();
				    double delta_phi = SortedJets[i].delta_phi_to(fjInputs[j]);
				    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
				    if( fabs(chargeList[j]) > 0.01 && delta_r <= JetRadius){//charged particle in jet cone
				      double z_jet = fjInputs[j].perp()/SortedJets[i].perp();
				      HistTempJetFF->Fill( z_jet ); // fill in z distribution
				    }
				  }
				
				//Jet Mass---
				double jet_e = 0.0, jet_px = 0.0, jet_py = 0.0, jet_pz = 0.0;
				for( int j = 0; j < constituents.size(); j++ )
				  {
				    double delta_eta = constituents[j].eta() - SortedJets[i].eta();
				    double delta_phi = SortedJets[i].delta_phi_to(constituents[j]);
				    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
				    if( delta_r <= JetRadius){// "ALL" particle in jet cone
				      jet_e += constituents[j].e();
				      jet_px += constituents[j].px();
				      jet_py += constituents[j].py();
				      jet_pz += constituents[j].pz();
				    }
				  }
				double jet_mass = TMath::Sqrt(jet_e*jet_e - jet_px*jet_px - jet_py*jet_py - jet_pz*jet_pz);
				HistTempJetMass->Fill( jet_mass ); // fill in mass distribution
				  
			      }
			  }
		      }
			fjInputs.resize(0);  // reset jet information and put the first particle into a new jet list
			chargeList.resize(0);
			//cout<<"Found a Jet \t "<<pTBinString<<"\t NetJetevents is \t"<<NetJetEvents<<endl;
			// note for Status>=0: here we only analyze positive (jet + recoil) partons in the particle list in Pb-Pb collisions, in principle, one should subtract negative (back-reaction) particle contribution in realistic calcualtions			
		  }
		
		if( MyString.compare(EventLabel)==0 && EtaVector.size()>0 )
		  {
		    for(int n=0; n<pTVector.size(); n++)
		      {for(int m=0; m<pTVector.size(); m++)
			  {double DeltaPhi = PhiVector[n] - PhiVector[m];
			    if(DeltaPhi < -0.5*PI){ DeltaPhi = DeltaPhi + 2.0*PI;}
			    if(DeltaPhi >  1.5*PI){ DeltaPhi = DeltaPhi - 2.0*PI;}
			    if(pTVector[n] > pTVector[m])
			      {HistTemp2PCSignal->Fill( EtaVector[n] - EtaVector[m],  DeltaPhi, 2.0/(pTVector.size()*(pTVector.size()-1.0)));
			      }
			  }
		      }
		    EtaVector.resize(0); PhiVector.resize(0); pTVector.resize(0);
		  }
		
		if(MyString.compare(EventLabel)!=0)
		  { //myfile >> PID >> pStat >> E >> Px >> Py >> Pz ;
		    myfile >> PID >> pStat >> E >> Px >> Py >> Pz >> Eta >> Phi;
		    cout<<" "<<PID<<" "<< pStat<< " " << E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi<<endl;                 
		    double PT = TMath::Sqrt( pow(Px ,2.0) + pow(Py , 2.0) );double ModP = sqrt(Px*Px + Py*Py + Pz*Pz);
		    double Phi2 = atan2(Py,Px);double Eta = 0.5*log((ModP+Pz)/(ModP-Pz));
		    SumE= SumE + E; SumPx = SumPx + Px; SumPy= SumPy + Py; SumPz = SumPz + Pz;		    
		    if(pStat>=0  && fabs(Eta) < DetectorEtaCut && PT>0.01  &&  PID!=12 && PID!=14 && PID!=16 && PID!=18 )
		      {
			//cout<<PID<<" "<< pStat<< " " << E <<" "<< Px<<" " << Py<<" " << Pz<<endl;
			fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
			chargeList.push_back( pythia.particleData.charge( PID ) );
		      }
		    // Add this particle into SingleHadron spectrum
		    //if(pStat>=0 && pStat==11 && fabs(Eta) < SingleHadronEtaCut && PT>0.01  && fabs(PID) > 100 &&  pythia.particleData.charge( PID )!=0   )
		      if(pStat>=0  && fabs(Eta) < SingleHadronEtaCut && PT>0.01  && fabs(PID)==5)
		      {
			//cout<<"PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
			HistTempSingleHadron->Fill(PT); // fill pT distribution of single hadron
		      }
		    // Two-particle correlation function
                    if( (pStat>=0 ) && pTTrigMin< PT && PT < pTTrigMax &&  fabs(Eta) < EtaTrackMax )
                      {PhiVector.push_back(Phi2); EtaVector.push_back(Eta); pTVector.push_back(PT);
                      }
		  }		
		if(MyString.compare(EventLabel)==0)
		  {
		    //for(int i=0;i<6;i++){myfile >> EventLabel;}// cout<<EventLabel<<" ";}
		    for(int i=0;i<1;i++){myfile >> EventLabel;cout<<EventLabel<<" ";}
		    Events[k]++; EventsPerFile++;
		    cout<<HadronFile<<endl;
		  }
		
	      }
	    EventA<<"File i="<<FileIndex<<" has total events="<<EventsPerFile<<", Sum event with pT>110GeV are = "<<FakeEventCount<<endl;
	    cout<<"File i="<<FileIndex<<" has total events="<<EventsPerFile<<", Sum event with pT>110GeV are = "<<FakeEventCount<<endl;
	    //read hard cross section ptHat	    
	    myfile2>>HardCrossSection[k]>>HardCrossSectionError[k];
	    // end of reading cross section
	    myfile.close();
	    myfile2.close();
	    Events[k] = Events[k] -1;
	  } // end of reading one file		
	
	//For jet spectrum, weighted by cross section and combined through multiple pTHatBins
	for(int j=0;j<NpTJetBin;j++)
	  {
	    dNdpTCountJet[k][j]= HistTempJet->GetBinContent(j+1);
	    if(dNdpTCountJet[k][j] > 0.0)
	      {
		pTHardBinJetBinError[k][j] = (dNdpTCountJet[k][j]*HardCrossSection[k]/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut))*TMath::Sqrt( (1/dNdpTCountJet[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
		cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTempJet->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountJet[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut)<<endl;
	      }
	    else
	      {
		pTHardBinJetBinError[k][j] = 0.0;
		cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTempJet->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
	      }
	  }
        
	    //For single Hadron spectrum
	    for(int j=0;j<NpTSingleHadronBin;j++)
	      {
		dNdpTCountSingleHadron[k][j]= HistTempSingleHadron->GetBinContent(j+1);
		if(dNdpTCountSingleHadron[k][j]> 0.0)
		  {
		    pTHardBinSingleHadronBinError[k][j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection[k]/(Events[k]*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut))*TMath::Sqrt( (1/dNdpTCountSingleHadron[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
		    cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountSingleHadron[k][j]*HardCrossSection[k])/(Events[k]*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut)<<endl;
		  }
		else
		  {
		    pTHardBinSingleHadronBinError[k][j] = 0.0;
		    cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
		  }
	      }    
    
	    //Write histogram into a root file
	    HistTempJet->Write();
	    HistTempJetShape->Write();
	    HistTempSingleHadron->Write();
	    HistTempJetFF->Write();
	    HistTempJetMass->Write();
	    HistTemp2PCSignal->Write();
	    
	    TVector EventInfo(3);
	    EventInfo[0] = HardCrossSection[k];
	    EventInfo[1] = HardCrossSectionError[k];
	    EventInfo[2] = Events[k];
	    EventInfo.Write("EventInfo");
	    
	    TVector TriggeredJetInfo(2);
	    TriggeredJetInfo[0] = JetpTCut;
	    TriggeredJetInfo[1] = TriggeredJetNumber;
	    TriggeredJetInfo.Write("TriggeredJetInfo");
	    /*
	    //Plots for jet spectrum
	    double DifferentialJetCrossSection[NpTJetBin],DifferentialJetCrossSectionError[NpTJetBin],JetpT[NpTJetBin],JetpTError[NpTJetBin];
	    TGraphErrors * GEJet;	    
        cout<<"For ptHardBin = "<<k+1<<"\t CrossSection is Below "<<endl;
        for(int j=0; j<NpTJetBin;j++)  // fill in x, y observables together with error
	  {
            DifferentialJetCrossSection[j] = (dNdpTCountJet[k][j]*HardCrossSection)/(Events*((JetpTBin[j]+JetpTBin[j+1])/2.0)*2.0*JetEtaCut);
            DifferentialJetCrossSectionError[j] = pTHardBinJetBinError[k][j];
            JetpT[j] = (JetpTBin[j] + JetpTBin[j+1])/2.0;
            JetpTError[j] = (JetpTBin[j+1] - JetpTBin[j])/2.0;
            cout<<JetpT[j]<<"\t"<<DifferentialJetCrossSection[j]<<"\t"<<DifferentialJetCrossSectionError[j]<<endl;
	  }        
        // put observables and error bars into a root graph
	GEJet = new TGraphErrors(NpTJetBin,JetpT,DifferentialJetCrossSection,JetpTError,DifferentialJetCrossSectionError);
        char MyGraphName[100];
        sprintf(MyGraphName,"DifferentialJetCrossSectionBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJet->SetName(MyGraphName);
        GEJet->Write();
	    
        // For charged Hadron spectrum
        double DifferentialSingleHadronYield[NpTSingleHadronBin],DifferentialSingleHadronYieldError[NpTSingleHadronBin],SingleHadronpT[NpTSingleHadronBin],SingleHadronpTError[NpTSingleHadronBin];
        TGraphErrors * GESingleHadron;
        
        cout<<"For ptHardBin = "<<k+1<<"\t SingleHadron differential yield is Below "<<endl;
        for(int j=0; j<NpTSingleHadronBin;j++)
	  {
            DifferentialSingleHadronYield[j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection)/(Events*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*2.0*SingleHadronEtaCut);
            DifferentialSingleHadronYieldError[j] = pTHardBinSingleHadronBinError[k][j];
            SingleHadronpT[j] = (SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0;
            SingleHadronpTError[j] = (SingleHadronpTBin[j+1]-SingleHadronpTBin[j])/2.0;
            cout<<SingleHadronpT[j]<<"\t"<<DifferentialSingleHadronYield[j]<<"\t"<<DifferentialSingleHadronYieldError[j]<<endl;
        }
        
        GESingleHadron = new TGraphErrors(NpTSingleHadronBin,SingleHadronpT,DifferentialSingleHadronYield,SingleHadronpTError,DifferentialSingleHadronYieldError);
        char MyGraphName2[100];
        sprintf(MyGraphName2,"DifferentialSingleHadronYieldBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GESingleHadron->SetName(MyGraphName2);
        GESingleHadron->Write();
	   
        // for jet shape
        HistTempJetShape->Scale( (1.0/TriggeredJetNumber), "width" ); // average over number of triggered jets
        double ErrorNormJetShape;
        double NormJetShape = HistTempJetShape->IntegralAndError( 1, NrJetShapeBin,
                                                                  ErrorNormJetShape,
                                                                  "width" );
        TH1D *Norm = (TH1D*)HistTempJetShape->Clone("Normalization");
        Norm->Sumw2();Norm->SetTitle("Normalization");
        int nbins = Norm->GetSize();
        for( int i=1; i < nbins-1; i++){
            Norm->SetBinContent(i, NormJetShape);
            Norm->SetBinError(i, ErrorNormJetShape);
        }
        HistTempJetShape->Divide(Norm); // jet shape is defined as a self-normalized observable
        delete Norm;

        double JetShape[NrJetShapeBin],JetShapeError[NrJetShapeBin],JetShapeR[NrJetShapeBin],JetShapeRError[NrJetShapeBin];
        TGraphErrors * GEJetShape;
        
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Shape is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NrJetShapeBin;j++)
        {
            JetShape[j] = HistTempJetShape->GetBinContent(j+1);
            JetShapeError[j] = HistTempJetShape->GetBinError(j+1);
            JetShapeR[j] = (JetShaperBin[j]+JetShaperBin[j+1])/2.0;
            JetShapeRError[j] = (JetShaperBin[j+1]-JetShaperBin[j])/2.0;
            cout<<JetShapeR[j]<<"\t"<<JetShape[j]<<"\t"<<JetShapeError[j]<<endl;
        }
        GEJetShape
        = new TGraphErrors(NrJetShapeBin,JetShapeR,JetShape,JetShapeRError,JetShapeError);
        char MyGraphName3[100];
        sprintf(MyGraphName3,"JetShapeBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJetShape->SetName(MyGraphName3);
        GEJetShape->Write();
        
        
        
        // for jet fragmentation function
        HistTempJetFF->Scale( (1.0/TriggeredJetNumber), "width" );
        double JetFF[NzJetFFBin], JetFFError[NzJetFFBin], JetZ[NzJetFFBin], JetZError[NzJetFFBin];
        TGraphErrors * GEJetFF;
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Fragmentation Function is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NzJetFFBin;j++){
            JetFF[j] = HistTempJetFF->GetBinContent(j+1);
            JetFFError[j] = HistTempJetFF->GetBinError(j+1);
            JetZ[j] = (JetFFzBin[j]+JetFFzBin[j+1])/2.0;
            JetZError[j] = (JetFFzBin[j+1]-JetFFzBin[j])/2.0;
            cout<<JetZ[j]<<"\t"<<JetFF[j]<<"\t"<<JetFFError[j]<<endl;
        }
        GEJetFF
        = new TGraphErrors(NzJetFFBin,JetZ,JetFF,JetZError,JetFFError);
        char MyGraphName4[100];
        sprintf(MyGraphName4,"JetFFBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJetFF->SetName(MyGraphName4);
        GEJetFF->Write();
        
        // for jet mass
        HistTempJetMass->Scale( (1.0/TriggeredJetNumber), "width" );
        double JetMassDist[NJetMassBin], JetMassDistError[NJetMassBin], JetMass[NJetMassBin], JetMassError[NJetMassBin];
        TGraphErrors * GEJetMass;
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Mass Distribution is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NJetMassBin;j++){
            JetMassDist[j] = HistTempJetMass->GetBinContent(j+1);
            JetMassDistError[j] = HistTempJetMass->GetBinError(j+1);
            JetMass[j] = (JetMassBin[j]+JetMassBin[j+1])/2.0;
            JetMassError[j] = (JetMassBin[j+1]-JetMassBin[j])/2.0;
            cout<<JetMass[j]<<"\t"<<JetMassDist[j]<<"\t"<<JetMassDistError[j]<<endl;
        }
        GEJetMass
        = new TGraphErrors(NJetMassBin,JetMass,JetMassDist,JetMassError,JetMassDistError);
        char MyGraphName5[100];
        sprintf(MyGraphName5,"JetMassDistBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJetMass->SetName(MyGraphName5);
        GEJetMass->Write();

        delete GEJet; delete GESingleHadron; delete GEJetShape; delete GEJetFF; delete GEJetMass;
	    */
    //delete outFile;
	EventA<<"(AvgE,AvgPx,AvgPy,AvgPz) = ("<<SumE/Events[k]<<","<< SumPx/Events[k]<<","<<SumPy/Events[k]<<","<<SumPz/Events[k]<<")"<<endl;
	EventR<<"For bin k="<<k<<", total events = "<<Events[k]<<endl; 
    } //k-loop ends here (pTHatBin loop)
    
    delete HistTempJet; delete HistTempSingleHadron;
    delete HistTempJetShape; delete HistTempJetFF; delete HistTempJetMass; 
    delete HistTemp2PCSignal;
    
    //Full spectrum Below
    //pTHat cross section vs pTHatBin
    double pTHatError[NpTHardBin], pTHatAvg[NpTHardBin];
    for(int k=0; k<NpTHardBin; k++)
      {
        pTHatError[k]=0.0;
        pTHatAvg[k]= ((pTHatMin[k] + pTHatMax[k])/2.0);
      }
    TGraphErrors *GHardCrossSection = new TGraphErrors(NpTHardBin,pTHatAvg,HardCrossSection,pTHatError,HardCrossSectionError);
    GHardCrossSection->SetName("HardCrossSection");
    GHardCrossSection->Write();
    //Plots for full jet spectrum
    double DifferentialJetCrossSection[NpTJetBin],DifferentialJetCrossSectionError[NpTJetBin],JetpT[NpTJetBin],JetpTError[NpTJetBin];
    TGraphErrors * GEJet;
    for(int j=0; j<NpTJetBin; j++)
      {
	TotalDifferentialJetCrossSection[j] = 0.0;
	TotalDifferentialJetCrossSectionError[j] = 0.0;
      }

    for(int k=0;k<NpTHardBin;k++)
      {cout<<"For ptHardBin = "<<k+1<<"\t CrossSection is Below "<<endl;
        for(int j=0; j<NpTJetBin;j++)
	  {
            DifferentialJetCrossSection[j] = (dNdpTCountJet[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut);
            DifferentialJetCrossSectionError[j] = pTHardBinJetBinError[k][j];
            TotalDifferentialJetCrossSection[j] = TotalDifferentialJetCrossSection[j] + DifferentialJetCrossSection[j];
            TotalDifferentialJetCrossSectionError[j] = TotalDifferentialJetCrossSectionError[j] + TMath::Power( pTHardBinJetBinError[k][j], 2.0);
            JetpT[j] = (JetpTBin[j]+JetpTBin[j+1])/2.0;
            JetpTError[j] = (JetpTBin[j+1]-JetpTBin[j])/2.0;
            cout<<JetpT[j]<<"\t"<<DifferentialJetCrossSection[j]<<"\t"<<DifferentialJetCrossSectionError[j]<<endl;
	  }
        GEJet = new TGraphErrors(NpTJetBin,JetpT,DifferentialJetCrossSection,JetpTError,DifferentialJetCrossSectionError);
        char MyGraphName[100];
        sprintf(MyGraphName,"DifferentialJetCrossSectionBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJet->SetName(MyGraphName);
        GEJet->Write();
      }
    cout<<"Final results for total crossSection"<<endl;
    cout<<"JetpT \t"<<"Total Differential CrossSection \t"<<"Error"<<endl;
    for(int j=0;j<NpTJetBin;j++)
      {
	TotalDifferentialJetCrossSectionError[j] =TMath::Sqrt(TotalDifferentialJetCrossSectionError[j]);
        cout<<JetpT[j]<<"\t"<<TotalDifferentialJetCrossSection[j]<<" +/- \t"<<TotalDifferentialJetCrossSectionError[j]<<endl;
      }
    cout<<"Create Final TGraphError plot "<<endl;
    //Final Plot and save also
    TGraphErrors *GEJetTotal = new TGraphErrors(NpTJetBin,JetpT,TotalDifferentialJetCrossSection,JetpTError,TotalDifferentialJetCrossSectionError);
    GEJetTotal->SetName("TotalDifferentialJetCrossSection");
    GEJetTotal->Write();

    //For single hadron spectrum
    double DifferentialSingleHadronYield[NpTSingleHadronBin],DifferentialSingleHadronYieldError[NpTSingleHadronBin],SingleHadronpT[NpTSingleHadronBin],SingleHadronpTError[NpTSingleHadronBin];
    TGraphErrors * GESingleHadron;
    for(int j=0; j<NpTSingleHadronBin;j++)
      {
        TotalDifferentialSingleHadronYield[j] = 0.0;
        TotalDifferentialSingleHadronYieldError[j] = 0.0;
      }

    for(int k=0;k<NpTHardBin;k++)
      {cout<<"For ptHardBin = "<<k+1<<"\t SingleHadron differential yield is Below "<<endl;
        for(int j=0; j<NpTSingleHadronBin;j++)
	  {
            DifferentialSingleHadronYield[j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection[k])/(Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*68*2*M_PI*((SingleHadronpTBin[j+1]+SingleHadronpTBin[j])/2.0)*2.0*SingleHadronEtaCut);
            DifferentialSingleHadronYieldError[j] = pTHardBinSingleHadronBinError[k][j];
            TotalDifferentialSingleHadronYield[j] = TotalDifferentialSingleHadronYield[j] + DifferentialSingleHadronYield[j];
            TotalDifferentialSingleHadronYieldError[j] = TotalDifferentialSingleHadronYieldError[j] + TMath::Power( pTHardBinSingleHadronBinError[k][j], 2.0);
            SingleHadronpT[j] = (SingleHadronpTBin[j+1] + SingleHadronpTBin[j])/2.0;
            SingleHadronpTError[j] = (SingleHadronpTBin[j+1] - SingleHadronpTBin[j])/2.0;
            cout<<SingleHadronpT[j]<<"\t"<<DifferentialSingleHadronYield[j]<<"\t"<<DifferentialSingleHadronYieldError[j]<<endl;
	  }

        GESingleHadron = new TGraphErrors(NpTSingleHadronBin,SingleHadronpT,DifferentialSingleHadronYield,SingleHadronpTError,DifferentialSingleHadronYieldError);
        char MyGraphName2[100];
        sprintf(MyGraphName2,"DifferentialSingleHadronYieldBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GESingleHadron->SetName(MyGraphName2);
        GESingleHadron->Write();
      }
    cout<<"Final results for total differential single hadron yield"<<endl;
    cout<<"SingleHadronpT \t"<<"total differential single hadron yield \t"<<"Error"<<endl;
    for(int j=0;j<NpTSingleHadronBin;j++)
      {
        TotalDifferentialSingleHadronYieldError[j] =TMath::Sqrt(TotalDifferentialSingleHadronYieldError[j]);
        cout<<SingleHadronpT[j]<<"\t"<<TotalDifferentialSingleHadronYield[j]<<" +/- \t"<<TotalDifferentialSingleHadronYieldError[j]<<endl;
      }
    cout<<"Create Final TGraphError plot "<<endl;
    // Final Plot and save also
    TGraphErrors *GESingleHadronTotal  = new TGraphErrors(NpTSingleHadronBin, SingleHadronpT, TotalDifferentialSingleHadronYield, SingleHadronpTError, TotalDifferentialSingleHadronYieldError);
    GESingleHadronTotal->SetName("TotalDifferentialSingleHadronYield");
    GESingleHadronTotal->Write();

    delete GEJet; delete GESingleHadron;
    delete GEJetTotal; delete GESingleHadronTotal;
    delete GHardCrossSection;
    
    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    EventR<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    EventA.close();
    EventR.close();
    return 0;
}
