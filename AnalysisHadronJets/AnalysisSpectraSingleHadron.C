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
#include "TROOT.h"
using namespace std;
using namespace Pythia8;
using namespace fjcore;
int main(int argc, char **argv){
  double PI=3.1415926;
    int StartTime = time(NULL);    
    //gROOT->SetBatch(); 
    //TApplication theApp("hist", &argc, argv);   
    //By default we want to look at only one pTHatBin spectrum
    //int pTHatMin[1]={11};
    //int pTHatMax[1]={13};
    //int NpTHardBin = sizeof(pTHatMin)/sizeof(pTHatMin[0]); //sizeof(pTHatMin[]); or sizeof(pTHatMax[]); // # of bin we analyze
    int *pTHatMin;
    int *pTHatMax;
    int NpTHardBin;
    //[1]InputDataFileDIR, [2]OutputDataFileName,  [3]=CMS or ATLAS or ALICE or STAR, [4]= 200, 2760 or 5020, 7000, [5]SingleHadronEtaCut,          [6]EtaCutFlag, [7]=ChargedFlag=0 or 1 (all or charged), [8]=Jetscape or Pythia, [9] Hadron or Parton [10] BGS
    char InputDataFileDIR[100000], OutputDataFileName[100000], JetscapeORPythia[100], HadronORParton[100];
    sprintf(InputDataFileDIR,"%s",argv[1]);
    sprintf(OutputDataFileName,"%s",argv[2]);
    sprintf(JetscapeORPythia,"%s",argv[8]);
    sprintf(HadronORParton,"%s",argv[9]);

    string ExpName=std::string(argv[3]);
    int Ecm = atoi(argv[4]);
    double SingleHadronEtaCut=atof(argv[5])/1.0;
    int EtaCutFlag =atoi(argv[6]); //0 or 1; 0 to use rapididty cut, 1 to use eta cut
    int ChargedFlag = atoi(argv[7]); // 0 or 1; 0 means all particle, 1 means only charged
    int BGS=atoi(argv[10]);
    
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

    if(  Ecm == 5020 )
      {
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


    int NumberOfFiles=1;
    int Events[NpTHardBin];
    double HardCrossSection[NpTHardBin], HardCrossSectionError[NpTHardBin];
    double *SingleHadronpTBin;
    int NpTSingleHadronBin;    
    //Variable for single hadron spectrum
    //double SingleHadronpTBin[26] = {0.5, 0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 12, 14, 16, 18, 20};
    //int NpTSingleHadronBin = 26-1;
    if(  Ecm == 5020 && ExpName=="CMS"  && SingleHadronEtaCut==1.0  && EtaCutFlag==1)
      {
	double Bin[38] ={0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6,  1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24,                        28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400};
	NpTSingleHadronBin = 38-1;
	SingleHadronpTBin = new double[NpTSingleHadronBin+1];
	for(int i=0; i<NpTSingleHadronBin+1; i++)
	  {
	    SingleHadronpTBin[i] = Bin[i];
	  }
      }

    double dNdpTCountSingleHadron[NpTHardBin][NpTSingleHadronBin];  //[ptHatBin] [Regular pt]
    double pTHardBinSingleHadronBinError[NpTHardBin][NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYield[NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYieldError[NpTSingleHadronBin];
    double InelasticCS = 70.0;//CMS 5.02TeV (70mb)


    // Histograms.1. Number of jets vs pT of the jet. 2. Number of charged hadron vs pT of the hadron
    TH1D *HistTempSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
    TH1D *HistTempSingleHadronPositive = new TH1D("SingleHadronSpectrumPositive", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin);
    TH1D *HistTempSingleHadronNegative = new TH1D("SingleHadronSpectrumNegative", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin);

    ofstream EventR, EventA, foutput; char VarFileName[100000];
    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/EventRecord%s.txt",OutputDataFileName);
    EventR.open(VarFileName,ios::out);
    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/CodeCheck%s.txt",OutputDataFileName);
    EventA.open(VarFileName,ios::out);
    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/%s.root",OutputDataFileName); // name of the output root file  
    TFile* outFile = new TFile( VarFileName, "RECREATE");

    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/FastJet/Files/%s.txt",OutputDataFileName);
    foutput.open(VarFileName,ios::out);
    cout<<"Final result will be store in  File="<<VarFileName<<endl;
    std::vector <int> chargeList;    
    Pythia8::Pythia pythia;//("",false);

    EventA<<"OutputFile Name is "<<VarFileName<<endl;
    cout<<"These are pTHat loops "<<endl;
   
    // For loop to open different pTHat bin files, in this default example, only one file
    for (int k = 0; k<NpTHardBin; ++k) //NpTHardBin
    {
        char pTBinString[100];        
        sprintf(pTBinString,"Current pTHatBin is %i (%i,%i) GeV",k,pTHatMin[k],pTHatMax[k]);
	EventA<<pTBinString<<endl;
        int  SN=0,PID=0,pStat=0;
        double Px, Py, Pz, E, Eta, Phi;
	int EventsPerFile=0;
        double SumE=0, SumPx=0, SumPy=0, SumPz=0;
	Events[k]=0;
        // Reset for each pTHardBin
        char HistName[10000];        
      
        HistTempSingleHadron->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadron->SetName(HistName); 
	
	HistTempSingleHadronPositive->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%iPositive",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadronPositive->SetName(HistName);

	HistTempSingleHadronNegative->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%iNegative",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadronNegative->SetName(HistName);

	for(int FileIndex=0; FileIndex<NumberOfFiles; FileIndex++)
	  {
	    char HadronFile[30000], HardCrossSectionFile[100000];
	    sprintf(HadronFile,"%s/%s%sListBin%i_%i.out",InputDataFileDIR,JetscapeORPythia,HadronORParton, pTHatMin[k],pTHatMax[k]);  // name of the input file  
	    sprintf(HardCrossSectionFile,"%s/SigmaHardBin%i_%i.out",InputDataFileDIR,pTHatMin[k],pTHatMax[k]);
	    sprintf(pTBinString,"Current Bin is %i pTHatBin(%i,%i) GeV and file is %i",k,pTHatMin[k],pTHatMax[k], FileIndex);        
	    EventA<<pTBinString<<endl;
	    EventA<<"Hadron File ="<<HadronFile<<endl;
	    EventA<<"HardCrossSectionFile="<<HardCrossSectionFile<<endl;
	    EventsPerFile =0;
	    // Read file
	    string EventLabel, MyString("#");
	    //string MyString("Event");
	    ifstream myfile;  myfile.open(HadronFile,ios::in);
	    ifstream myfile2; myfile2.open(HardCrossSectionFile,ios::in);
	    //myfile >>EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel;
	    while (myfile >>EventLabel )
	      {
		
		if(MyString.compare(EventLabel)!=0)
		  { //myfile >> PID >> pStat >> E >> Px >> Py >> Pz ;
		    myfile >> PID >> pStat >> E >> Px >> Py >> Pz >> Eta >> Phi;
		    cout<<" "<<PID<<" "<< pStat<< " " << E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi<<endl;                 
		    double PT = TMath::Sqrt( pow(Px ,2.0) + pow(Py , 2.0) );double ModP = sqrt(Px*Px + Py*Py + Pz*Pz);
		    double Phi2 = atan2(Py,Px);double Eta = 0.5*log((ModP+Pz)/(ModP-Pz));
		    double Rapidity= 0.5*log((E+Pz)/(E-Pz));
		    
		    if(EtaCutFlag==0){Eta = Rapidity;}
		    // Add this particle into SingleHadron spectrum
		    if(pStat>=0 && fabs(Eta) < SingleHadronEtaCut && PT>0.01)
		      {
			//cout<<"PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
			if(ChargedFlag==1 && pythia.particleData.charge( PID )!=0 && fabs(PID) > 100)
			  {
			    HistTempSingleHadronPositive->Fill(PT); // fill pT distribution of single hadron
			  }

			if(ChargedFlag==0)
                          {
                            HistTempSingleHadronPositive->Fill(PT); // fill pT distribution of single hadron  
			  }
			SumE= SumE + E; SumPx = SumPx + Px; SumPy= SumPy + Py; SumPz = SumPz + Pz;
		      }

		    if(pStat<0 && fabs(Eta) < SingleHadronEtaCut) 
		    {
		      //cout<<"PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
		      HistTempSingleHadronNegative->Fill(PT); // fill pT distribution of single hadron       
		      SumE= SumE - E; SumPx = SumPx - Px; SumPy= SumPy - Py; SumPz = SumPz - Pz;
		    }
		    
		  }
		
		if(MyString.compare(EventLabel)==0)
		  {
		    for(int i=0;i<8;i++) { myfile >> EventLabel;cout<<EventLabel<<" ";}
		    //for(int i=0;i<1;i++){myfile >> EventLabel;cout<<EventLabel<<" ";}
		    Events[k]++; EventsPerFile++;
		    cout<<"\n"<<HadronFile<<endl;
		  }
	      }
	    EventA<<"File k="<<k<<",i="<<FileIndex<<" has total events="<<EventsPerFile<<endl;
	    cout<<"File k="<<k<<",i="<<FileIndex<<" has total events="<<EventsPerFile<<endl;
	    //read hard cross section ptHat	    
	    myfile2>>HardCrossSection[k]>>HardCrossSectionError[k];
	    // end of reading cross section
	    myfile.close();
	    myfile2.close();
	  }//end of reading  files for same pTHatBins 
	if(BGS==1)
	  {
	    HistTempSingleHadron->Add(HistTempSingleHadronPositive, HistTempSingleHadronNegative,1.0,-1.0);
	  }
	else
	  {
	    HistTempSingleHadron->Add(HistTempSingleHadronPositive, HistTempSingleHadronNegative,1.0,0.0);
	  }
	    //For single Hadron spectrum
	    for(int j=0;j<NpTSingleHadronBin;j++)
	      {
		dNdpTCountSingleHadron[k][j]= HistTempSingleHadron->GetBinContent(j+1);
		if(dNdpTCountSingleHadron[k][j]> 0.0)
		  {
		    //Yield
		    //pTHardBinSingleHadronBinError[k][j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection[k]/(Events[k]*InelasticCS*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut))*TMath::Sqrt( (1/dNdpTCountSingleHadron[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
		    //cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountSingleHadron[k][j]*HardCrossSection[k])/(Events[k]*InelasticCS*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut)<<endl;
		    //cross section		    
		    pTHardBinSingleHadronBinError[k][j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection[k]/(Events[k]*InelasticCS*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut))*TMath::Sqrt( (1/dNdpTCountSingleHadron[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
                    cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountSingleHadron[k][j]*HardCrossSection[k])/(Events[k]*InelasticCS*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut)<<endl;
               
		  }
		else
		  {
		    pTHardBinSingleHadronBinError[k][j] = 0.0;
		    cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
		  }
	      }    
    
	    //Write histogram into a root file
	    if(BGS==1)
	      {
		HistTempSingleHadron->Write();
		HistTempSingleHadronPositive->Write();
		HistTempSingleHadronNegative->Write();
	      }
	    else
	      {
		HistTempSingleHadron->Write();
	      }

	    TVector EventInfo(3);
	    EventInfo[0] = HardCrossSection[k];
	    EventInfo[1] = HardCrossSectionError[k];
	    EventInfo[2] = Events[k];
	    EventInfo.Write("EventInfo");	    
 
    //delete outFile;
	EventA<<"(AvgE,AvgPx,AvgPy,AvgPz) = ("<<SumE/Events[k]<<","<< SumPx/Events[k]<<","<<SumPy/Events[k]<<","<<SumPz/Events[k]<<")"<<endl;
	EventR<<"For bin k="<<k<<",pTHat("<< pTHatMin[k]<<","<<pTHatMax[k] <<")GeV"<<"\t total events = "<<Events[k]<<endl; 
    } //k-loop ends here (pTHatBin loop)
    
    delete HistTempSingleHadron;
    delete HistTempSingleHadronPositive, HistTempSingleHadronNegative;
       
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
            DifferentialSingleHadronYield[j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection[k])/(Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*InelasticCS*2*M_PI*((SingleHadronpTBin[j+1]+SingleHadronpTBin[j])/2.0)*2.0*SingleHadronEtaCut);
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
    foutput<<"#HadronpT \t"<<"DifferentialYieldORCrossSection \t"<<"JetpTBinWidth/2.0 \t"<<"DCSorYieldError \t"<<"DCSorYieldErr%"<<endl;
    for(int j=0;j<NpTSingleHadronBin;j++)
      {
        TotalDifferentialSingleHadronYieldError[j] =TMath::Sqrt(TotalDifferentialSingleHadronYieldError[j]);
        cout<<SingleHadronpT[j]<<"\t"<<TotalDifferentialSingleHadronYield[j]<<" +/- \t"<<TotalDifferentialSingleHadronYieldError[j]<<endl;
	foutput<<SingleHadronpT[j]<<"\t"<<TotalDifferentialSingleHadronYield[j]<<"\t"<<SingleHadronpTError[j]<<"\t"<<TotalDifferentialSingleHadronYieldError[j]<<"\t"<<TotalDifferentialSingleHadronYieldError[j]*100/TotalDifferentialSingleHadronYield[j]<<endl;
      }
    cout<<"Create Final TGraphError plot "<<endl;
    // Final Plot and save also
    TGraphErrors *GESingleHadronTotal  = new TGraphErrors(NpTSingleHadronBin, SingleHadronpT, TotalDifferentialSingleHadronYield, SingleHadronpTError, TotalDifferentialSingleHadronYieldError);
    GESingleHadronTotal->SetName("TotalDifferentialSingleHadronYield");
    GESingleHadronTotal->Write();
    foutput.close();
    
    delete GESingleHadron;
    delete GESingleHadronTotal;
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
    outFile->Close();
        
    return 0;
}
