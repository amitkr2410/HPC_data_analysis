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
  
    //[1]InputDataFileDIR, [2]OutputDataFileName,  [3]=CMS or ATLAS or ALICE or STAR, [4]= 200, 2760 or 5020, 7000, [5]SingleHadronEtaCut,          [6]EtaCutFlag, [7]=ChargedFlag=0 or 1 (all or charged), [8]=Jetscape or Pythia, [9] Hadron or Parton [10] BGS  [11] Centrality =0-10  [12] SoftV2Flag=0 or 1 (off, on)
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
    string Centrality=std::string(argv[11]);
    int SoftV2Flag = atoi(argv[12]);  //0= without soft v2, 1=with soft v2
    cout<<"Working 1"<<endl;

    //By default we want to look at only one pTHatBin spectrum
    //int pTHatMin[1]={500};
    //int pTHatMax[1]={550};
    //int NpTHardBin = sizeof(pTHatMin)/sizeof(pTHatMin[0]); //sizeof(pTHatMin[]); or sizeof(pTHatMax[]); // # of bin we analyze                                                                                            
    
      int *pTHatMin;
      int *pTHatMax;
      int NpTHardBin;
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
	int Bin1[54] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120,                                     	130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400,                                            450, 500,   550, 600, 700, 800, 900, 1000 };
	int Bin2[54] = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130,                  			140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450,                                            500, 550, 600, 700, 800, 900, 1000, -1 };
	NpTHardBin = 54;
	pTHatMin = new int[NpTHardBin];
	pTHatMax = new int[NpTHardBin];
	for(int i=0; i<NpTHardBin; i++)
	  {
	    pTHatMin[i] = Bin1[i];       pTHatMax[i] = Bin2[i];
	  }
      }
    
    if(  Ecm == 5020 )
      {                 //{ 1, 2, 3, 4, 5}
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
    if(  Ecm == 2760 && ExpName=="CMS"  && SingleHadronEtaCut==1.0  && EtaCutFlag==1)
      {
        double Bin[16] ={5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24,                              			28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6};
        NpTSingleHadronBin = 16-1;
        SingleHadronpTBin = new double[NpTSingleHadronBin+1];
        for(int i=0; i<NpTSingleHadronBin+1; i++)
          {
            SingleHadronpTBin[i] = Bin[i];
          }
      }


    if(  Ecm == 5020 && ExpName=="CMS"  && SingleHadronEtaCut==1.0  && EtaCutFlag==1 && BGS>=0)
      {
        double Bin[21] ={5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24,                        28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400};
        NpTSingleHadronBin = 21-1;
	SingleHadronpTBin = new double[NpTSingleHadronBin+1];
        for(int i=0; i<NpTSingleHadronBin+1; i++)
          {
            SingleHadronpTBin[i] = Bin[i];
          }
      }

    if(  Ecm == 5020 && ExpName=="Atlas"  && SingleHadronEtaCut==2.5  && EtaCutFlag==1 && BGS>=0)
      {
	double Bin[15]={5, 5.5, 6, 6.5, 7.5, 9.5, 12, 16, 20, 25, 30, 35, 40, 50, 60};
        NpTSingleHadronBin = 15-1;
        SingleHadronpTBin = new double[NpTSingleHadronBin+1];
        for(int i=0; i<NpTSingleHadronBin+1; i++)
          {
            SingleHadronpTBin[i] = Bin[i];
          }
      }

    double dNdpTCountSingleHadronPositive[NpTHardBin][NpTSingleHadronBin], dNdpTCountSingleHadronNegative[NpTHardBin][NpTSingleHadronBin];
    double dNCosinedpTSingleHadronPositive[NpTHardBin][NpTSingleHadronBin];
    double dNCosinedpTSingleHadronNegative[NpTHardBin][NpTSingleHadronBin];
    double SumSquareCosineWeight[NpTHardBin][NpTSingleHadronBin];
    double SumSquareSoftV2[NpTHardBin];    
    double pTHardBinSingleHadronBinError[NpTHardBin][NpTSingleHadronBin];
    double V2Plus[NpTSingleHadronBin], V2PlusError[NpTSingleHadronBin];
    double V2PlusMinus[NpTSingleHadronBin], V2PlusMinusError[NpTSingleHadronBin];
    double V2Minus[NpTSingleHadronBin], V2MinusError[NpTSingleHadronBin];
    double InelasticCS =1.0;//CMS 5.02TeV (70mb)


    // Histograms.1. Number of jets vs pT of the jet. 2. Number of charged hadron vs pT of the hadron

    TH1D *HistTempWeightedCosinePositive = new TH1D("WeightedCosineSpectrumPositive", "WeightedCosine plus Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin);
    TH1D *HistTempWeightedCosineNegative = new TH1D("WeightedCosineSpectrumNegative", "WeightedCosine minus Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin);

    TH1D *HistTempSingleHadronPositive = new TH1D("SingleHadronSpectrumPositive", "SingleHadron plus Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin);
    TH1D *HistTempSingleHadronNegative = new TH1D("SingleHadronSpectrumNegative", "SingleHadron minus Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin);
    TH1D *HistTemp = new TH1D("EventByEventWeightFactor","Weight Factor single event",NpTSingleHadronBin, SingleHadronpTBin);

    ofstream EventR, EventA, foutputPlus, foutputPlusMinus; char VarFileName[100000];
    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/FilesV2/EventRecord%s.txt",OutputDataFileName);
    EventR.open(VarFileName,ios::out);
    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/FilesV2/CodeCheck%s.txt",OutputDataFileName);
    EventA.open(VarFileName,ios::out);
    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/FilesV2/%s.root",OutputDataFileName); // name of the output root file  
    TFile* outFile = new TFile( VarFileName, "RECREATE");

    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/FilesV2/%sPlus.txt",OutputDataFileName);
    foutputPlus.open(VarFileName,ios::out);
    sprintf(VarFileName,"/wsu/home/fy/fy41/fy4125/RUN/Analysis/FilesV2/%sPlusMinus.txt",OutputDataFileName);
    foutputPlusMinus.open(VarFileName,ios::out);
    
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
        double Px, Py, Pz, E, Eta, Phi, EventPlaneAngle=0;
	double soft_v2, soft_v2_error, psi2, psi2_error, soft_v3, soft_v3_error, psi3, psi3_error;
	int EventsPerFile=0;
        double SumE=0, SumPx=0, SumPy=0, SumPz=0;
	Events[k]=0;
        // Reset for each pTHardBin
        char HistName[10000];        
	for(int j=0;j<NpTSingleHadronBin;j++)
	  {
	    SumSquareCosineWeight[k][j] = 0;
	  }
	SumSquareSoftV2[k] =0;
	
	HistTempWeightedCosinePositive->Reset();
        sprintf(HistName,"CountWeightedCosineVspTSingleHadronSpectrumBin%i_%iPositive",pTHatMin[k],pTHatMax[k]);
        HistTempWeightedCosinePositive->SetName(HistName);

	HistTempWeightedCosineNegative->Reset();
        sprintf(HistName,"CountWeightedCosineVspTSingleHadronSpectrumBin%i_%iNegative",pTHatMin[k],pTHatMax[k]);
        HistTempWeightedCosineNegative->SetName(HistName);	

	HistTempSingleHadronPositive->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%iPositive",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadronPositive->SetName(HistName);

	HistTempSingleHadronNegative->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%iNegative",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadronNegative->SetName(HistName);

	HistTemp->Reset();

	for(int FileIndex=0; FileIndex<NumberOfFiles; FileIndex++)
	  {
	    char HadronFile[30000], HardCrossSectionFile[100000], SoftV2Coefficient[10000];
	    sprintf(HadronFile,"%s/%s%sListBin%i_%i.out",InputDataFileDIR,JetscapeORPythia,HadronORParton, pTHatMin[k],pTHatMax[k]);  // name of the input file  
	    sprintf(HardCrossSectionFile,"/wsu/home/groups/maj-shen/AAPaperData/pTHatCrossSection_%iGeV/SigmaHardBin%i_%i.out",  Ecm, pTHatMin[k],pTHatMax[k]);
	    //sprintf(HardCrossSectionFile,"%s/SigmaHardBin%i_%i.out",InputDataFileDIR,pTHatMin[k],pTHatMax[k]);
	    sprintf(pTBinString,"Current Bin is %i pTHatBin(%i,%i) GeV and file is %i",k,pTHatMin[k],pTHatMax[k], FileIndex);        
	    EventA<<pTBinString<<endl;
	    EventA<<"Hadron File ="<<HadronFile<<endl;
	    EventA<<"HardCrossSectionFile="<<HardCrossSectionFile<<endl;
	    EventsPerFile =0;
	    // Read file
	    string EventLabel, MyString("#"), SoftV2Label, DummyS;
	    //string MyString("Event");
	    ifstream myfile;  myfile.open(HadronFile,ios::in);
	    ifstream myfile2; myfile2.open(HardCrossSectionFile,ios::in);
	    ifstream myfile3; 
	    if(SoftV2Flag==1)
	      { 
		sprintf(SoftV2Coefficient,"/wsu/home/groups/maj-shen/AAPaperData/MATTER_LBT_RunningAlphaS_Q2qhat/SoftAnisotropy_Coefficients_5TeV_%s.txt",Centrality.c_str());
		myfile3.open(SoftV2Coefficient,ios::in);
		myfile3>> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label >> SoftV2Label;
		EventA<<"Soft v2 File ="<<SoftV2Coefficient<<endl;
	      }
	    //myfile >>EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel>> EventLabel;
	    while (myfile >>EventLabel  )
	      {
		
		if(MyString.compare(EventLabel)!=0)
		  { //myfile >> PID >> pStat >> E >> Px >> Py >> Pz ;
		    myfile >> PID >> pStat >> E >> Px >> Py >> Pz >> Eta >> Phi;
		    double PT = TMath::Sqrt( pow(Px ,2.0) + pow(Py , 2.0) );double ModP = sqrt(Px*Px + Py*Py + Pz*Pz);
		    double Phi2 = atan2(Py,Px) ;double Eta = 0.5*log((ModP+Pz)/(ModP-Pz));
		    double Rapidity= 0.5*log((E+Pz)/(E-Pz));
		    double COSINE=0;
		    COSINE =  soft_v2*std::cos(2.0*(Phi2-EventPlaneAngle));

		    //cout<<" "<<PID<<" "<< pStat<< " " << E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi*180.0/PI<<" "<<Phi2*180.0/PI<<" "<<EventPlaneAngle*180.0/PI<<endl;
		    if(EtaCutFlag==0){Eta = Rapidity;}
		    // Add this particle into SingleHadron spectrum
		    if(pStat>=0 && fabs(Eta) < SingleHadronEtaCut && PT>0.01)
		      {
			//cout<<"PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
			if(ChargedFlag==1 && pythia.particleData.charge( PID )!=0 )
			  {
			    HistTempWeightedCosinePositive->Fill(PT, COSINE);
			    HistTempSingleHadronPositive->Fill(PT); // fill pT distribution of single hadron	
			    HistTemp->Fill(PT,COSINE);
			    
			  }
			
			if(ChargedFlag==0)
                          {HistTempWeightedCosinePositive->Fill(PT, COSINE);
                            HistTempSingleHadronPositive->Fill(PT); // fill pT distribution of single hadron  
			    HistTemp->Fill(PT,COSINE);
			    
			  }
			SumE= SumE + E; SumPx = SumPx + Px; SumPy= SumPy + Py; SumPz = SumPz + Pz;
		      }

		    if(pStat<0 && fabs(Eta) < SingleHadronEtaCut && pythia.particleData.charge( PID )!=0 ) 
		    {
		      //cout<<"PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
		      HistTempWeightedCosineNegative->Fill(PT, COSINE);
		      HistTempSingleHadronNegative->Fill(PT); // fill pT distribution of single hadron       
		      SumE= SumE - E; SumPx = SumPx - Px; SumPy= SumPy - Py; SumPz = SumPz - Pz;
		    }
		    
		  }
		
		if(MyString.compare(EventLabel)==0)
		  {
		    myfile >> EventPlaneAngle; cout<<EventPlaneAngle<<" ";
		    for(int i=0;i<7;i++) { myfile >> EventLabel;cout<<EventLabel<<" ";}
		    //for(int i=0;i<7;i++){myfile >> EventLabel;cout<<EventLabel<<" ";}
		    cout<<endl;
		    soft_v2 = 1.0;
		    if(SoftV2Flag==1)
		      {
			myfile3>>SoftV2Label >> soft_v2 >> DummyS >> psi2 >> DummyS >> DummyS >> DummyS >> DummyS >> DummyS;
			cout<<SoftV2Label <<" "<< soft_v2 <<" "<< DummyS <<" "<< psi2 <<" "<< DummyS <<" "<<DummyS <<" "<< DummyS <<" "<< DummyS <<" "<< DummyS <<" "<<Events[k]+1<<endl;
			if( int((EventPlaneAngle-psi2)*10000) !=0 ) 
			  {cout<<"\n Error soft coefficient are mismatched:("<<EventPlaneAngle<<","<<psi2<<"), pTHat= "<<pTHatMin[k]<<", Events="<<Events[k]+1<<endl; 
			    EventA<<" \n Error soft coefficent are mismatched:("<<EventPlaneAngle<<","<<psi2<<"), pTHat= "<<pTHatMin[k]<<", Events="<<Events[k]+1<<endl;//cin>>DummyS; 
			    soft_v2 = 1.0;
			  }
		      }
		    for(int j=0;j<NpTSingleHadronBin;j++)
		      {			
			SumSquareCosineWeight[k][j] = SumSquareCosineWeight[k][j] + pow(HistTemp->GetBinContent(j+1),2);
			cout<<"For Event="<<Events[k]<<",j="<<j<<", one event BinContent is "<<HistTemp->GetBinContent(j+1)<<endl;
		      }
		    SumSquareSoftV2[k] = SumSquareSoftV2[k] + pow(soft_v2, 2);
		    HistTemp->Reset();
		    Events[k]++; EventsPerFile++;
		    cout<<"\n"<<HadronFile<<endl;
		  }
	      }
	    EventA<<"File k="<<k<<",i="<<FileIndex<<" has total events="<<EventsPerFile<<endl;
	    cout<<"File k="<<k<<",i="<<FileIndex<<" has total events="<<EventsPerFile<<endl;
	    //read hard cross section ptHat	    
	    myfile2>>HardCrossSection[k]>>HardCrossSectionError[k];	    
	    EventA<<"HardCross section = "<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl;
	    cout<<"HardCross section = "<<HardCrossSection[k]<<"+/-"<<HardCrossSectionError[k]<<endl;
	    // end of reading cross section
	    for(int j=0;j<NpTSingleHadronBin;j++)
	      {
		SumSquareCosineWeight[k][j] = SumSquareCosineWeight[k][j] + pow(HistTemp->GetBinContent(j+1),2);
	      }
	    
	    myfile.close();
	    myfile2.close();
	    myfile3.close();
	  }//end of reading  files for same pTHatBins 
	
	    //For single Hadron spectrum
	    for(int j=0;j<NpTSingleHadronBin;j++)
	      {
		dNdpTCountSingleHadronPositive[k][j] = HistTempSingleHadronPositive->GetBinContent(j+1);
		dNdpTCountSingleHadronNegative[k][j] = HistTempSingleHadronNegative->GetBinContent(j+1);
	
		dNCosinedpTSingleHadronPositive[k][j]= HistTempWeightedCosinePositive->GetBinContent(j+1);
		dNCosinedpTSingleHadronNegative[k][j]= HistTempWeightedCosineNegative->GetBinContent(j+1);
		if(dNdpTCountSingleHadronPositive[k][j]> 0.0)
		  {
		    //cross section		    
		    pTHardBinSingleHadronBinError[k][j] = (dNdpTCountSingleHadronPositive[k][j]*HardCrossSection[k]/(Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])))*TMath::Sqrt( (1/dNdpTCountSingleHadronPositive[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
                    cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadronPositive->GetBinContent(j+1)<<"\t Scaled Value = "<<((dNdpTCountSingleHadronPositive[k][j]*HardCrossSection[k])/(Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])))<<endl;
		    cout<<"For SingleHadronV2Bin j = "<<j<<" \t BinContent = \t"<<HistTempWeightedCosinePositive->GetBinContent(j+1)<<"\t Scaled Value = "<<((dNCosinedpTSingleHadronPositive[k][j]*HardCrossSection[k])/(Events[k]*sqrt(SumSquareSoftV2[k]/Events[k])*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])))<<endl;
		    cout<<"Inside: For k="<<k<<",j="<<j<<", Single hadron cross section error="<<pTHardBinSingleHadronBinError[k][j]<<", HardCrossSection[k]"<< HardCrossSection[k]<<endl;
		  }
		else
		  {
		    pTHardBinSingleHadronBinError[k][j] = 0.0;
		    cout<<"Zero Count For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadronPositive->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
		  }
		//cout<<"OutSide: For k="<<k<<",j="<<j<<", Single hadron cross section error="<<pTHardBinSingleHadronBinError[k][j]<<endl;
	      }    
	    //Write histogram into a root file
		HistTempSingleHadronPositive->Write();
		HistTempSingleHadronNegative->Write();
		HistTempWeightedCosinePositive->Write();
		HistTempWeightedCosineNegative->Write();


	    TVector EventInfo(3);
	    EventInfo[0] = HardCrossSection[k];
	    EventInfo[1] = HardCrossSectionError[k];
	    EventInfo[2] = Events[k];
	    EventInfo.Write("EventInfo");	    
 
    //delete outFile;
	EventA<<"(AvgE,AvgPx,AvgPy,AvgPz) = ("<<SumE/Events[k]<<","<< SumPx/Events[k]<<","<<SumPy/Events[k]<<","<<SumPz/Events[k]<<")"<<endl;
	EventR<<"For bin k="<<k<<",pTHat("<< pTHatMin[k]<<","<<pTHatMax[k] <<")GeV"<<"\t total events = "<<Events[k]<<endl; 
    } //k-loop ends here (pTHatBin loop)
    
    delete  HistTempSingleHadronPositive, HistTempSingleHadronNegative;
    delete  HistTempWeightedCosinePositive, HistTempWeightedCosineNegative;
    delete  HistTemp;
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
    double  DifferentialSingleHadronYieldPositive[NpTSingleHadronBin], DifferentialSingleHadronYieldNegative[NpTSingleHadronBin];
    double DifferentialSingleHadronYieldError[NpTSingleHadronBin],SingleHadronpT[NpTSingleHadronBin],SingleHadronpTError[NpTSingleHadronBin];
    double  TotalDifferentialSingleHadronYieldPositive[NpTSingleHadronBin], TotalDifferentialSingleHadronYieldNegative[NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYieldError[NpTSingleHadronBin];
    TGraphErrors * GESingleHadron;
    for(int j=0; j<NpTSingleHadronBin;j++)
      {
        TotalDifferentialSingleHadronYieldPositive[j] = 0.0; TotalDifferentialSingleHadronYieldNegative[j] = 0.0;
        TotalDifferentialSingleHadronYieldError[j] = 0.0;
      }
    
    
    for(int k=0;k<NpTHardBin;k++)
      {cout<<"For ptHardBin = "<<k+1<<"\t SingleHadron differential yield is Below "<<endl;
	for(int j=0; j<NpTSingleHadronBin;j++)
	  {
            
	    DifferentialSingleHadronYieldPositive[j] = (dNdpTCountSingleHadronPositive[k][j]*HardCrossSection[k])/(Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j]));
	    DifferentialSingleHadronYieldNegative[j] = (dNdpTCountSingleHadronNegative[k][j]*HardCrossSection[k])/(Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j]));

	    TotalDifferentialSingleHadronYieldPositive[j] = TotalDifferentialSingleHadronYieldPositive[j] + DifferentialSingleHadronYieldPositive[j];
	    TotalDifferentialSingleHadronYieldNegative[j] = TotalDifferentialSingleHadronYieldNegative[j] + DifferentialSingleHadronYieldNegative[j];
            TotalDifferentialSingleHadronYieldError[j] = TotalDifferentialSingleHadronYieldError[j] + TMath::Power( pTHardBinSingleHadronBinError[k][j], 2.0);
            SingleHadronpT[j] = (SingleHadronpTBin[j+1] + SingleHadronpTBin[j])/2.0;
            SingleHadronpTError[j] = (SingleHadronpTBin[j+1] - SingleHadronpTBin[j])/2.0;
            cout<<SingleHadronpT[j]<<"\t"<<DifferentialSingleHadronYieldPositive[j]<<"\t"<<pTHardBinSingleHadronBinError[k][j]<<endl;
	  }
	//        GESingleHadron = new TGraphErrors(NpTSingleHadronBin,SingleHadronpT,DifferentialSingleHadronYield,SingleHadronpTError,DifferentialSingleHadronYieldError);
        char MyGraphName2[100];
        //sprintf(MyGraphName2,"DifferentialSingleHadronYieldBin%i_%i",pTHatMin[k],pTHatMax[k]);
        //GESingleHadron->SetName(MyGraphName2);
        //GESingleHadron->Write();
      }

    //Error
    for(int j=0;j<NpTSingleHadronBin;j++)
      {
        TotalDifferentialSingleHadronYieldError[j] =TMath::Sqrt(TotalDifferentialSingleHadronYieldError[j]);
      }
    //Compute Error in Numerator of v2 calculation
    double ErrorAverageCosineTerm[NpTHardBin][NpTSingleHadronBin], AverageSquareV2Soft[NpTHardBin];
    //Calculating V2 now
    for(int j=0; j<NpTSingleHadronBin;j++)
      { V2Plus[j]=0.0; V2PlusMinus[j]=0.0;  V2PlusError[j]=0;
	for(int k=0;k<NpTHardBin;k++)
	  {
	    ErrorAverageCosineTerm[k][j] = (SumSquareCosineWeight[k][j] - Events[k]*pow(dNCosinedpTSingleHadronPositive[k][j]/Events[k],2) )/(Events[k] -1);
	    ErrorAverageCosineTerm[k][j] = sqrt(ErrorAverageCosineTerm[k][j]/Events[k])/(SingleHadronpTBin[j+1]-SingleHadronpTBin[j]);
	    AverageSquareV2Soft[k] = sqrt(SumSquareSoftV2[k]/(1.0*Events[k]));
	    
	    V2Plus[j] = V2Plus[j] + ((dNCosinedpTSingleHadronPositive[k][j]*HardCrossSection[k])/(AverageSquareV2Soft[k]*Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*TotalDifferentialSingleHadronYieldPositive[j]));
	    V2PlusMinus[j] = V2PlusMinus[j] + ((dNCosinedpTSingleHadronPositive[k][j] - dNCosinedpTSingleHadronNegative[k][j])*HardCrossSection[k]/(AverageSquareV2Soft[k]*Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*(TotalDifferentialSingleHadronYieldPositive[j]-TotalDifferentialSingleHadronYieldNegative[j])));	    
	    V2PlusError[j] =  V2PlusError[j] + pow(ErrorAverageCosineTerm[k][j]*HardCrossSection[k]/(AverageSquareV2Soft[k]*TotalDifferentialSingleHadronYieldPositive[j]), 2 ) ;
	    if(j==0){cout<<"at k="<<k<<", AverageSquareV2Soft[k]="<<AverageSquareV2Soft[k]<<endl; }
	  }
	V2PlusError[j] = V2PlusError[j] + pow(V2Plus[j]*TotalDifferentialSingleHadronYieldError[j]/TotalDifferentialSingleHadronYieldPositive[j], 2) ;
	V2PlusError[j] =  sqrt( V2PlusError[j]);
	V2PlusMinusError[j] = V2PlusError[j];
	cout<<"For pT hadron bin j="<<j<<",V2Plus[j]="<<V2Plus[j]<<",+- "<<V2PlusError[j]<<", V2PlusMinus[j]="<<V2PlusMinus[j]<<",+- "<<V2PlusMinusError[j]<<endl;
      }
    cout<<"Final results for total differential single hadron yield"<<endl;
    cout<<"SingleHadronpT \t"<<"v2 \t"<<"Error"<<endl;
    foutputPlus<<"#HadronpT \t"<<"v2 \t"<<"pTBinWidth/2.0 \t"<<"v2Error \t"<<"v2Err%"<<endl;
    foutputPlusMinus<<"#HadronpT \t"<<"v2 \t"<<"pTBinWidth/2.0 \t"<<"v2Error \t"<<"v2Err%"<<endl;
    for(int j=0;j<NpTSingleHadronBin;j++)
      {
        cout<<SingleHadronpT[j]<<"\t"<<V2PlusMinus[j]<<" +/- \t"<<V2PlusMinusError[j]<<endl;
	foutputPlus<<SingleHadronpT[j]<<"\t"<<V2Plus[j]<<"\t"<<SingleHadronpTError[j]<<"\t"<<V2PlusError[j]<<"\t"<<V2PlusError[j]*100.0/V2Plus[j]<<endl;
	foutputPlusMinus<<SingleHadronpT[j]<<"\t"<<V2PlusMinus[j]<<"\t"<<SingleHadronpTError[j]<<"\t"<<V2PlusMinusError[j]<<"\t"<<V2PlusMinusError[j]*100.0/V2PlusMinus[j]<<endl;
      }
    cout<<"Create Final TGraphError plot "<<endl;
    // Final Plot and save also
    TGraphErrors *GESingleHadronV2PlusMinusTotal  = new TGraphErrors(NpTSingleHadronBin, SingleHadronpT, V2PlusMinus, SingleHadronpTError, V2PlusMinusError);
    GESingleHadronV2PlusMinusTotal->SetName("V2PlusMinus");
    GESingleHadronV2PlusMinusTotal->Write();    
    foutputPlus.close();foutputPlusMinus.close();

    delete GESingleHadronV2PlusMinusTotal;
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
