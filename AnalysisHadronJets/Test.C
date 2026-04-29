#include "string"
#include <iostream>
#include <fstream>
#include "iomanip"

#include <memory>
#include <chrono>
#include <thread>
#include "/nfs/zfs2/Software/pythia8309/include/Pythia8/Pythia.h"
#include "/nfs/zfs2/Software/pythia8309/include/Pythia8/FJcore.h"

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
int main(int argc, char **argv)
{  
  Pythia8::Pythia pythia;
  double charge = pythia.particleData.charge( 2212 );
  cout<<"Charge of proton is "<<charge<<endl;
  TH1D *HistTempSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", 10,0.0,10.0);
  cout<<"Spectrum is created "<<endl;
  return 0;
}
