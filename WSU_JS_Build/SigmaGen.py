import sys
import os
import subprocess
import math
from decimal import Decimal

#argvs[1]=Full path of Input directory where test_out_*.dat files are located
#argvs[2]=Full path of output directory where you want to store pTHatCrossSection
#argvs[3]=lower pTHatBin Index 
#argvs[4]=Upper pTHatBin Index

def Main(argc,argvs):

    print('argvs[1]:', argvs[1])
    print('argvs[2]:', argvs[2])
    
    InputDir = argvs[1]
    OutputDir = argvs[2]
    
    pthat_bins = GetPtHatBins( 5020 )

    for BinIndex in range(int(argvs[3]),int(argvs[4])):
        CurrentSigma=Decimal(0.0)
        CurrentSigmaError=Decimal(0.0)
        SumSigma=Decimal(0.0)
        SumRelativeError=Decimal(0.0)
        SumRun=Decimal(0.0)
        this_bin = pthat_bins[BinIndex]
        print('Bin=', this_bin[0],' ,',this_bin[1])
        for run in range(0,12):
            TestOutFileName = InputDir + '/TestOutBin'+str(this_bin[0])+'_'+str(this_bin[1])+'_Run'+str(run) +'.dat'
            print ("File exists:"+str(os.path.exists(TestOutFileName)))
            if  os.path.exists(TestOutFileName) == True:
                tail_command = 'tail -50000  '+ TestOutFileName +' | grep sigma ' +' >         dummy.txt' 
                print('Command: ', tail_command)
                os.system(tail_command)
            
                in_file = open("dummy.txt", "r")
                flineAll = in_file.readlines()
                for fline in flineAll:
                    if fline.split()[1] == 'sigmaGen':
                        CurrentSigma = Decimal(fline.split()[2])
                    if fline.split()[1] == 'sigmaErr':
                        CurrentSigmaError = Decimal(fline.split()[2])
                in_file.close()
            SumSigma = SumSigma + CurrentSigma
            SumRelativeError = SumRelativeError + pow(CurrentSigmaError/CurrentSigma,2)
            SumRun = SumRun + Decimal(1.0)
            print('Total Run=',SumRun,', sigmaGen=', CurrentSigma,', sigmaError=',CurrentSigmaError)
        SumRelativeError = math.sqrt(SumRelativeError)
        print("Avg SigmaGen=",SumSigma/SumRun ,',Avg SigmaErr=',SumSigma*Decimal(SumRelativeError)/SumRun) 
        SigmaFileName = OutputDir + '/SigmaHardBin'  + str(this_bin[0]) + '_' + str(this_bin[1])  + '.out' 
        out_file = open(SigmaFileName,'w')
        out_file.write(str(SumSigma/SumRun))
        out_file.write("  ")
        out_file.write( str(SumSigma*Decimal(SumRelativeError)/SumRun) + " \n ")
        out_file.close()
#                for fWord in fline.split():
#                   print(fWord)
                    
                    ##  char InFileName[100000], OutFileName[10000];
                    ##  ifstream fIn; ofstream fOut;
                    ##  sprintf(InFileName, "", );
                    
                    
                    
                    
def GetPtHatBins(eCM):

    if eCM == 5020:

        pthat_bins = [                       
                      [2400,2510], [2200,2400], [2000,2200], [1900,2000], [1800,1900], [1700,1800], [1600,1700], [1500,1600], [1400,1500], [1300,1400], [1200,1300], [1100,1200], [1000,1100], [900, 1000], [800,900], [700, 800], [600, 700], [550, 600], [500, 550], [450, 500], [400, 450], [350, 400], [300, 350], [290, 300], [280,290], [270, 280], [260, 270], [250, 260], [240, 250], [230, 240], [220, 230], [210, 220], [200, 210], [190, 200], [180,190], [170, 180], [160, 170], [150, 160], [140, 150], [130, 140], [120, 130], [110, 120],  [100, 110], [90,100],  [80,90], [70,80], [60, 70], [55, 60], [50, 55],  [45, 50], [40, 45], [35, 40], [30, 35], [25, 30], [20,25], [17,20], [15,17], [13,15], [11,13], [9,11], [7,9], [5,7], [4, 5], [3, 4], [2,3], [1,2] ]

        return pthat_bins
    else:
        print('error')
        exit() 



if __name__ == '__main__':
    argvs = sys.argv
    argc = len(argvs)
    Main(argc,argvs)
