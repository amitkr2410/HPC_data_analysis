#!/usr/bin/env python
import os
def Main():
    JETSCAPEDIR="/scratch/amitkr/JSWork/JETSCAPE"
    TestOutDIR="/scratch/amitkr/JSWork/RawFile"
    FinalOutputDIR="/scratch/amitkr/JSWork/FinalData" 
    CodeDIR="/home/amitkr"
    Exec=CodeDIR +"/RUN.sh"
    UserXML="/home/amitkr/jetscape_user_pbpb_grid.xml"
    
    for i in range(0,23):
        for run in range(0,20):        
            Bin = ptHardBins(i)    
            print(Bin[0]," ", Bin[1])
            LastName="Bin"+	str(Bin[0])+"_"+str(Bin[1])+ "_Run" + str(run)
            LogFile=TestOutDIR + "/ScreenLogIndex" + str(i) + "_Run" +str(run) + "_"+ LastName
            JobLogFile=TestOutDIR + "/JobLogIndex" + str(i) + "_Run" +str(run) + "_"+ LastName
            JobErrFile=TestOutDIR + "/JobErrIndex" + str(i) + "_Run" +str(run) + "_"+ LastName
            JobName="I" + str(i) + "_Run" +str(run) + "_"+ LastName
            XML=TestOutDIR + "/jetscape_user_pbpb_" + LastName + ".xml"
            TestOutFile=TestOutDIR + "/TestOut" + LastName 
            HadronFile=TestOutDIR + "/JetscapeHadronList" + LastName + ".out"
            PartonFile=TestOutDIR + "/JetscapePartonList" + LastName + ".out"
            Command="sbatch   -N 1 -n 1 --mem=4G --time=10:00:00 --account=rrg-jeon-ac  --job-name " + JobName
            Command = Command +  "  -o " + JobLogFile + " -e " +  JobErrFile + " -- " +  Exec + " " + UserXML
            Command = Command +  " "+ XML + " " + TestOutFile + " " + HadronFile +" "+ PartonFile + " " + LogFile
            Command = Command +  " "+ str(Bin[0]) + " " + str(Bin[1])
            Command = Command + " "+ JETSCAPEDIR +" "+  TestOutDIR + " "+ FinalOutputDIR + " "+ CodeDIR 
            print(Command)
            
            os.system(Command)


def ptHardBins(i):

    Bin =[1, 2, 3, 4, 5, 7,  9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]
    
    return Bin[i], Bin[i+1]



if __name__ == "__main__":
    Main()
    print("Hello, World!")
