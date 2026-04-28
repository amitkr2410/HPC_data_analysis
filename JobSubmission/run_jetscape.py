#!/usr/bin/env python
import os
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
import SetupXML

def Main():
    JETSCAPEDIR="/nfs/zfs2/RawData/C-SCAPE-COMP"          #Change
    RawOutputDIR="/nfs/zfs2/RawData/RunPP_5020GeV"        #Change
    TestOutDIR= RawOutputDIR + "/TestOut"
    FinalOutputDIR= RawOutputDIR + "/FinalData"
    CodeDIR="/home/amit/CloudComputing/JobSubmission"    #Change this line
    Exec=CodeDIR +"/RUN.sh"
    UserXML="/nfs/zfs2/RawData/jetscape_pp_jetsmith_5020.xml"   #Change

    Command="  mkdir  " + TestOutDIR
    if os.path.isdir(TestOutDIR)==False:
        os.system(Command)

    Command="  mkdir  " + FinalOutputDIR
    if os.path.isdir(FinalOutputDIR)==False:
        os.system(Command)
    
    
    for i in range(0,66):
        for run in range(1,2):
            Bin = ptHardBins(i)    
            print(Bin[0]," ", Bin[1])
            LastName="Index_" + str(i) + "_Bin" + str(Bin[0]) + "_" + str(Bin[1]) + "_Run" + str(run)
            Seed= str(i) + "000" + str(run)
            print("Seed is ",Seed)
            LogFile=TestOutDIR + "/ScreenLog" + LastName
            JobLogFile=TestOutDIR + "/JobLog" + LastName
            JobErrFile=TestOutDIR + "/JobErr" + LastName
            JobName="I" + str(i) + "_R" +str(run) + "_" + LastName
            XML=TestOutDIR + "/jetscape_user_pp_" + LastName + ".xml"
            TestOutFile=TestOutDIR + "/TestOut" + LastName 
            HadronFile=TestOutDIR + "/JetscapeHadronList" + LastName + ".out"
            PartonFile=TestOutDIR + "/JetscapePartonList" + LastName + ".out"
            
            Command = " echo  " + UserXML + "  " + XML +  " " + TestOutFile + " " + HadronFile + " " + PartonFile
            Command = Command + " " + LogFile + " "+ str(Bin[0]) + " " + str(Bin[1])  + " " + JETSCAPEDIR +" "+ TestOutDIR
            Command = Command + " " + FinalOutputDIR + " " + CodeDIR + " "+ Seed   + " "
            Command = Command + " " + JobLogFile + " " + JobErrFile + " "+ JobName + " "
            Command = Command  + " >> " + LogFile
            print(Command)
            os.system(Command)

            Command="cd " + CodeDIR
            os.system(Command)
            XMLARGV=[UserXML, XML, TestOutFile, Bin[0], Bin[1], int(Seed)]
            print(XMLARGV)
            SetupXML.GenerateXML(XMLARGV)
            
            Command = "sbatch   -N 1 -n 1 --mem=4G --time=7-00:00:00 --nodelist=i48c  --job-name " + JobName
            Command = Command +  "  -o " + JobLogFile + " -e " +  JobErrFile + " -- " +  Exec + " " + JETSCAPEDIR 
            Command = Command +  "  " + XML + " " + TestOutFile + " " + HadronFile + " " + PartonFile + " " + LogFile
            print(Command)
            
            os.system(Command)

            Command = Command  + " >> " + LogFile
            os.system(Command)

def ptHardBins(i):
    #Ecm=5.02TeV; size of bins is 67
    Bin =[1, 2, 3, 4, 5, 7, 9,  11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140,                               150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510]

    #Ecm=2.76TeV; size of bins is 55
    #Bin = [1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1380]
    
    #Ecm=200GeV; size of bins is 24
    #Bin =[1, 2, 3, 4, 5, 7,  9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]
    
    return Bin[i], Bin[i+1]



if __name__ == "__main__":
    Main()
    print("Hello, World!")
