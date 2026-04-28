import os
import sys
#import SetupXML
from pathlib import Path
#import pandas as pd
import glob
import xml.etree.ElementTree as ET
from datetime import datetime

def OrganizeData():
    TestOutDIR="/nfs/zfs2/RawData/RunPP_5020GeV/TestOut"
    FinalOutputDIR="/nfs/zfs2/RawData/RunPP_5020GeV/FinalData"    
    # 5020GeV (0,66)
    Command= TestOutDIR +"/*.xml"
    xml_files = glob.glob(Command)
    for xml in xml_files:
        print(xml)
        XML=xml
        LastName = XML.split("/jetscape_user_pp_")[1].split('.xml')[0]
        print(LastName)
        #XML=TestOutDIR + "/jetscape_user_pp_" + LastName + ".xml"
        TestOutFile=TestOutDIR + "/TestOut" + LastName 
        HadronFile=TestOutDIR + "/JetscapeHadronList" + LastName + ".out"
        PartonFile=TestOutDIR + "/JetscapePartonList" + LastName + ".out"

        FlagFile=Path(HadronFile).exists()
        if FlagFile==False:
            print("file ", HadronFile, ' does not exist')
            continue
        
        mytree = ET.parse(XML)
        myroot = mytree.getroot()
        for x in myroot.iter("nEvents"):
            NeventRequested = x.text
            
        Command="grep Event " +  HadronFile  + " > testFile.txt"
        print(Command)
        os.system(Command)
        Command="wc -l testFile.txt > count.txt"
        print(Command)
        os.system(Command)
        #NeventGenerated=list("count.txt")[0]
        f = open("count.txt")
        lines = f.readlines()
        f.close()
        print(lines)
        print(lines[0])
        NeventGenerated=lines[0]
        #df = pd.read_csv(HadronFile, sep="\s+", header=None)
        #df_new = df.loc[df[0] == "#"]    
        #NeventGenerated = len(df_new)
        Command="echo Number of events matches " + str(NeventRequested) + " " + str(NeventGenerated)
        #print(Command)        
        #if int(NeventRequested) == int(NeventGenerated):
        #        Command = "cp " + HadronFile + " " + FinalOutputDIR + "/"                
        #        os.system(Command)
        #        Command = "cp " + PartonFile + " " + FinalOutputDIR + "/"
        #        os.system(Command)

    EndTime=datetime.now()
    Command = "echo starttime is " + str(StartTime) + " and endtime is "+ str(EndTime)
    

def ptHardBins(i):
    #Ecm=5.02TeV; size of BinArray is 67                                                                                      
    Bin =[1, 2, 3, 4, 5, 7, 9,  11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510]

    #Ecm=2.76TeV; size of bins is 55                                                                                      
    #Bin = [1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1380]                                                                                                          

    #Ecm=200GeV; size of bins is 24                                                                                       
    #Bin =[1, 2, 3, 4, 5, 7,  9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]                     

    return Bin[i], Bin[i+1]

if __name__== "__main__":
    #argv =sys.argv
    OrganizeData()


