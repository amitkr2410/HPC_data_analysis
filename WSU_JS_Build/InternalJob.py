import os
import sys
import SetupXML
import pandas as pd
import xml.etree.ElementTree as ET
from datetime import datetime

def GenerateData(argv):
    StartTime=datetime.now()
    
    InputXML=argv[1]
    OutputXML=argv[2]
    TestOutName=argv[3]
    HadronFile=argv[4]
    PartonFile=argv[5]
    LogFile=argv[6]
    pTHard=(argv[7],argv[8])
    JETSCAPEDIR=argv[9]
    TestOutDIR=argv[10]
    FinalOutputDIR=argv[11]
    CodeDIR=argv[12]
    BuildDIR=argv[13]

    for x in argv:
        Command="echo " + x + " >> " +  LogFile
        os.system(Command)

    Command="cd " + CodeDIR
    os.system(Command)
    XMLARGV=[InputXML, OutputXML, TestOutName, pTHard[0], pTHard[1]]    
    SetupXML.GenerateXML(XMLARGV)

    MainXML=JETSCAPEDIR + "/" + "config/jetscape_main.xml"
    Command="mkdir " + BuildDIR +" ;"
    Command=Command + " cd " + BuildDIR + " ;"
    Command=Command + " cmake " + JETSCAPEDIR + " -DUSE_MUSIC=ON -DUSE_ISS=ON ;" + " make -j 8 ;" 
    #Command="cd " + JETSCAPEDIR + "/build ;"
    Command=Command + "./runJetscape " + OutputXML + " " + MainXML  + " >> " + LogFile + " ;"
    Command=Command + "./FinalStateHadrons " + TestOutName + ".dat  " + HadronFile + " >> " + LogFile + " ;"
    Command=Command + "./FinalStatePartons " + TestOutName + ".dat  " + PartonFile + " >> " + LogFile    + " ;"
    print(Command)
    os.system(Command)

    mytree = ET.parse(OutputXML)
    myroot = mytree.getroot()
    for x in myroot.iter("nEvents"):
        NeventRequested = x.text
    
    df = pd.read_csv(HadronFile, sep="\s+", header=None)
    df_new = df.loc[df[0] == "#"]    
    NeventGenerated = len(df_new)
    Command="echo Number of events matches " + str(NeventRequested) + " " + str(NeventGenerated) + " >> " + LogFile
    os.system(Command)
    if int(NeventRequested) == int(NeventGenerated):
        Command = "cp " + HadronFile + " " + FinalOutputDIR + "/"
        os.system(Command)

    EndTime=datetime.now()
    Command = "echo starttime is " + str(StartTime) + " and endtime is "+ str(EndTime)  + " >> " + LogFile
    os.system(Command)
    
if __name__== "__main__":
    argv =sys.argv
    GenerateData(argv)


