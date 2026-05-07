
import os

def move():

    TestOutDIR="/nfs/zfs2/RawData/RUN7_PP_5020GeV/TestOut"
    FinalOutputDIR="/nfs/zfs2/RawData/RUN7_PP_5020GeV/FinalData"

    for i in range(0, 66):
        Bin = ptHardBins(i)
        LastName="Index_" + str(i) + "_Bin" + str(Bin[0]) + "_" + str(Bin[1]) + "_Run"
        HP=['Hadron', 'Parton']
        for Particle in HP:
            InputName=TestOutDIR + "/Jetscape" + Particle + "List" + LastName +"*.out"
            OutputName=FinalOutputDIR + "/Jetscape" + Particle + "ListBin" +str(Bin[0])+ "_" + str(Bin[1]) + ".out"   
            Command="cat " + InputName + " > " + OutputName
            print(Command)
            os.system(Command)
        #Testoutfile and sigmaGen
        TestOutFile=TestOutDIR + "/TestOut" + LastName + "1.dat"
        SigmaFile=FinalOutputDIR + "/SigmaHardBin" + str(Bin[0])+ "_" + str(Bin[1]) + ".out"
        LineNumbers=['5000','10000','40000']
        SigmaGen=''
        SigmaGenErr=''
        for Numbers in LineNumbers:
            Command= "tail -" + Numbers + " " + TestOutFile + " >  " + "Dummy.txt"
            os.system(Command)
            with open("Dummy.txt",'r') as file:
                for line in file:
                    lineNoN = line.strip()
                    LineList = lineNoN.split()
                    if "writersigmaGen" in LineList:
                        SigmaGen = LineList[3]
                    if "writersigmaErr" in LineList:
                        SigmaGenErr = LineList[3]

            if SigmaGen=='':
                break

        print("SigmaGen=",SigmaGen, ", Err=", SigmaGenErr)
        print(SigmaFile)
        with open(SigmaFile,"w") as file:
            Command=SigmaGen + "\t" + SigmaGenErr
            file.write(Command)
              






def ptHardBins(i):
    #Ecm=5.02TeV; size of BinArray is 67
                                                                                                                    
    Bin =[1, 2, 3, 4, 5, 7, 9,  11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130,\
 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600,\
 700, 800, 900,  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2510]

    #Ecm=2.76TeV; size of bins is 55
    #Bin = [1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1380]

    #Ecm=200GeV; size of bins is 24
    #Bin =[1, 2, 3, 4, 5, 7,  9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]
    return Bin[i], Bin[i+1]

if __name__== "__main__":
    #argv =sys.argv                                                                                                 
    move()








