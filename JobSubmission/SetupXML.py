import xml.etree.ElementTree as ET

def GenerateXML(args):
    InputXML = args[0]
    OutputXML  = args[1]
    TestOutName = args[2]
    pMin=args[3]
    pMax=args[4]
    Seed=args[5]

    print("Amit:",args[0]," ", args[1]," ", args[2])
    print("Amit:",args[3]," ", args[4]," ", args[5])
    
    mytree = ET.parse(InputXML)
    myroot = mytree.getroot()
    print(myroot.tag, "\t",myroot.attrib)
    for x in myroot.iter("outputFilename"):
        x.text = TestOutName
        print(x.tag, "\t",x.attrib, "\t", x.text)

    for x in myroot.iter("seed"):
        x.text = str(Seed);
        print(x.tag, "\t",x.attrib, "\t", x.text)
        
    for x in myroot.iter("pTHatMin"):
        x.text = str(pMin);    
        print(x.tag, "\t",x.attrib, "\t", x.text)   
        
    for x in myroot.iter("pTHatMax"):
        x.text = str(pMax);
        print(x.tag, "\t",x.attrib, "\t", x.text)    
        
    mytree.write(OutputXML)
