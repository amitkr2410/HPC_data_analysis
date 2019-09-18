#include <iostream>
#include <fstream>
#include <string>

int GenerateXML()
{

  int ValuepTHatMin[37] = {40, 45, 50, 55, 60, 65, 70, 75,    80, 90,  100,  110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600};
  int ValuepTHatMax[37] = {45, 50, 55, 60, 65, 70, 75, 80,    90, 100, 110,  120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 2760};

  ifstream ReadFile;
  ofstream OutFile;
  char Name1[100], Name2[100], XMLFileName[100];
  
  for(int i=0;i<3;i++)
    {
      sprintf(XMLFileName,"InputFile/jetscape_initBin%i_%i.xml",ValuepTHatMin[i] ,ValuepTHatMax[i] );
      ReadFile.open("JETSPCAE/examples/jetscape_init.xml");
      OutFile.open(XMLFileName,ios::out);
      sprintf(Name1, "\t\t<pTHatMin>%i</pTHatMin>",ValuepTHatMin[i]);
      sprintf(Name2, "     <pTHatMax>%i</pTHatMax>",ValuepTHatMax[i]);
      stringstream SS1,SS2;
      string Replace1, Replace2;
      SS1 << Name1; SS1 >> Replace1; 
      SS2 << Name2; SS2 >> Replace2; 
      
      string ReadOut="";
      string Search1="       <pTHatMin>110</pTHatMin>";
      string Search2="       <pTHatMax>120</pTHatMax>";
      
      while(getline(ReadFile,ReadOut))
	{
	  if(ReadOut == Search1)
	    {
	      OutFile << Replace1<<endl;
	      //cout<<ReadOut<<endl;
	    }
	  else
	    {
	      if(ReadOut == Search2)
		{
		  OutFile << Replace2<<endl;
		  //cout<<ReadOut<<endl;
		}
	      else
		{
		  OutFile << ReadOut<<endl;
		  //cout<<ReadOut<<endl;
		}
	    }
	}
      OutFile.close();ReadFile.close();
    }
  return 0;
  
}
