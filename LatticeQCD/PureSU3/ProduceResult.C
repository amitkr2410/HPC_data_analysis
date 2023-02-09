#include"iostream"
#include"fstream"

void ProduceResult(int Nt,int Ns)
{
  //43, 33
  double Beta[33]={4.0, 4.2, 4.4, 4.6, 4.8, 4.9, 4.95, 5.0, 5.05, 5.1, 5.15, 5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 5.5, 5.55, 5.6, 5.65, 5.7, 5.75, 5.8, 5.85, 5.9, 5.95, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0};// 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0};

  char FileName[5000];
  ofstream fout;
  sprintf(FileName,"UnfinishedResult2/FFCorrelatorLO_NLO_Nt%i_Ns%i_wrtBeta.txt",Nt,Ns);
  fout.open(FileName,ios::out);
  fout<<"# Nt="<<Nt<<", Ns="<<Ns<<", SU3 calculation for Beta"<<endl;
  fout<<"#Beta   \t F3iF3iMinusF4iF4i.Re  \t  ErrorF3iF3iMinusF4iF4i.Re  \t     F3iF3iMinusF4iF4i.Im  \t  ErrorF3iF3iMinusF4iF4i.Im  \t F3iF4iPlusF4iF3i.Re \t ErrorF3iF4iPlusF4iF3i.Re  \t      F3iF4iPlusF4iF3i.Im \t    ErrorF3iF4iPlusF4iF3i.Im  \t      F3iDzF3iMinusF4iDzF4i.Re  \t    ErrorF3iDzF3iMinusF4iDzF4i.Re \t  F3iDzF3iMinusF4iDzF4i.Im  \t   ErrorF3iDzF3iMinusF4iDzF4i.Im \t  F3iDzF4iPlusF4iDzF3i.Re  \t     ErrorF3iDzF4iPlusF4iDzF3i.Re  \t  F3iDzF4iPlusF4iDzF3i.Im  \t       ErrorF3iDzF4iPlusF4iDzF3i.Im"<<endl;

  for(int i=0; i<33; i++)
    {
      ifstream f1;
      vector <double> g1, g2, g3, g4, g5, g6, g7, g8;
      double Mean1=0.0, Mean2=0.0, Mean3=0.0, Mean4=0.0, Mean5=0.0, Mean6=0.0, Mean7=0.0,Mean8=0.0;
      double Err1=0.0, Err2=0.0, Err3=0.0, Err4=0.0, Err5=0.0, Err6=0.0, Err7=0.0, Err8=0.0;
      double a=0; cout<<"Started working on file for beta = "<<Beta[i]<<endl;
      for(int Run=1; Run<7; Run++)
	{
	  sprintf(FileName,"ExtraRun24_24/Run%i/FFCorrelatorEachConfigurationDataLO_NLO_Nt%i_Ns%i_BetaTimes100_is_%i.txt",Run,Nt,Ns,int(Beta[i]*100) );
	  //sprintf(FileName,"UnfinishedData/FFCorrelatorEachConfigurationDataLO_NLO_Nt%i_Ns%i_BetaTimes100_is_%i.txt",Nt,Ns,int(Beta[i]*100) );
	  f1.open(FileName,ios::in);
	  string MyString;
	  getline (f1, MyString);
	  getline (f1, MyString);      
	  while(f1>>a)
	    {      
	      f1>>a;g1.push_back(a); f1>>a;g2.push_back(a); f1>>a;g3.push_back(a); f1>>a;g4.push_back(a);
	      f1>>a;g5.push_back(a); f1>>a;g6.push_back(a); f1>>a;g7.push_back(a); f1>>a;g8.push_back(a);
	    }
	  f1.close();
	}
      cout<<"Size of g1="<<int(g1.size())<<" and g8 ="<<int(g8.size())<<endl;
      //Compute Average
      for(int k=0; k< g1.size(); k++)
	{
	  Mean1 = Mean1 + (g1[k]/double(g1.size()));
	  Mean2 = Mean2 + (g2[k]/double(g2.size()));
	  Mean3 = Mean3 + (g3[k]/double(g3.size()));
	  Mean4 = Mean4 + (g4[k]/double(g4.size()));

	  Mean5 = Mean5 + (g5[k]/double(g5.size()));
	  Mean6 = Mean6 + (g6[k]/double(g6.size()));
	  Mean7 = Mean7 + (g7[k]/double(g7.size()));
	  Mean8 = Mean8 + (g8[k]/double(g8.size()));
	}

      //Mean deviation estimate
      for(int k=0; k< int(g1.size()); k++)
	{
	  Err1 =Err1 + pow(Mean1 - g1[k]  ,2.0);
	  Err2 =Err2 + pow(Mean2 - g2[k]  ,2.0);
	  Err3 =Err3 + pow(Mean3 - g3[k]  ,2.0);
	  Err4 =Err4 + pow(Mean4 - g4[k]  ,2.0);

	  Err5 =Err5 + pow(Mean5 - g5[k]  ,2.0);
	  Err6 =Err6 + pow(Mean6 - g6[k]  ,2.0);
	  Err7 =Err7 + pow(Mean7 - g7[k]  ,2.0);
	  Err8 =Err8 + pow(Mean8 - g8[k]  ,2.0);
	}

      for(int k=0; k< g1.size(); k++)
	{
	  Err1 =sqrt( Err1/(g1.size()-1) )/sqrt(g1.size());
	  Err2 =sqrt( Err2/(g2.size()-1) )/sqrt(g2.size());
	  Err3 =sqrt( Err3/(g3.size()-1) )/sqrt(g3.size());
	  Err4 =sqrt( Err4/(g4.size()-1) )/sqrt(g4.size());

	  Err5 =sqrt( Err5/(g5.size()-1) )/sqrt(g5.size());
	  Err6 =sqrt( Err6/(g6.size()-1) )/sqrt(g6.size());
	  Err7 =sqrt( Err7/(g7.size()-1) )/sqrt(g7.size());
	  Err8 =sqrt( Err8/(g8.size()-1) )/sqrt(g8.size());
	}

      g1.resize(0);g2.resize(0);g3.resize(0);g4.resize(0);   
      g5.resize(0);g6.resize(0);g7.resize(0);g8.resize(0);

      fout<<Beta[i]<<"\t"<<Mean1<<"\t"<<Err1<<"\t"<<Mean2<<"\t"<<Err2<<"\t"<<Mean3<<"\t"<<Err3<<"\t"<<Mean4<<"\t"<<Err4<<"\t";
      fout<<Mean5<<"\t"<<Err5<<"\t"<<Mean6<<"\t"<<Err6<<"\t"<<Mean7<<"\t"<<Err7<<"\t"<<Mean8<<"\t"<<Err8<<endl;
    }
  fout.close();
}
