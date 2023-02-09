void Mean(TGraph *g, double *Ymean, double *Yerr, int FirstPoint, int SkipN)
{
  double *y= g->GetY(); int n= g->GetN();
  double Ysum=0; int N=0;
  double mean=0.0; double err=0.0;
 
  for(int i=FirstPoint; i<n-1; i++)
    {
      Ysum = Ysum + y[i];//cout<<"i="<<i<<endl;
      N= N+1; i = i+SkipN;
    }
  mean = Ysum*1.0/N;
  for(int i=FirstPoint; i<n-1; i++)
    {
      err = err + pow(y[i]-mean, 2.0);
      i = i+SkipN;
    }
  err = sqrt(err/(N-1))/sqrt(N);
  *Ymean=mean; *Yerr=err;
  //  std::setprecision(6);
  cout<<"Ymean = "<<mean<<"+/-"<<err<<"\t \t total points n="<<n<<"\t total effective points N="<<N<<endl;
}

void JackknifeMean(TGraph *g, double *mean, double *err, int FirstRun, int SkipN)
{
  int n= g->GetN();  double *y= g->GetY(); 
  int BinSize=SkipN+1; int Index=0; const int N= n - FirstRun -5;
  double yall[N];
  const int NewSize=N/BinSize; double Y[NewSize], YSum;
  cout<<"n="<<n<<", FirstRun="<<FirstRun<<", N="<<N<<", BinSize="<<BinSize<<", NewSize="<<NewSize<<"\n";

  for(int k=0; k<N; k++)
    {yall[k]=y[k+FirstRun];
    }

  for(int i=0; i<NewSize; i++)
    {  Y[i]=0;
      for(int j=0; j<BinSize;j++)
	{ Index = j +(i*BinSize);
	  Y[i]= Y[i] + yall[Index];
	}
      Y[i] = Y[i]/BinSize;      
    }
    
  YSum=0; 
  for(int i=0; i<NewSize; i++)
    {  YSum = YSum +  Y[i];
    }
  
  double NaiveMean=0, NaiveError=0; NaiveMean = YSum/NewSize;
  for(int i=0; i<NewSize; i++)
    { NaiveError = NaiveError + TMath::Power(Y[i] - NaiveMean, 2.0);
    }
  NaiveError = TMath::Sqrt(NaiveError/(NewSize*(NewSize-1)));
  cout<<"Naive LO (Mean, Error)=("<<NaiveMean<<","<<NaiveError<<") \n";

  double JackknifedMeanY[NewSize]; 
  double JackknifedMean=0, JackknifedMeanError=0;
  
  for(int i=0; i<NewSize; i++)
    {
      JackknifedMeanY[i]= (YSum - Y[i])*1.0/(NewSize-1);  
      JackknifedMean = JackknifedMean + JackknifedMeanY[i];
    }
  JackknifedMean = JackknifedMean*1.0/NewSize;

  for(int i=0; i<NewSize; i++)
    {
      JackknifedMeanError = JackknifedMeanError + TMath::Power(JackknifedMeanY[i] - JackknifedMean, 2.0);
    }      
  JackknifedMeanError = TMath::Sqrt((NewSize-1)*JackknifedMeanError/NewSize);

  *mean=JackknifedMean; *err=JackknifedMeanError;
  cout<<" Jackknife Average: (Mean,Error)=("<<mean[0]<<","<<err[0]<<")\n";
}

//Ends here

void JackknifeMeanD4z(TGraph *g, TGraph *g1, TGraph *g2, double TadpoleFactorU0, double *mean, double *err, int FirstRun, int SkipN)
{
  int n0= g->GetN(); int n1= g1->GetN(); int n2=g2->GetN(); double *y= g->GetY(); double *y1= g1->GetY(); double *y2= g2->GetY();  
  int n=n0<n1 ? n0:n1;
  n = n<n2 ? n:n2;

  int BinSize=SkipN+1; int Index=0; const int N= n - FirstRun -5;
  double yall[3][N];
  const int NewSize=N/BinSize; double Y[3][NewSize], YSum[3];
  cout<<"n="<<n<<", FirstRun="<<FirstRun<<", N="<<N<<", BinSize="<<BinSize<<", NewSize="<<NewSize<<"\n";
  cout<<"n0="<<n0<<",n1="<<n1<<", n2="<<n2<<"\n";

  for(int k=0; k<N; k++)
    {yall[0][k]=y[k+FirstRun];    yall[1][k]=y1[k+FirstRun];    yall[2][k]=y2[k+FirstRun];
    }

  for(int i=0; i<NewSize; i++)
    {  for(int m=0; m<3; m++)	
	{ Y[m][i]=0;
	  for(int j=0; j<BinSize;j++)
	    { Index= j +(i*BinSize);	
	      Y[m][i]= Y[m][i] + yall[m][Index];
	    }
	  Y[m][i] = Y[m][i]/BinSize;      
	}
    }

  for(int m=0; m<3; m++)
    { YSum[m]=0; 
      for(int i=0; i<NewSize; i++)
	{  
	  YSum[m] = YSum[m] +  Y[m][i]; 
	}
    }
  
  double NaiveMean=0, NaiveError=0; NaiveMean = YSum[0]/NewSize;
  for(int i=0; i<NewSize; i++)
    { NaiveError = NaiveError + TMath::Power(Y[0][i] - NaiveMean, 2.0);
    }
  NaiveError = TMath::Sqrt(NaiveError/(NewSize*(NewSize-1)));
  cout<<"Naive LO (Mean, Error)=("<<NaiveMean<<","<<NaiveError<<") \n";
  
  double JackknifedMeanY[3][NewSize]; 
  double JackknifedFmunu[3][NewSize], JackknifedFmunuMean[3], JackknifedFmunuError[3];
  
  for(int m=0; m<3; m++)
    { for(int i=0; i<NewSize; i++)
	{
	  JackknifedMeanY[m][i]= (YSum[m] - Y[m][i])*1.0/(NewSize-1);  
	}
      JackknifedFmunuMean[m]=0.0;
      JackknifedFmunuError[m]=0.0;
    }

  for(int i=0; i<NewSize; i++)
    {
      JackknifedFmunu[0][i]= JackknifedMeanY[0][i];
      JackknifedFmunu[1][i]= -2.0*JackknifedMeanY[0][i] + (JackknifedMeanY[1][i]/TMath::Power(TadpoleFactorU0,2));
      JackknifedFmunu[2][i]= 6.0*JackknifedMeanY[0][i] + (-4.0*JackknifedMeanY[1][i]/TMath::Power(TadpoleFactorU0,2)) + (JackknifedMeanY[2][i]/TMath::Power(TadpoleFactorU0,4));
      for(int m=0; m<3; m++)
	{ JackknifedFmunuMean[m] = JackknifedFmunuMean[m] + JackknifedFmunu[m][i];
	}
    }
  
  for(int m=0; m<3; m++)
    {JackknifedFmunuMean[m] = JackknifedFmunuMean[m]/NewSize;
      for(int i=0; i<NewSize; i++)
	{
	  JackknifedFmunuError[m] = JackknifedFmunuError[m] + TMath::Power(JackknifedFmunu[m][i] - JackknifedFmunuMean[m], 2.0);
	}      
      JackknifedFmunuError[m] = TMath::Sqrt((NewSize-1)*JackknifedFmunuError[m]/NewSize);
    }
  

  for(int m=0; m<3; m++)
    { mean[m] = JackknifedFmunuMean[m];
      err[m]  = JackknifedFmunuError[m];
    }  
  cout<<" Jackknife Average: FmunuLO, FmunuNNLO, FmunuNNNNLO=("<<mean[0]<<","<<err[0]<<"),\t("<<mean[1]<<","<<err[1]<<"),\t("<<mean[2]<<","<<err[2]<<")\n";
}

TGraphErrors * ScaleGraph(TGraphErrors *g, double PreFactor)
{

  const int N=g->GetN(); double *X=g->GetX();double *Y=g->GetY();double *EY=g->GetEY(); double *EX=g->GetEX();

  double x[N], y[N], ex[N], ey[N];
  for(int i=0;i<N;i++)
    {
      x[i]=X[i];
      ex[i]=EX[i];//cout<<"X[i] = "<<endl;
      y[i]  = Y[i]*PreFactor;
      ey[i] = EY[i]*PreFactor;
    }
  TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
  return gr;
}

TGraphErrors * SubtractTGraphErrors(TGraphErrors *g1, TGraphErrors *g2)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  double Diff[N1], DiffError[N1],Zero[N1];

  if(N1==N2)
    {   for(int i=0;i<N1;i++)
	{cout<<"Amit X1,X2 = "<<X1[i]<<","<<X2[i]<<endl;
	  Diff[i] = Y1[i]-Y2[i];cout<<"Subtract Y1,Y2 = "<<Y1[i]<<","<<Y2[i]<<", Diff="<<Diff[i]<<endl;
	  DiffError[i] = TMath::Sqrt( pow(EY1[i], 2.0) + pow(EY2[i], 2.0) );
	  Zero[i]=0.0;
	}
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;}
  cout<<"\n\n";
  TGraphErrors *G = new TGraphErrors(N1,X1,Diff,Zero,DiffError);
  return G;
}


TGraphErrors * DivideTGraphErrors(TGraphErrors* Numerator, double powNum, TGraphErrors* Denominator, double powDen)
{
  int N1=Numerator->GetN(); double *X1=Numerator->GetX();double *Y1=Numerator->GetY();double *EY1=Numerator->GetEY(); double *EX1=Numerator->GetEX();
  int N2=Denominator->GetN(); double *X2=Denominator->GetX();double *Y2=Denominator->GetY();double *EY2=Denominator->GetEY(); double *EX2=Denominator->GetEX();
  double Ratio[N1], RatioError[N1],Zero[N1];
 
  if(N1==N2)
    {   for(int i=0;i<N1;i++)
          {cout<<"Amit X1,X2 = "<<X1[i]<<","<<X2[i]<<endl;
	    if(Y2[i]!=0)
	      {
	    Ratio[i] = TMath::Power(Y1[i], powNum)/TMath::Power(Y2[i], powDen);cout<<"Divide Y1,Y2 = pow("<<Y1[i]<<","<<powNum<<"), pow("<<Y2[i]<<","<<powDen<<"), Ratio="<<Ratio[i]<<endl;
	    RatioError[i] = Ratio[i]*TMath::Sqrt( pow(powNum*EY1[i]/Y1[i],2.0) + pow(powDen*EY2[i]/Y2[i],2.0) );	    
	      }
	    else 
	      {
		Ratio[i] = 0;  cout<<"Division by zero XX"<<endl;
		RatioError[i] =0;
	      }
	    Zero[i]=0.0;
	  }      
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
  cout<<"\n\n";
  TGraphErrors *G = new TGraphErrors(N1,X2,Ratio,Zero,RatioError);
  return G;

}

void WriteToFile(TGraphErrors *g1, TGraphErrors *g2, string OutputFileS)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();

  ofstream f; f.open(OutputFileS, std::ios::out);
  if(N1==N2)
    {
      for(int i=0; i<N1; i++)
        {
          f<<X1[i]<<"\t"<<Y1[i]<<"\t"<<EY1[i]<<"\t"<<Y2[i]<<"\t"<<EY2[i]<<endl;
        }
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
  f.close();
}

void WriteToFile(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, string OutputFileS)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  ofstream f; f.open(OutputFileS, std::ios::out);
  if(N1==N2 && N1==N3)
    {
      for(int i=0; i<N1; i++)
        {
          f<<X1[i]<<"\t"<<Y1[i]<<"\t"<<EY1[i]<<"\t"<<Y2[i]<<"\t"<<EY2[i]<< "\t"<<Y3[i]<<"\t"<<EY3[i] <<endl;
        }
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
  f.close();
}

void WriteToFile(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, string OutputFileS)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  int N4=g4->GetN(); double *X4=g4->GetX();double *Y4=g4->GetY();double *EY4=g4->GetEY(); double *EX4=g4->GetEX();

  ofstream f; f.open(OutputFileS, std::ios::out);
  if(N1==N2 && N3==N4)
    {
      for(int i=0; i<N1; i++)
	{
	  f<<X1[i]<<"\t"<<Y1[i]<<"\t"<<EY1[i]<<"\t"<<Y2[i]<<"\t"<<EY2[i]<<"\t"<<Y3[i]<<"\t"<<EY3[i]<<"\t"<<Y4[i]<<"\t"<<EY4[i]<<endl;
	}
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
  f.close();
}

void WriteToFile(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, string OutputFileS)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  int N4=g4->GetN(); double *X4=g4->GetX();double *Y4=g4->GetY();double *EY4=g4->GetEY(); double *EX4=g4->GetEX();
  int N5=g5->GetN(); double *X5=g5->GetX();double *Y5=g5->GetY();double *EY5=g5->GetEY(); double *EX5=g5->GetEX();
  int N6=g6->GetN(); double *X6=g6->GetX();double *Y6=g6->GetY();double *EY6=g6->GetEY(); double *EX6=g6->GetEX();
  ofstream f; f.open(OutputFileS, std::ios::out);
  if(N1==N2 && N3==N4)
    {
      for(int i=0; i<N1; i++)
        {
          f<<X1[i]<<"\t"<<Y1[i]<<"\t"<<EY1[i]<<"\t"<<Y2[i]<<"\t"<<EY2[i]<<"\t"<<Y3[i]<<"\t"<<EY3[i]<<"\t"<<Y4[i]<<"\t"<<EY4[i]<<"\t"<<Y5[i]<<"\t"<<EY5[i]<<"\t"<<Y6[i]<<"\t"<<EY6[i] <<endl;
	}
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
}

void WriteToFile(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, TGraphErrors *g7, TGraphErrors *g8,  string OutputFileS)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  int N4=g4->GetN(); double *X4=g4->GetX();double *Y4=g4->GetY();double *EY4=g4->GetEY(); double *EX4=g4->GetEX();
  int N5=g5->GetN(); double *X5=g5->GetX();double *Y5=g5->GetY();double *EY5=g5->GetEY(); double *EX5=g5->GetEX();
  int N6=g6->GetN(); double *X6=g6->GetX();double *Y6=g6->GetY();double *EY6=g6->GetEY(); double *EX6=g6->GetEX();
  int N7=g7->GetN(); double *X7=g7->GetX();double *Y7=g7->GetY();double *EY7=g7->GetEY(); double *EX7=g7->GetEX();
  int N8=g8->GetN(); double *X8=g8->GetX();double *Y8=g8->GetY();double *EY8=g8->GetEY(); double *EX8=g8->GetEX();
  ofstream f; f.open(OutputFileS, std::ios::out);
  if(N1==N2 && N3==N4)
    {
      for(int i=0; i<N1; i++)
        {
          f<<X1[i]<<"\t"<<Y1[i]<<"\t"<<EY1[i]<<"\t"<<Y2[i]<<"\t"<<EY2[i]<<"\t"<<Y3[i]<<"\t"<<EY3[i]<<"\t"<<Y4[i]<<"\t"<<EY4[i]<<"\t"<<Y5[i]<<"\t"<<EY5[i]<<"\t"<<Y6[i]<<"\t"<<EY6[i]<<"\t"<<Y7[i]<<"\t"<<EY7[i]<<"\t"<<Y8[i]<<"\t"<<EY8[i] <<endl;
        }
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
  f.close();
}
double * FindMinMax8(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, TGraphErrors *g7, TGraphErrors *g8)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();
  int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();
  int N4=g4->GetN(); double *X4=g4->GetX();double *Y4=g4->GetY();
  int N5=g5->GetN(); double *X5=g5->GetX();double *Y5=g5->GetY();
  int N6=g6->GetN(); double *X6=g6->GetX();double *Y6=g6->GetY();
  int N7=g7->GetN(); double *X7=g7->GetX();double *Y7=g7->GetY();
  int N8=g8->GetN(); double *X8=g8->GetX();double *Y8=g8->GetY();
  double Min=0, Max=0;static double MinMax[2];
  Min=Y1[0]; Max=Y1[0];
  for(int i=0; i<N1; i++)
    { if(Min>Y1[i]){Min=Y1[i];}
      if(Max<Y1[i]){Max=Y1[i];}
    }

  for(int i=0; i<N2; i++)
    { if(Min>Y2[i]){Min=Y2[i];}
      if(Max<Y2[i]){Max=Y2[i];}
    }

  for(int i=0; i<N3; i++)
    { if(Min>Y3[i]){Min=Y3[i];}
      if(Max<Y3[i]){Max=Y3[i];}
    }

  for(int i=0; i<N4; i++)
    { if(Min>Y4[i]){Min=Y4[i];}
      if(Max<Y4[i]){Max=Y4[i];}
    }
  for(int i=0; i<N5; i++)
    { if(Min>Y5[i]){Min=Y5[i];}
      if(Max<Y5[i]){Max=Y5[i];}
    }
  for(int i=0; i<N6; i++)
    { if(Min>Y6[i]){Min=Y6[i];}
      if(Max<Y6[i]){Max=Y6[i];}
    }

  for(int i=0; i<N7; i++)
    { if(Min>Y7[i]){Min=Y7[i];}
      if(Max<Y7[i]){Max=Y7[i];}
    }
  for(int i=0; i<N8; i++)
    { if(Min>Y8[i]){Min=Y8[i];}
      if(Max<Y8[i]){Max=Y8[i];}
    }
  cout<<"Min,Max="<<Min<<","<<Max<<endl;
  MinMax[0]=Min; MinMax[1]=Max; return MinMax;
}

double * FindMinMax6(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  int N4=g4->GetN(); double *X4=g4->GetX();double *Y4=g4->GetY();double *EY4=g4->GetEY(); double *EX4=g4->GetEX();
  int N5=g5->GetN(); double *X5=g5->GetX();double *Y5=g5->GetY();double *EY5=g5->GetEY(); double *EX5=g5->GetEX();
  int N6=g6->GetN(); double *X6=g6->GetX();double *Y6=g6->GetY();double *EY6=g6->GetEY(); double *EX6=g6->GetEX();
  double Min=0, Max=0;static double MinMax[2];
  Min=Y1[0]; Max=Y1[0];
  for(int i=0; i<N1; i++)
    { if(Min>Y1[i]){Min=Y1[i];}
      if(Max<Y1[i]){Max=Y1[i];}
    }
  
  for(int i=0; i<N2; i++)
      { if(Min>Y2[i]){Min=Y2[i];}
	if(Max<Y2[i]){Max=Y2[i];}
      }
    
    for(int i=0; i<N3; i++)	 
      { if(Min>Y3[i]){Min=Y3[i];}
	if(Max<Y3[i]){Max=Y3[i];}
      }
    
    for(int i=0; i<N4; i++)
      { if(Min>Y4[i]){Min=Y4[i];}
	if(Max<Y4[i]){Max=Y4[i];}
      }

    for(int i=0; i<N5; i++)
      { if(Min>Y5[i]){Min=Y5[i];}
        if(Max<Y5[i]){Max=Y5[i];}
      }
    for(int i=0; i<N6; i++)
      { if(Min>Y6[i]){Min=Y6[i];}
        if(Max<Y6[i]){Max=Y6[i];}
      }
    cout<<"Min,Max="<<Min<<","<<Max<<endl;
    MinMax[0]=Min; MinMax[1]=Max; return MinMax;
}

double * FindMinMax4(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  int N3=g3->GetN(); double *X3=g3->GetX();double *Y3=g3->GetY();double *EY3=g3->GetEY(); double *EX3=g3->GetEX();
  int N4=g4->GetN(); double *X4=g4->GetX();double *Y4=g4->GetY();double *EY4=g4->GetEY(); double *EX4=g4->GetEX();
  double Min=0, Max=0;static double MinMax[2];
  Min=Y1[0]; Max=Y1[0];
  for(int i=0; i<N1; i++)
    { if(Min>Y1[i]){Min=Y1[i];}
      if(Max<Y1[i]){Max=Y1[i];}
    }

  for(int i=0; i<N2; i++)
    { if(Min>Y2[i]){Min=Y2[i];}
      if(Max<Y2[i]){Max=Y2[i];}
    }

  for(int i=0; i<N3; i++)
    { if(Min>Y3[i]){Min=Y3[i];}
      if(Max<Y3[i]){Max=Y3[i];}
    }

  for(int i=0; i<N4; i++)
    { if(Min>Y4[i]){Min=Y4[i];}
      if(Max<Y4[i]){Max=Y4[i];}
    }
  cout<<"Min,Max="<<Min<<","<<Max<<endl;
  MinMax[0]=Min; MinMax[1]=Max; return MinMax;
}

double * FindMinMax2(TGraphErrors *g1, TGraphErrors *g2)
{
  int N1=g1->GetN(); double *X1=g1->GetX();double *Y1=g1->GetY();double *EY1=g1->GetEY(); double *EX1=g1->GetEX();
  int N2=g2->GetN(); double *X2=g2->GetX();double *Y2=g2->GetY();double *EY2=g2->GetEY(); double *EX2=g2->GetEX();
  double Min=0, Max=0;static double MinMax[2];
  Min=Y1[0]; Max=Y1[0];
  for(int i=0; i<N1; i++)
    { if(Min>Y1[i]){Min=Y1[i];}
      if(Max<Y1[i]){Max=Y1[i];}
    }

  for(int i=0; i<N2; i++)
    { if(Min>Y2[i]){Min=Y2[i];}
      if(Max<Y2[i]){Max=Y2[i];}
    }

  cout<<"Min,Max="<<Min<<","<<Max<<endl;
  MinMax[0]=Min; MinMax[1]=Max; return MinMax;
}
