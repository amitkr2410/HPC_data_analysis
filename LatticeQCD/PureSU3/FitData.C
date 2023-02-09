void FitData()
{

  TGraphErrors *g = new TGraphErrors("DataTMinusV6_24.dat","%lf %lf %lf");
  g->Draw("A*");
  

  //TF1 *myfit = new TF1("myfit","[0]+[1]*x + [2]*x*x+[3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x ");
  //TF1 *myfit = new TF1("myfit","[0]+[1]*x + [2]*x*x + [3]*x*x*x + [4]*exp(-0.5*pow(((x-[5])/[6]),2.0))",6,8.0,7);
  TF1 *myfit = new TF1("myfit","[0]+[1]*x + [2]*x*x ",6.0,8.0,3);
  myfit->SetParameter(0, 0.2);
  myfit->SetParameter(1, 0.05);
  myfit->SetParameter(2, 10.2);
  // myfit->SetParameter(3, 0.05);
  // myfit->SetParameter(4, 5.2);
  //myfit->SetParameter(5, 6.2);
  //myfit->SetParameter(6, 1.0);

  g->Fit("myfit","R");
    /*
    double *a = myfit->GetParameters(); 
    double Beta[33]={4.0, 4.2, 4.4, 4.6, 4.8, 4.9, 4.95, 5.0, 5.05, 5.1, 5.15, 5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 5.5, 5.55, 5.6, 5.65, 5.7, 5.75, 5.8, 5.85, 5.9, 5.95, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0};
    for(int i=0; i<33; i++)
      {
	double x=Beta[i];
	double c=((rand()%10)/2.0) + ((rand()%1)/1.0) + ((rand()%1)/1.0);
	double y= a[0]+a[1]*x + a[2]*x*x+a[3]*x*x*x + a[4]*x*x*x*x + a[5]*x*x*x*x*x;// + a[6]*x*x*x*x*x*x + a[7]*x*x*x*x*x*x*x + a[8]*x*x*x*x*x*x*x*x;
	cout<<Beta[i]<<"\t"<<y<<"\t"<<(6.3-c)*pow(10,-5.0)<<endl;
      }	
    */

	
}
