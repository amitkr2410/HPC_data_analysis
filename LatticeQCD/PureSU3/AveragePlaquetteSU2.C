//This program calculates the average plaquette for SU(2) gauge field for N*N*N*N Lattice with Periodic BC's
//  Following is lattice structure for 3*4 (Notice the missing link at boundary;they are to incorporate the periodic BC )
//
//     |_|_|_|_
//     |_|_|_|_
//     |_|_|_|_
//
//We use cpp function "rand()" to generate random number
//We use 5 variables to store the coordinate of the link (first four for t,x,y,z and last-one for its direction; We use 0 to point in t-dir, 1 (x-dir).. ) 
//We use Heat-bath algorithm to update each link
//To compile this code (for beta =2.3) in terminal type:  root -l "AveragePlaquetteSU2.C(2.3)" 
//Contact Amit Kumar for questions (amitkr2410@gmail.com or kumar.amit@wayne.edu)

#include"iostream"
#include"fstream"
#include"TMath.h"
//#include"stdlib.h" //for srand and rand
#include"time.h"  //for time() function
#include"TCanvas.h"
#include"TGraph.h"
#include"TComplex.h"
#include"cmath"

using namespace std;
#define register      // Deprecated in C++11. Needed for "MersenneTwister.h"                                                                                 
#include"MersenneTwister.h" //for improved rand() function; Has a class "MTRand" and inside it has a function rand()( [0,1]);Creat ObjectMTRand to use


class LatticeCalculation
{
public:

  double PI; //Value of PI=3.1415926  
  int Nt; //Discretizing each dimension in N-points
  int Nx; //Discretizing each dimension in N-points 
  int Ny; //Discretizing each dimension in N-points 
  int Nz; //Discretizing each dimension in N-points 
  int LatticePoints; //Nt*Nx*Ny*Nz 
  int TotalLinks;   // 4 space-time dimension
  int NLinksX;      // Nt*Nx*Ny*Nz Links in X-direction
  int NLinksY;      // Nt*Nx*Ny*Nz Links in Y-direction
  int NLinksZ;      // Nt*Nx*Ny*Nz Links in Z-direction                  
  int NLinksT;      // Nt*Nx*Ny*Nz Links in T-direction
  int Run;          //  Number of iteration to get stable gauge field configuration
  int* LatticeSize;  //  LatticeSize[4]={Nt,Nx,Ny,Nz};
  double beta;

  MTRand ObjectMTRand;//Define a object to access rand() function defined inside class MTRand in file MersenneTwister.h

  double** LinksContent;   //LinksContent[TotalLinks][4];
  int** LinksCoordinate;   //LinksCoordinate[TotalLinks][4];
  int* LinksDirection;     //LinksDirection[TotalLinks];
  int CurrentLink, Iteration;
   
  TComplex ID[2][2],Pauli[3][2][2],CPauli[3][2][2],IYOTA,CIYOTA; //Complex Matrix
  
  TComplex SumProduct3Us[2][2];  //Complex matrix Sigma U2*U3*U4

  void DefineDynamicArrays();
  void LabelLinksCoordinate();
  Int_t OutSideBoundary(Int_t,Int_t ,Int_t, Int_t,Int_t);
  void UpdateLink(Int_t);
  Int_t SearchLinkNumber(Int_t, Int_t, Int_t, Int_t, Int_t);
  int Modulo(int,int);
  void DefineComplexMatrix();
  double TracePlaquette(int,int,int,int);
};

void LatticeCalculation::DefineDynamicArrays()
{

  LinksCoordinate = new int* [TotalLinks];
  LinksDirection = new int [TotalLinks];
  LinksContent = new double* [TotalLinks];

  LatticeSize = new int [4];

  for(int i=0; i<TotalLinks; i++)
    {
      LinksCoordinate[i] = new int [4];
      LinksContent[i] = new double [4];
    }

  LatticeSize[0]=Nt;
  LatticeSize[1]=Nx;
  LatticeSize[2]=Ny;
  LatticeSize[3]=Nz;
  //cout<<"Inside function _DefineArrayMatric__Lattice Size = "<<  LatticeSize[0]<<","<<LatticeSize[1]<<","<<LatticeSize[2]<<","<<LatticeSize[3]<<endl;
}


void LatticeCalculation::LabelLinksCoordinate()
{
  CurrentLink = 0;
  for(int a = 0; a < Nt ; a++)
    {
      for(Int_t b= 0; b < Nx ; b++)
        {
          for(Int_t c= 0; c < Ny ; c++)
            {
              for(Int_t d = 0; d< Nz ; d++)
                {
                  //z-direction                                                           
                  {
                    LinksCoordinate[CurrentLink][0]=a;
                    LinksCoordinate[CurrentLink][1]=b;
                    LinksCoordinate[CurrentLink][2]=c;
                    LinksCoordinate[CurrentLink][3]=d;
                    LinksDirection[CurrentLink]=3;
		    //std::cout<<LinksCoordinate[CurrentLink][0]<<LinksCoordinate[CurrentLink][1]<<LinksCoordinate[CurrentLink][2]<<LinksCoordinate[CurrentLink][3]<<endl;
                  }
		  
		  //y-dir
		  {
                    LinksCoordinate[CurrentLink + NLinksZ][0]=a;
                    LinksCoordinate[CurrentLink + NLinksZ][1]=b;
                    LinksCoordinate[CurrentLink + NLinksZ][2]=c;
                    LinksCoordinate[CurrentLink + NLinksZ][3]=d;
                    LinksDirection[CurrentLink + NLinksZ]=2;
		    // std::cout<<LinksCoordinate[CurrentLink+NLinks][0]<<LinksCoordinate[CurrentLink][1]<<LinksCoordinate[CurrentLink][2]<<LinksCoordinate[CurrentLink][3]<<endl;
                  }

                  //x-direction                                                          
                  {
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][0]=a;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][1]=b;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][2]=c;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][3]=d;
                    LinksDirection[CurrentLink + NLinksZ + NLinksY]=1;
                  }

                  //t-direction                                                           
                  {
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][0]=a;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][1]=b;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][2]=c;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][3]=d;
                    LinksDirection[CurrentLink + NLinksZ + NLinksY + NLinksX]=0;
                  }
		  CurrentLink = CurrentLink +1;
                }
            }
        }
    }

  CurrentLink =0;
}

int LatticeCalculation::OutSideBoundary(int t,int x, int y,int z, int d) // "i" is Link index
{
  Int_t State = 1; //Chosse zero if inside
  Int_t xt[4];
  xt[0]=t; xt[1]=x; xt[2]=y; xt[3]=z;
  //cout<<"txyz =" <<t<<x<<y<<z<<"\t"<<d<<endl;
  
  if( (t ==-1 ||x==-1 ||y==-1|| z==-1) || (t== Nt || x==Nx || y==Ny || z==Nz)   )    
    {
      //cout<<"Outside"<<endl;
      return 1;
    }
  else 
    {
      //cout<<"Inside"<<endl;
      return 0;
    }
  
}

void LatticeCalculation::UpdateLink(int  i)
{
  int StapleCoordinate[3][6][5]={0};
  int CoordinateLinkI[5]={0};
  int LinkNumberStaple[6][3]={0}; //6 loop; each with three links

  CurrentLink = i;
  //cout<<"\nCurrent Link"<<CurrentLink;cout<<"\n"<<"Current link Coordinate are txyz = ";
  for (Int_t m =0 ; m<4; m++)
    {
      CoordinateLinkI[m] = LinksCoordinate[CurrentLink][m];
  //cout<<CoordinateLinkI[m];
    }
  //cout<<"\n\n";
  CoordinateLinkI[4] = LinksDirection[CurrentLink];
  int d = LinksDirection[CurrentLink];

  //cout<<"Identifying Two-loop plane"<<endl;
  int p = -1; //p signifyies which two-loop plane we are in; p =0,1,2
  int count = -1; //
  for(int m = 0; m<4; m++)
    {
      
    	  if(m!=d)
	    {
	      p = p+1;
	      for(int ra=0 ; ra<6; ra++)
		{
		  for(int rb=0; rb<4 ; rb++ )
		    {
		      StapleCoordinate[p][ra][rb] = CoordinateLinkI[rb];
		    }
		}
	      
		  StapleCoordinate[p][0][4]= m;
		  StapleCoordinate[p][1][4]= d;
		  StapleCoordinate[p][2][4]= m;
		  StapleCoordinate[p][3][4]= m;
		  StapleCoordinate[p][4][4]= d;
		  StapleCoordinate[p][5][4]= m;
		  
		  StapleCoordinate[p][1][m]=StapleCoordinate[p][1][m] + 1;
                  StapleCoordinate[p][2][d]=StapleCoordinate[p][2][d] + 1;

		  StapleCoordinate[p][3][m]=StapleCoordinate[p][3][m] - 1;
                  StapleCoordinate[p][4][m]=StapleCoordinate[p][4][m] - 1;
		  StapleCoordinate[p][5][m]=StapleCoordinate[p][5][m] - 1;
		  StapleCoordinate[p][5][d]=StapleCoordinate[p][5][d] + 1;	      


		  
		  /*cout<<"Before PC"<<endl;
		  for(Int_t j=0; j<6; j++)
		   {for(Int_t k=0; k<5; k++)
		    {cout<<StapleCoordinate[p][j][k]<<"\t";}cout<<"\n"<<endl;}cout<<"\n"<<endl;
		  */
		  
		  for(Int_t g =0; g<6; g++)
		    {
		      Int_t Status = OutSideBoundary(StapleCoordinate[p][g][0],StapleCoordinate[p][g][1],StapleCoordinate[p][g][2],StapleCoordinate[p][g][3],StapleCoordinate[p][g][4]);
		      
		      if(Status == 1) //Status =1 means outside boundary
			{
			  
			  if(StapleCoordinate[p][g][m] == LatticeSize[m] )
			    {
			      StapleCoordinate[p][g][m] = StapleCoordinate[p][g][m] - (LatticeSize[m]);
			    }
			  
			  if(StapleCoordinate[p][g][d] == LatticeSize[d] )
                            {
                              StapleCoordinate[p][g][d] = StapleCoordinate[p][g][d] - (LatticeSize[d]);
                            }
			  if(StapleCoordinate[p][g][m] == -1 )
			    {
			      StapleCoordinate[p][g][m] = StapleCoordinate[p][g][m] + (LatticeSize[m]);
			    }
			  if(StapleCoordinate[p][g][d] == -1 )
                            {
                              StapleCoordinate[p][g][d] = StapleCoordinate[p][g][d] + (LatticeSize[d]);
                            }

			}
		    }
		  
		  /* cout<<"After PC"<<endl;
		  for(Int_t j=0; j<6; j++)
		  {for(Int_t k=0; k<5; k++)
		  {cout<<StapleCoordinate[p][j][k]<<"\t";}cout<<"\n"<<endl;}cout<<"\n"<<endl;
		  */
		  
		  
	    }
	  
    }  
  
  //cout<<"Done"<<endl;
  //Search Linknumber using t,x,y,z,d
  //cout<<"Search LinkNumber around link with txyzd ="<<CoordinateLinkI[0]<<CoordinateLinkI[1]<<CoordinateLinkI[2]<<CoordinateLinkI[3]<<CoordinateLinkI[4]<<endl;  
  count = -1;

  for(Int_t p = 0;  p <3; p++ )
    {
      for( Int_t ra = 0 ; ra< 6;ra++)
	{
	  if(ra == 0 || ra == 3)
	    {
	      count = count + 1;
	    }
	      Int_t i2 = ra%3;
	      //cout<<"Searching Link around "<<endl;
	      LinkNumberStaple[count][i2] = SearchLinkNumber( StapleCoordinate[p][ra][0], StapleCoordinate[p][ra][1], StapleCoordinate[p][ra][2], StapleCoordinate[p][ra][3], StapleCoordinate[p][ra][4] );
	      //cout<<"LinkNumberFound "<<LinkNumberStaple[count][i2]<<endl;
	      // cout<<"Content is a0 a1 a2 a3 "<<LinksContent[LinkNumberStaple[count][i2]][0]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][1]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][2]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][3]<<endl;
	}
    }
  
  //cout<<"Six Loops identified"<<endl;
  //Lets do matrix multiplication and summation ( summation U2*U3*U4), first initialize SumProduct3Us to zero
  //cout<<"Starting matrixmultiplication and summation"<<endl;

  SumProduct3Us[0][0] = TComplex(0.0,0.0);
  SumProduct3Us[1][0] = TComplex(0.0,0.0);
  SumProduct3Us[0][1] = TComplex(0.0,0.0);
  SumProduct3Us[1][1] = TComplex(0.0,0.0);


  TComplex A,B,C; 
  A=TComplex(0.0,0.0);
  B=TComplex(0.0,0.0);
  C=TComplex(0.0,0.0);

  for( Int_t p = 0; p<6; p++)
    {
      //cout<<"LinksContent ="<<LinksContent[LinkNumberStaple[p][2]]<<endl;
      for(Int_t i = 0; i<2; i++)
	{
	  for(Int_t j = 0; j<2; j++)
	    {
	      for(Int_t m = 0; m<2; m++)
		{
		  for(Int_t n = 0; n<2; n++)
		    {
		      if(p==0 || p==2 || p==4)
			{
			  
		       A = LinksContent[LinkNumberStaple[p][2]][0]*ID[i][m] + IYOTA*( LinksContent[LinkNumberStaple[p][2]][1]*Pauli[0][i][m] + LinksContent[LinkNumberStaple[p][2]][2]*Pauli[1][i][m] + LinksContent[LinkNumberStaple[p][2]][3]*Pauli[2][i][m]  );

		      B = LinksContent[LinkNumberStaple[p][1]][0]*ID[n][m] + CIYOTA*( LinksContent[LinkNumberStaple[p][1]][1]*CPauli[0][n][m] + LinksContent[LinkNumberStaple[p][1]][2]*CPauli[1][n][m] + LinksContent[LinkNumberStaple[p][1]][3]*CPauli[2][n][m] );

		      C = LinksContent[LinkNumberStaple[p][0]][0]*ID[j][n] + CIYOTA*( LinksContent[LinkNumberStaple[p][0]][1]*CPauli[0][j][n] + LinksContent[LinkNumberStaple[p][0]][2]*CPauli[1][j][n] + LinksContent[LinkNumberStaple[p][0]][3]*CPauli[2][j][n] );
			}
		      
		      if(p==1 || p==3 || p==5)
			{
			  A = LinksContent[LinkNumberStaple[p][2]][0]*ID[m][i] + CIYOTA*( LinksContent[LinkNumberStaple[p][2]][1]*CPauli[0][m][i] + LinksContent[LinkNumberStaple[p][2]][2]*CPauli[1][m][i] + LinksContent[LinkNumberStaple[p][2]][3]*CPauli[2][m][i] );

			  B = LinksContent[LinkNumberStaple[p][1]][0]*ID[n][m] + CIYOTA*( LinksContent[LinkNumberStaple[p][1]][1]*CPauli[0][n][m] + LinksContent[LinkNumberStaple[p][1]][2]*CPauli[1][n][m] + LinksContent[LinkNumberStaple[p][1]][3]*CPauli[2][n][m] );

			  C = LinksContent[LinkNumberStaple[p][0]][0]*ID[n][j] + IYOTA*( LinksContent[LinkNumberStaple[p][0]][1]*Pauli[0][n][j] + LinksContent[LinkNumberStaple[p][0]][2]*Pauli[1][n][j] + LinksContent[LinkNumberStaple[p][0]][3]*Pauli[2][n][j] );
			}
		      //cout<<"LinkNumber = "<<LinkNumberStaple[p][0]<<"\tContent = "<<LinksContent[LinkNumberStaple[p][0]][0]<<LinksContent[LinkNumberStaple[p][0]][1]<<LinksContent[LinkNumberStaple[p][0]][2]<<LinksContent[LinkNumberStaple[p][0]][3]<<"\n";
		      //cout<<"A = "<<A<<"\t B= "<<B<<"\t C="<<C<<"\t";
		      SumProduct3Us[i][j] = SumProduct3Us[i][j] + A*B*C ;
		      //cout<<"A = "<<A<<"\tB ="<<B<<"\tC="<<C<<endl;
		    } 
		}
	      //cout<<SumProduct3Us[i][j]<<"\t";
	    }
	  //cout<<"\n";
	}
      //cout<<"LinkNumber = "<<LinkNumberStaple[p][0]<<"\tContent = "<<LinksContent[LinkNumberStaple[p][0]][0]<<LinksContent[LinkNumberStaple[p][0]][1]<<LinksContent[LinkNumberStaple[p][0]][2]<<LinksContent[LinkNumberStaple[p][0]][3]<<"\n";
      //cout<<"LinkNumber = "<<LinkNumberStaple[p][1]<<"\tContent = "<<LinksContent[LinkNumberStaple[p][1]][0]<<LinksContent[LinkNumberStaple[p][1]][1]<<LinksContent[LinkNumberStaple[p][1]][2]<<LinksContent[LinkNumberStaple[p][1]][3]<<"\n";
      //cout<<"LinkNumber = "<<LinkNumberStaple[p][2]<<"\tContent = "<<LinksContent[LinkNumberStaple[p][2]][0]<<LinksContent[LinkNumberStaple[p][2]][1]<<LinksContent[LinkNumberStaple[p][2]][2]<<LinksContent[LinkNumberStaple[p][2]][3]<<"\n";


    }
  
  //cout<<"Finding determinant"<<endl;
	  //Determinant of SumProduct3Us matrix
	  TComplex D =  SumProduct3Us[0][0]*SumProduct3Us[1][1] - SumProduct3Us[0][1]*SumProduct3Us[1][0];
	  Double_t SqrtDet = TMath::Sqrt(TMath::Power( D.Re()*D.Re() + D.Im()*D.Im() , 0.5));
	  //if(CurrentLink==4101){cout<<"Determinant is "<<D<<endl;}
	  //cout<<"Sqrt of determinant = "<<SqrtDet<<endl;
	  //Create pdf and generate a0 a1 a2 a3 with apporpiate weight 
	  
	  Double_t K = SqrtDet, a0, a1, a2, a3;
	  
	  /////////////////////............Differential Sampling for a0.......////////////////////////
	  
	    Double_t a = 1.0;
	    Double_t Zmin, Zmax, Z, dpdz, dpdzmax=0.0, dz;
	    Double_t r1 , r2;
	    Int_t Success =0;
	    Zmin = TMath::Exp(-a*beta*K);	    
	    Zmax = TMath::Exp(a*beta*K);
	    
	    //dz = (Zmax-Zmin)/(100000.0-1.0);
	    //for(Int_t i=0; i< 100000; i++)
	    //  {Z = Zmin + (i*dz);dpdz = sqrt(1 - (pow(beta*K,-2)*pow(log(Z),2)) );
	    //if( dpdzmax < dpdz ){dpdzmax = dpdz;}}
	    dpdzmax =1.0;
	    Int_t SearchCount = 0;//cout<<"Search for a0 started "<<endl;
	    while(Success==0)
	      {		    
		Double_t Dummy1 = ObjectMTRand.rand();
		Double_t Dummy2 = ObjectMTRand.rand();
		r1 = Dummy1;r2 = Dummy2;
		Z = Zmin + (Zmax - Zmin)*r1; 
		dpdz = TMath::Sqrt(1 - (TMath::Power(beta*K,-2)*TMath::Power(log(Z),2)) );
		
		SearchCount = SearchCount + 1;
		//cout<<"Search count ="<<SearchCount<<endl;
		if( dpdz > (dpdzmax*r2) )
		  {
		    Success =1;
		  }
		else
		  {Success =0;}
	      }

	    a0 = log(Z)/(beta*K);
	  
	    /////////////.......End of differential sampling algorithm for a0...///////////

                
            /////////////////////Using Integral sampling method for a0.......///////////////
	    /*
	  
            Double_t x0 =-1.0, fval1=0.,fval2=0.,fmax=0., maxarea = 0., area = 0.;
            Int_t n =1000;
            Double_t dx0 = 2.0/(n-1);
            double *normalizedarea = new double[n];normalizedarea[0]=0.0;
            double *x = new double[n];

            x0 = -1.0;
            fval1 = 0.0;
            area = 0.0;
            maxarea = 0.0;

            for(Int_t i = 0 ; i< n-1; i++)
              {
                x[i] = x0 + i*dx0;
                fval1 = TMath::Sqrt(1- TMath::Power(x0+i*dx0,2) )*TMath::Exp(beta*K*(x0+i*dx0));
                fval2 = TMath::Sqrt(1- TMath::Power(x0+(i+1)*dx0,2) )*TMath::Exp(beta*K*(x0+(i+1)*dx0));
                maxarea = maxarea + dx0*0.5*(fval1 + fval2);
              }

            x[n-1]=1.0;

            for(Int_t i = 0 ; i< n-1; i++)
              {
                fval1 = TMath::Sqrt(1- TMath::Power(x0+i*dx0,2) )*TMath::Exp(beta*K*(x0+i*dx0));
                fval2 = TMath::Sqrt(1- TMath::Power(x0+(i+1)*dx0,2) )*TMath::Exp(beta*K*(x0+(i+1)*dx0));
                area = area + dx0*0.5*(fval1 + fval2);
                normalizedarea[i+1] = area/maxarea;
              }
	    
	    TGraph *GPdf = new TGraph(n,normalizedarea,x);
	    Double_t Dummy1 = ObjectMTRand.rand();
            Double_t Area = Dummy1;
	    a0 = GPdf->Eval(Area); 
            */
	    ////////..........End of Integral method of sampling a0..........///////////////////
	  
	    ///////////////A different Algorithm for a0 generation.../////////
	  /*
	  Double_t lambda0 =0.0, fval1=0.,fval2=0.,fmax=0., maxarea = 0., area = 0.;
	  Int_t n =1000;
	  Double_t dlambda0 = 1.0/(n-1);
	  double *normalizedarea = new double[n];normalizedarea[0]=0.0;
	  
	  double *lambda = new double[n];
	  
                                                                                                                     
	  lambda0 = 0.0;
	  fval1 = 0.0;
	  area = 0.0; 
	  maxarea = 0.0; 
	  
	  for(Int_t i = 0 ; i< n-1; i++)
	    {                                                                                                      
	      lambda[i] = lambda0 + i*dlambda0;
	      fval1 = TMath::Power((lambda0+ i*dlambda0),2)*TMath::Exp(-2*K*beta*TMath::Power((lambda0+ i*dlambda0),2));
	      fval2 = TMath::Power((lambda0+ (i+1)*dlambda0),2)*TMath::Exp(-2*K*beta*TMath::Power((lambda0+ (i+1)*dlambda0),2));
	      maxarea = maxarea + dlambda0*0.5*(fval1 + fval2);                                                         
	    }                           
	  lambda[n-1]=1.0;                                                                                                                                                                                                                   
	  for(Int_t i = 0 ; i< n-1; i++)                                                                           
	    {                                                                                                      
	      fval1 = TMath::Power((lambda0+ i*dlambda0),2)*TMath::Exp(-2*K*beta*TMath::Power((lambda0+ i*dlambda0),2));
              fval2 = TMath::Power((lambda0+ (i+1)*dlambda0),2)*TMath::Exp(-2*K*beta*TMath::Power((lambda0+ (i+1)*dlambda0),2));
	      area = area + dlambda0*0.5*(fval1 + fval2);                                                               
	      normalizedarea[i+1] = area/maxarea;                                                                  
	    }                                                                                                      
                                                                                                                     
	  TGraph *GPdf = new TGraph(n,normalizedarea,lambda);                                        
	  Double_t lambdaTrial = 0.0; Int_t Success = 0; 	  

	  while(Success==0)
	    {
	      Double_t Dummy1 = ObjectMTRand.rand();
	      Double_t Area = Dummy1;	      
	      lambdaTrial = GPdf->Eval(Area);                      
	      Double_t Dummy2 = ObjectMTRand.rand();
	      if(Dummy2 <= TMath::Sqrt(1- TMath::Power(lambdaTrial,2)))
		{
		  Success =1;
		}
	      else
		{Success =0;}
	    }
	  
	  a0 = 1.0 - 2.0*TMath::Power(lambdaTrial,2);
          */
	  //////////////...End of different algorithm..............///////////

	  ///////////..........Algorithm given in Book "QCD Lattice"........//////////
	  /*
	  Double_t r1,r2,r3,r; Double_t lambdaTrialS=0.0;
	  Int_t Success=0;

	  while(Success==0)
	    {
	      r1 = 1.0 - ObjectMTRand.rand();
	      r2 = 1.0 - ObjectMTRand.rand();
	      r3 = 1.0 - ObjectMTRand.rand();
	      
	      lambdaTrialS = - ( ( log(r1) +( TMath::Power(cos(2*PI*r2),2)*log(r3)) )/(2*K*beta) ) ;
	      r = ObjectMTRand.rand();
	      if(r <= TMath::Sqrt(1- lambdaTrialS))
		{
		  Success = 1;   a0= (1.0 - 2.0*lambdaTrialS); 
		}
	      else{ Success = 0;}
	    }
	  */	  
	  ///////////...End of algorithm given in Book "QCD Lattice".........////////////


	  // {cout<<"\n random number a0 = "<<a0<<endl;}
	    
	    Double_t Dummy = (ObjectMTRand.rand())*2.0;
	    Double_t Theta =  TMath::ACos( Dummy - 1.0);
	    Dummy = ObjectMTRand.rand();                                                                                                                   
           
	    Double_t Phi = Dummy*2.0*PI;
	    Double_t aVec = TMath::Sqrt(1.0- a0*a0 );                                                                                                                   
	    a1 = aVec*cos(Theta);
	    a2 = aVec*sin(Theta)*cos(Phi);
	    a3 = aVec*sin(Theta)*sin(Phi);
	    //cout<<"Theta =\t "<<Theta*180.0/PI;
	    //cout<<"\t Modulus of 4-vector a = "<<(a0*a0 + a1*a1 + a2*a2 + a3*a3)<<endl;
          
	 

	  
	  TComplex UCurrentLink[2][2]; 

	  //cout<<"\n\nCreated a0 a1 a2 a3 \t "<<a0<<"\t"<<a1<<"\t"<<a2<<"\t"<<a3<<endl;
	  UCurrentLink[0][0] = TComplex(0.0,0.0);
          UCurrentLink[1][0] = TComplex(0.0,0.0);
          UCurrentLink[0][1] = TComplex(0.0,0.0);
          UCurrentLink[1][1] = TComplex(0.0,0.0);
	  
	  for(Int_t i=0; i<2; i++)
            {
	      for(Int_t j =0; j<2; j++)
		{
		  for(Int_t m=0; m<2; m++)
		    {
		      UCurrentLink[i][j] = UCurrentLink[i][j] +( (a0*ID[i][m] + IYOTA*(a1*Pauli[0][i][m] + a2*Pauli[1][i][m] + a3*Pauli[2][i][m]) )*(TComplex::Conjugate(SumProduct3Us[j][m]))/K);
		      
		    }
		  
		}   //if(CurrentLink ==5000){cout<<"\n";}
	    }
	  
	  //cout<<(SumProduct3Us[0][0]*SumProduct3Us[1][1] - SumProduct3Us[0][1]*SumProduct3Us[1][0])/(K*K)<<endl;
	      LinksContent[CurrentLink][0] = UCurrentLink[0][0].Re();
	      LinksContent[CurrentLink][3] = UCurrentLink[0][0].Im();
	      LinksContent[CurrentLink][1] = UCurrentLink[0][1].Im();
	      LinksContent[CurrentLink][2] = UCurrentLink[0][1].Re();
              
	      //cout<<"Modulus of 4-vector after = "<<pow(UCurrentLink[0][0].Re(),2)+pow(UCurrentLink[0][0].Im(),2) + pow(UCurrentLink[0][1].Re(),2) + pow(UCurrentLink[0][1].Im(),2)<<endl;
	      
              
	      
}


int LatticeCalculation::SearchLinkNumber(int t, int x, int y, int z, int d)
{
  //cout<<"\t\t ..... Search Function started"<<endl;
  //cout<<"\t \t .... txyzd = "<<t<<","<<x<<","<<y<<","<<z<<","<<d<<endl;
  Int_t found =0; //1 means found
  Int_t Link = 0;
  Int_t ra = 0;
    
  while(ra< TotalLinks && found == 0)
    {
      if(d == LinksDirection[ra])
	{
	  if(t == LinksCoordinate[ra][0] && x == LinksCoordinate[ra][1] && y == LinksCoordinate[ra][2] && z == LinksCoordinate[ra][3])
	    {
	      found = 1;
	      Link = ra;
	      
	    }
	  else {ra = ra+1;}
	}
      else {ra = ra + 1;}
    }
  
  if(found ==1)
    {
      //cout<<"\t\t .... End using searchFunction with Link number = "<<Link<<endl;
      return Link;
    }
  else{
    //cout<<"\t Alert Link not found"<<endl;
    return -1;}

  //cout<<"\t\t .... End using searchFunction with Link number = "<<Link<<endl;
} 

double LatticeCalculation::TracePlaquette(int n1 ,int n2, int n3, int n4)
{
  //cout<<"n1 n2 n3 n4"<<n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<n4<<endl;
  TComplex U1,U2,U3,U4,TraceU;

  /*
  cout<<"Coordinate of Link n1 = "<<LinksCoordinate[n1][0]<<LinksCoordinate[n1][1]<<LinksCoordinate[n1][2]<<LinksCoordinate[n1][3]<<LinksDirection[n1]<<"\n";
  cout<<"LinksContent = "<<LinksContent[n1][0]<<"\t"<<LinksContent[n1][1]<<"\t"<<LinksContent[n1][2]<<"\t"<<LinksContent[n1][3]<<"\t Modulus = "<<pow(LinksContent[n1][0],2) + pow(LinksContent[n1][1],2) + pow(LinksContent[n1][2],2) + pow(LinksContent[n1][3],2)<<endl;

  cout<<"Coordinate of Link n2 = "<<LinksCoordinate[n2][0]<<LinksCoordinate[n2][1]<<LinksCoordinate[n2][2]<<LinksCoordinate[n2][3]<<LinksDirection[n2]<<"\n";
  cout<<"LinksContent = "<<LinksContent[n2][0]<<"\t"<<LinksContent[n2][1]<<"\t"<<LinksContent[n2][2]<<"\t"<<LinksContent[n2][3]<<"\t Modulus = "<<pow(LinksContent[n2][0],2) + pow(LinksContent[n2][1],2) + pow(LinksContent[n2][2],2) + pow(LinksContent[n2][3],2)<<endl;

  cout<<"Coordinate of Link n1 = "<<LinksCoordinate[n3][0]<<LinksCoordinate[n3][1]<<LinksCoordinate[n3][2]<<LinksCoordinate[n3][3]<<LinksDirection[n3]<<"\n";
  cout<<"LinksContent = "<<LinksContent[n3][0]<<"\t"<<LinksContent[n3][1]<<"\t"<<LinksContent[n3][2]<<"\t"<<LinksContent[n3][3]<<"\t Modulus = "<<pow(LinksContent[n3][0],2) + pow(LinksContent[n3][1],2) + pow(LinksContent[n3][2],2) + pow(LinksContent[n3][3],2)<<endl;

  cout<<"Coordinate of Link n1 = "<<LinksCoordinate[n4][0]<<LinksCoordinate[n4][1]<<LinksCoordinate[n4][2]<<LinksCoordinate[n4][3]<<LinksDirection[n4]<<"\n";
  cout<<"LinksContent = "<<LinksContent[n4][0]<<"\t"<<LinksContent[n4][1]<<"\t"<<LinksContent[n4][2]<<"\t"<<LinksContent[n4][3]<<"\t Modulus = "<<pow(LinksContent[n4][0],2) + pow(LinksContent[n4][1],2) + pow(LinksContent[n4][2],2) + pow(LinksContent[n4][3],2)<<endl;
  */

	  U1=TComplex(0.0,0.0); U2=TComplex(0.0,0.0);
	  U3=TComplex(0.0,0.0); U4=TComplex(0.0,0.0);TraceU=TComplex(0.0,0.0);
	
	  for(Int_t i=0; i<2; i++)
	    {
	      for(Int_t m=0; m<2; m++)
		{
		  for(Int_t n=0; n<2; n++)
		    {
		      for(Int_t p=0; p<2; p++)
			{
			  U1 = LinksContent[n1][0]*ID[i][m] + IYOTA*(LinksContent[n1][1]*Pauli[0][i][m] + LinksContent[n1][2]*Pauli[1][i][m] + LinksContent[n1][3]*Pauli[2][i][m]);

			  U2 = LinksContent[n2][0]*ID[m][n] + IYOTA*(LinksContent[n2][1]*Pauli[0][m][n] + LinksContent[n2][2]*Pauli[1][m][n] + LinksContent[n2][3]*Pauli[2][m][n]);

			  U3 = LinksContent[n3][0]*ID[p][n] + CIYOTA*(LinksContent[n3][1]*CPauli[0][p][n] + LinksContent[n3][2]*CPauli[1][p][n] + LinksContent[n3][3]*CPauli[2][p][n]);

			  U4 = LinksContent[n4][0]*ID[i][p] + CIYOTA*(LinksContent[n4][1]*CPauli[0][i][p] + LinksContent[n4][2]*CPauli[1][i][p] + LinksContent[n4][3]*CPauli[2][i][p]);

			  TraceU = TraceU + U1*U2*U3*U4;
			}
		    }
		}
	    }

  //cout<<"Trace U = "<<TraceU<<endl;
  //if(TraceU.Re() < 0.0){TraceU=TComplex(0.0,0.0);}
  return (TraceU.Re());
}

void LatticeCalculation::DefineComplexMatrix()
{

//Define Identity and Pauli Matrices
   
  ID[0][0]= TComplex(1.0,0.0);  ID[0][1]= TComplex(0.0,0.0);
  ID[1][0]= TComplex(0.0,0.0);  ID[1][1]= TComplex(1.0,0.0);
   
  Pauli[0][0][0] = TComplex(0.0,0.0); Pauli[0][0][1] = TComplex(1.0,0.0);
  Pauli[0][1][0] = TComplex(1.0,0.0); Pauli[0][1][1] = TComplex(0.0,0.0);
  
  Pauli[1][0][0] = TComplex(0.0,0.0); Pauli[1][0][1] = TComplex(0.0,-1.0);
  Pauli[1][1][0] = TComplex(0.0,1.0); Pauli[1][1][1] = TComplex(0.0,0.0);

  Pauli[2][0][0] = TComplex(1.0,0.0); Pauli[2][0][1] = TComplex(0.0,0.0);
  Pauli[2][1][0] = TComplex(0.0,0.0); Pauli[2][1][1] = TComplex(-1.0,0.0);
  
  IYOTA = TComplex(0.0,1.0);

//Define Conjugate Pauli matrices and IYOta

  CPauli[0][0][0] = TComplex(0.0,0.0); CPauli[0][0][1] = TComplex(1.0,0.0);
  CPauli[0][1][0] = TComplex(1.0,0.0); CPauli[0][1][1] = TComplex(0.0,0.0);

  CPauli[1][0][0] = TComplex(0.0,0.0); CPauli[1][0][1] = TComplex(0.0,1.0);
  CPauli[1][1][0] = TComplex(0.0,-1.0);CPauli[1][1][1] = TComplex(0.0,0.0);

  CPauli[2][0][0] = TComplex(1.0,0.0); CPauli[2][0][1] = TComplex(0.0,0.0);
  CPauli[2][1][0] = TComplex(0.0,0.0); CPauli[2][1][1] = TComplex(-1.0,0.0);

  CIYOTA = TComplex(0.0,-1.0);

}

void AveragePlaquetteSU2(double Beta)
{
  int StartTime = time(NULL);

  //const double PI = 3.1415926; //Value of pi
  //const int NIteration = 10;
  //srand(time(NULL));
  LatticeCalculation LC;

  LC.PI = 3.1415926;
  LC.Nt = 8; //Discretizing each dimension in N-points
  LC.Nx = 8;
  LC.Ny = 8;
  LC.Nz = 8;
  LC.LatticePoints = LC.Nt*LC.Nx*LC.Ny*LC.Nz; //
  LC.TotalLinks = 4*LC.Nt*LC.Nx*LC.Ny*LC.Nz; // 4 space-time dimension
  LC.NLinksX = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in X-direction
  LC.NLinksY = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Y-direction 
  LC.NLinksZ = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Z-direction           
  LC.NLinksT = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in T-direction           
  LC.Run = 10;                //Number of iteration to get stable gauge field configuration  


  LC.DefineDynamicArrays();
  LC.DefineComplexMatrix();
  //.............End of Input Parameter.....Lattice

  LC.beta = Beta;
  //LC.beta = 0.2 + Beta*0.2;
  cout<<"Program started"<<endl;
  cout<<"Dimension of Lattice used is  \t" <<LC.Nt<<"*"<<LC.Nx<<"*"<<LC.Ny<<"*"<<LC.Nz<<", So Total Links are "<<LC.TotalLinks<<endl;
  cout<<"Total Plaquette should be 6*Nt*Nx*Ny*Nz "<<6*LC.Nt*LC.Nx*LC.Ny*LC.Nz<<endl;
  
  //ofstream f1;
  //f1.open("Coordinate.txt",ios::out);
  //f1.open("AveragePlaquetteForBeta3_2.txt",ios::out);
  //f1.open("Plaquette_wrt_beta_N6_I30.txt",ios::app);
  //f1.open("Plaquette_wrt_beta_N3_12_12_12_I30.txt",ios::app);


//Labeling each link with its 4 space-time coordinate and a direction
  LC.LabelLinksCoordinate();     
    
//Setting all link variables to Unity
  for(Int_t i = 0; i<LC.TotalLinks ; i++)
    { 
      LC.LinksContent[i][0]=1.0;      //a0
      LC.LinksContent[i][1]=0.0;      //a1
      LC.LinksContent[i][2]=0.0;      //a2
      LC.LinksContent[i][3]=0.0;      //a3
    }
 

  int TotalPlaquette = 6*LC.Nt*LC.Nx*LC.Ny*LC.Nz;

  double *ValueEachPlaquette = new double[TotalPlaquette];  
  double *AveragePlaquette = new double[LC.Run];
  double *AveragePlaquetteError = new double[LC.Run];
  double SumDeviationSquare; 
  double *SumPlaquette = new double[LC.Run];
  double *RUN = new double[LC.Run];
  double *RUNError = new double[LC.Run];
  for(Int_t in=0; in<LC.Run; in++)
    { AveragePlaquette[in] = 0.0;SumPlaquette[in]=0.0;RUN[in]=in;AveragePlaquetteError[in]=0.0;RUNError[in]=0.0;}
  
for(Int_t Iter=0; Iter<LC.Run; Iter++)
  {

    LC.Iteration = Iter;    
    cout<<"Performing Iteration "<<Iter<<endl;  
    //Touching each internal link with heat bath 
    for(Int_t ic = 0 ; ic< LC.TotalLinks ; ic++)
      {
	//cout<<"Processing Link"<<ic<<endl;
	LC.UpdateLink(ic); //Apply heat bath algorithm to each internal link and apply periodic boundary condition
	//cout<<"Content of Link "<<LC.LinksContent[ic][0]<<"\t"<<LC.LinksContent[ic][1]<<"\t"<<LC.LinksContent[ic][2]<<"\t"<<LC.LinksContent[ic][3]<<endl;   
      }
   
  
   
  //Evaluating plaquette
   Int_t xt[4] ={0,0,0,0};
   Int_t p1=0,p2=0;
   Int_t PlaquetteCount=0;
   SumDeviationSquare = 0.0;
   PlaquetteCount =0;
   
   for(Int_t d1 = 0; d1 <3; d1++)
     {
       for(Int_t d2 =d1+1; d2< 4; d2++)
	 {
	   // cout<<"("<<d1<<","<<d2<<")"<<endl;
	   if(d1==0 && d2==1){p1=2;p2=3;}
	   if(d1==0 && d2==2){p1=1;p2=3;}
	   if(d1==0 && d2==3){p1=1;p2=2;}
	   if(d1==1 && d2==2){p1=0;p2=3;}
	   if(d1==1 && d2==3){p1=0;p2=2;}
	   if(d1==2 && d2==3){p1=0;p2=1;}
	      
	  for(Int_t i=0; i<LC.LatticeSize[p1]; i++)  //axis perpendicular to V-axis and H-axis
	    {
	      for(Int_t j=0; j<LC.LatticeSize[p2]; j++)//axis perpendicular to V-axis and H-axis
		{
		  for(Int_t k=0; k<LC.LatticeSize[d1]; k++)  //Vertical-axis
		    {
		      for(Int_t m=0; m<LC.LatticeSize[d2]; m++) //Horizontal axis
			{
			  xt[d1]=k;
			  xt[d2]=m;
			  xt[p1]=i;
			  xt[p2]=j;
		
			 
			      Int_t Link4 = LC.SearchLinkNumber(xt[0],xt[1],xt[2],xt[3],d1);
			      Int_t Link1 = LC.SearchLinkNumber(xt[0],xt[1],xt[2],xt[3],d2);


			      xt[d2]=xt[d2]+1;
			      
			      if(xt[d2] == LC.LatticeSize[d2])
				{
				  xt[d2]=0;
				}			      				  
			      Int_t Link2 = LC.SearchLinkNumber(xt[0],xt[1],xt[2],xt[3],d1);
			      xt[d2] = m;
			      
			      xt[d1] = xt[d1]+1;
			      if( xt[d1] == LC.LatticeSize[d1])
                                {xt[d1]=0;
                                }
			      Int_t Link3 = LC.SearchLinkNumber(xt[0],xt[1],xt[2],xt[3],d2);
			      
			      PlaquetteCount = PlaquetteCount + 1;
			      ValueEachPlaquette[PlaquetteCount-1] = 1.0 - (0.5*LC.TracePlaquette(Link1,Link2,Link3,Link4));
			      SumPlaquette[Iter] = SumPlaquette[Iter] + ValueEachPlaquette[PlaquetteCount-1];
			      //cout<<"Value of current plaquette is "<<ValueEachPlaquette[PlaquetteCount-1]<<endl;
			      
			}
		      
			   
			  
		    }
		}
	    }
	}
    }
  
   //cout<<"Total number of Plaquette = "<<PlaquetteCount<<endl;
   //cout<<"SumPlaquette "<<SumPlaquette[Iter]<<endl;
    AveragePlaquette[Iter] = SumPlaquette[Iter]/PlaquetteCount;
    
    //Calculating Error in AveragePlaquette
    for(int i=0; i<PlaquetteCount; i++)
    {
      SumDeviationSquare = SumDeviationSquare + TMath::Power(AveragePlaquette[Iter] - ValueEachPlaquette[i] ,2.0);
    }
    
    AveragePlaquetteError[Iter] = TMath::Sqrt(SumDeviationSquare/(PlaquetteCount-1))/TMath::Sqrt(PlaquetteCount); 
    
    
    //cout<<Iter<<"\t"<<AveragePlaquette[Iter]<<endl;
    
    }
    
    for(Int_t id= 0; id<LC.Run; id++)
    {
    cout<<"AveragePlaquette with error at run "<<id<<"\t is \t"<<AveragePlaquette[id]<<"\t +-\t"<<AveragePlaquetteError[id]<<"\t"<<AveragePlaquetteError[id]*100.0/AveragePlaquette[id]<<"%"<<endl;
    //f1<<id<<"\t"<<AveragePlaquette[id]<<endl;
   }
   
   cout<<"Ending program"<<endl;
   
   
   //Plotting the AveragePlaquette vs Iteration
   
   TCanvas *c1 = new TCanvas("C1","c1");    
   TGraphErrors *GAvg = new TGraphErrors(LC.Run,RUN,AveragePlaquette,RUNError,AveragePlaquetteError);
   
   GAvg->Draw("AP");
   GAvg->SetMarkerStyle(6);
   GAvg->SetMarkerSize(1.0);
   GAvg->GetXaxis()->SetRangeUser(-1.,LC.Run);
   GAvg->GetYaxis()->SetRangeUser(0.2,1.0);
   GAvg->SetTitle("Average Plaquette (Lattice Size=3^4, #beta=2.3, Identity as I.C.) ");
   GAvg->GetYaxis()->SetTitle("Average Plaquette");
   GAvg->GetXaxis()->SetTitle("Iteration");
   GAvg->GetYaxis()->CenterTitle();
   GAvg->GetXaxis()->CenterTitle();
   //c1->SaveAs("Plaquette_N3_I30_Beta2_3.png");
   
   //f1<<LC.beta<<"\t"<<AveragePlaquette[LC.Run-1]<<"\t"<<AveragePlaquetteError[LC.Run-1]<<endl;
   //f1.close();
   
}
