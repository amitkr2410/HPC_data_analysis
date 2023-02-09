//This program calculates the Field-strength-Field-Strength correlation for SU(2) gauge field configuration for Nt*Nx*Ny*Nz Lattice with Periodic BC's
//  Following is lattice structure for Lattice size 3*4  (Notice the missing link at boundary; they are to incorporate the periodic BC )
//
//     |_|_|_|_
//     |_|_|_|_
//     |_|_|_|_
//
//We use Cpp function "rand()" to generate random number
//We use 5 variables to store the coordinate of the link (first four for t,x,y,z and last-one for its direction; We use 0 to point in t-dir, 1 (x-dir).. ) 
//We use Heat-bath algorithm to update each link and generate gauge field configuration
//Statistical Uncertainty are considered in estimated in error
//Program is written in C++ and uses ROOT-CERN library
//To compile this code for Beta=4/g^2 = 2.3,  in terminal type:  root -l "AverageCorrelatorF3iF3iMinusF4iF4iSU2.C(2.3)" 
//Contact Amit Kumar for questions (amitkr2410@gmail.com or kumar.amit@wayne.edu)

#include"iostream"
#include"fstream"
#include"TMath.h"
//#include"stdlib.h" //for srand and rand
#include"time.h"  //for time() function
#include"TCanvas.h"
#include"TGraph.h"
#include"TComplex.h"
#include"math.h"

using namespace std;
#define register      // Deprecated in C++11. Needed for "MersenneTwister.h"                                                                                 
#include"MersenneTwister.h" //for improved rand() function; Has a class "MTRand" and inside it has a function rand()( [0,1]);Creat ObjectMTRand to use


class LatticeCalculation
{
public:

  double PI;//Value of PI=3.1415926  
  int Nt; //Discretizing each dimension in N-points
  int Nx;
  int Ny;
  int Nz;
  int LatticePoints; // 
  int TotalLinks; // 4*Nt*Nx*Ny*Nz space-time dimension
  int NLinksX;      // Nt*Nx*Ny*Nz Links in X-direction
  int NLinksY;      // Nt*Nx*Ny*Nz Links in Y-direction
  int NLinksZ;      // Nt*Nx*Ny*Nz Links in Z-direction                 
  int NLinksT;      // Nt*Nx*Ny*Nz Links in T-direction
  int Run;                //Number of iteration to get stable gauge field configuration
  int TotalConfiguration; //Number of gauge field configuartion over which average of (F3iF3i-F4iF4i) is computed 
  int TotalTemperatureRun; //Number of different values of Temperature to see the dependence of correlation function  
  int* LatticeSize;//LatticeSize[4] = {Nt,Nx,Ny,Nz}
  double lambdaLattice; //(in MeV) A parameter that goes in LatticeSpacing formula using renormalization group equation    

  double beta, temperature, LatticeSpacing, BareCouplingConstant;

  MTRand ObjectMTRand;//MTRand Class is defined in "MTwister.h" which has rand() function defined
  
  int CurrentLink, Iteration;

  int** LinksCoordinate; //To define LinksCoordinate[TotalLinks][4];
  int* LinksDirection;  //To define LinksDirection[TotalLinks];
  double** LinksContent;//To define LinksContent[TotalLinks][4];

  int***** LinksNumber; //LinkNumber[t,x,y,z,d]

  TComplex ID[2][2],Pauli[3][2][2],CPauli[3][2][2],IYOTA,CIYOTA; //Complex Matrix
  
  TComplex SumProduct3Us[2][2];  //Complex matrix Sigma U2*U3*U4

  void DefineDynamicArrays();
  void LabelLinksCoordinate();
  Int_t OutSideBoundary(Int_t,Int_t ,Int_t, Int_t,Int_t);
  void UpdateLink(Int_t);
  //Int_t SearchLinkNumber(Int_t, Int_t, Int_t, Int_t, Int_t);
  void DefineComplexMatrix();

  double TracePlaquette(Int_t,Int_t,Int_t,Int_t);
  double TraceFieldStrengthSquare(int ,int ,int , int );

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

  LinksNumber = new int**** [Nt];
  for(int i=0; i<Nt; i++){LinksNumber[i] = new int*** [Nx];}

  for(int i=0; i<Nt; i++){for(int j=0; j<Nx; j++){ LinksNumber[i][j] = new int** [Ny];}}

  for(int i=0; i<Nt; i++) {for(int j=0; j<Nx; j++) {for(int k=0; k<Ny; k++) { LinksNumber[i][j][k] = new int* [Nz];}}}

  for(int i=0; i<Nt; i++) {for(int j=0; j<Nx; j++) {for(int k=0; k<Ny; k++) {for(int l=0; l<Nz; l++) { LinksNumber[i][j][k][l] = new int [4];}}}}



}



void LatticeCalculation::LabelLinksCoordinate()
{
  CurrentLink = 0;//cout<<"Value of total links = "<<TotalLinks<<"\t"<<LinksCoordinate[1][0]<<endl;
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

		    LinksNumber[a][b][c][d][3] = CurrentLink;
		    //std::cout<<LinksCoordinate[CurrentLink][0]<<LinksCoordinate[CurrentLink][1]<<LinksCoordinate[CurrentLink][2]<<LinksCoordinate[CurrentLink][3]<<endl;
                  }
		  
		  //y-dir
		  {
                    LinksCoordinate[CurrentLink + NLinksZ][0]=a;
                    LinksCoordinate[CurrentLink + NLinksZ][1]=b;
                    LinksCoordinate[CurrentLink + NLinksZ][2]=c;
                    LinksCoordinate[CurrentLink + NLinksZ][3]=d;
                    LinksDirection[CurrentLink + NLinksZ]=2;
		    
		    LinksNumber[a][b][c][d][2] = CurrentLink +  NLinksZ ;
		    //std::cout<<LinksCoordinate[CurrentLink+NLinks][0]<<LinksCoordinate[CurrentLink][1]<<LinksCoordinate[CurrentLink][2]<<LinksCoordinate[CurrentLink][3]<<endl;
                  }

                  //x-direction                                                          
                  {
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][0]=a;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][1]=b;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][2]=c;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY][3]=d;
                    LinksDirection[CurrentLink + NLinksZ + NLinksY]=1;

		    LinksNumber[a][b][c][d][1] = CurrentLink +  NLinksZ + NLinksY;
                  }

                  //t-direction                                                           
                  {
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][0]=a;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][1]=b;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][2]=c;
                    LinksCoordinate[CurrentLink + NLinksZ + NLinksY + NLinksX][3]=d;
                    LinksDirection[CurrentLink + NLinksZ + NLinksY + NLinksX]=0;

		    LinksNumber[a][b][c][d][0] = CurrentLink +  NLinksZ + NLinksY + NLinksX;
                  }
		  CurrentLink = CurrentLink +1;
		  //cout<<"CurrentLink = "<<CurrentLink<<endl;
                }
            }
        }
    }

  for(int i=0; i<TotalLinks; i++)
    {
      LinksContent[i][0]=0.0;
      LinksContent[i][1]=0.0;
      LinksContent[i][2]=0.0;
      LinksContent[i][3]=0.0;

    }

  CurrentLink =0;
}


int LatticeCalculation::OutSideBoundary(int t,int x, int y,int z, int d) // "i" is Link index
{
  Int_t State = 1; //Chosse zero if inside
  Int_t xt[4];
  xt[0]=t; xt[1]=x; xt[2]=y; xt[3]=z;
  //cout<<"txyz =" <<t<<x<<y<<z<<"\t"<<d<<endl;
  

  if( (t ==-1 ||x==-1 ||y==-1|| z==-1) || (t== Nt || x== Nx || y== Ny || z== Nz)   )    
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
  //  cout<<CoordinateLinkI[m];
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
			  //cout<<"Link is outside "<<endl;
			  //cout<<"Lattice Size = "<<LatticeSize[0]<<"\t"<<LatticeSize[1]<<"\t"<<LatticeSize[2]<<"\t"<<LatticeSize[3]<<endl;
			  
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
		  /*
		  cout<<"After PC"<<endl;
		  for(Int_t j=0; j<6; j++)
		  {for(Int_t k=0; k<5; k++)
		  {cout<<StapleCoordinate[p][j][k]<<"\t";}cout<<"\n"<<endl;}cout<<"\n"<<endl;
		  */
		  
		  
	    }
	  
    }  
    
  //cout<<"Done"<<endl;
  //Search Linknumber using t,x,y,z,d
  //cout<<"Search LinkNumber around link with txyzd =";
  //cout<<CoordinateLinkI[0]<<CoordinateLinkI[1]<<CoordinateLinkI[2]<<CoordinateLinkI[3]<<CoordinateLinkI[4]<<endl;  
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
	      LinkNumberStaple[count][i2] = LinksNumber[StapleCoordinate[p][ra][0]][ StapleCoordinate[p][ra][1]][ StapleCoordinate[p][ra][2]][ StapleCoordinate[p][ra][3]][StapleCoordinate[p][ra][4] ];
	      //cout<<"LinkNumberFound "<<LinkNumberStaple[count][i2]<<endl;
	      //cout<<"Content is a0 a1 a2 a3 "<<LinksContent[LinkNumberStaple[count][i2]][0]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][1]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][2]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][3]<<endl;
	}
    }
  
  //cout<<"Six Loops identified"<<endl;
     //Lets do matrix multiplication and summation ( summation U2*U3*U4), first initialize SumProduct3Us to zero
  //cout<<"Starting matrixmultiplication and summation"<<endl;


  SumProduct3Us[0][0] = TComplex(0.0,0.0);
  SumProduct3Us[1][0] = TComplex(0.0,0.0);
  SumProduct3Us[0][1] = TComplex(0.0,0.0);
  SumProduct3Us[1][1] = TComplex(0.0,0.0);
	  //cout<<SumProduct3Us[i][j]<<"\t";

	  //cout<<"\n";
    
  
  TComplex A,B,C; 
  A = TComplex(0.0,0.0);
  B = TComplex(0.0,0.0);
  C = TComplex(0.0,0.0);
  
  for( Int_t p = 0; p<6; p++)
    {
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
	    //  {Z = Zmin + (i*dz);dpdz = TMath::Sqrt(1 - (TMath::Power(beta*K,-2)*TMath::Power(log(Z),2)) );
	    //if( dpdzmax < dpdz ){dpdzmax = dpdz;}}
	    dpdzmax =1.0;
	    
	    while(Success==0)
	      {		    
		double Dummy1 = ObjectMTRand.rand();
		double Dummy2 = ObjectMTRand.rand();
		r1 = Dummy1;r2 = Dummy2;
		Z = Zmin + (Zmax - Zmin)*r1; 
		dpdz = TMath::Sqrt(1.0 - (TMath::Power(beta*K,-2.0)*TMath::Power(log(Z),2.0)) );
		
		//cout<<Dummy1<<"\t"<<Dummy2<<endl;
		//cout<<"Searching"<<endl;
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
	    double Dummy1 = ObjectMTRand.rand();
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
	  Int_t Success =0;
          Double_t lambdaTrial=0.0; 	  
	  while(Success==0)
	    {
	      double Dummy1 = ObjectMTRand.rand();
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
	      r1 = 1.0 - (ObjectMTRand.rand());
	      r2 = 1.0 - (ObjectMTRand.rand());
	      r3 = 1.0 - (ObjectMTRand.rand());
	      
	      lambdaTrialS = - ( ( log(r1) +( TMath::Power(cos(2*PI*r2),2)*log(r3)) )/(2*K*beta) ) ;
	      r = ObjectMTRand.rand();
	      if(r <= TMath::Sqrt(1- lambdaTrialS))
		{
		  Success = 1;   a0= (1.0 - 2*lambdaTrialS); 
		}
	      else{ Success = 0;}
	    }
	  */	  
	  ///////////...End of algorithm given in Book "QCD Lattice".........////////////
	    
	    	    
	    //{cout<<"\n random number a0 = "<<a0<<endl;}
	    
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
          TComplex ConjugateU;   
	  
	  //cout<<"\n\nCreated a0 a1 a2 a3 \t "<<a0<<"\t"<<a1<<"\t"<<a2<<"\t"<<a3<<endl;

	  UCurrentLink[0][0] = TComplex(0.0,0.0);
	  UCurrentLink[1][0] = TComplex(0.0,0.0);
	  UCurrentLink[0][1] = TComplex(0.0,0.0);
	  UCurrentLink[1][1] = TComplex(0.0,0.0);
	  ConjugateU = TComplex(0.0,0.0);

	  
	  for(Int_t i=0; i<2; i++)
            {
	      for(Int_t j =0; j<2; j++)
		{
		  for(Int_t m=0; m<2; m++)
		    {		      		      
		      ConjugateU = TComplex::Conjugate(SumProduct3Us[j][m]);		      
		      //cout<<"Conjugate U "<<ConjugateU<<endl;
		      UCurrentLink[i][j] = UCurrentLink[i][j] +( (a0*ID[i][m] + IYOTA*(a1*Pauli[0][i][m] + a2*Pauli[1][i][m] + a3*Pauli[2][i][m]) )*ConjugateU/K);
		      //OlderVersion
		      //UCurrentLink[i][j] = UCurrentLink[i][j] + ( (a0*ID[i][m] + IYOTA*(a1*Pauli[0][i][m] + a2*Pauli[1][i][m] + a3*Pauli[2][i][m]) )*(TComplex(SumProduct3Us[j][m].Re(),-SumProduct3Us[j][m].Im()))/K);
		    }
		} 
	    }
	  
	  
	  //cout<<(SumProduct3Us[0][0]*SumProduct3Us[1][1] - SumProduct3Us[0][1]*SumProduct3Us[1][0])/(K*K)<<endl;
	      LinksContent[CurrentLink][0] = UCurrentLink[0][0].Re();
	      LinksContent[CurrentLink][3] = UCurrentLink[0][0].Im();
	      LinksContent[CurrentLink][1] = UCurrentLink[0][1].Im();
	      LinksContent[CurrentLink][2] = UCurrentLink[0][1].Re();
	      
		
	      //cout<<"After updating link--Content is a0 a1 a2 a3= "<<LinksContent[CurrentLink][0]<<"\t"<<LinksContent[CurrentLink][1]<<"\t"<<LinksContent[CurrentLink][2]<<"\t"<<LinksContent[CurrentLink][3]<<endl;   
	      //cout<<"Modulus of 4-vector after = "<<TMath::Power(UCurrentLink[0][0].Re(),2.0)+TMath::Power(UCurrentLink[0][0].Im(),2.0) + TMath::Power(UCurrentLink[0][1].Re(),2.0) + TMath::Power(UCurrentLink[0][1].Im(),2.0)<<endl;
	      
	      
 
}
	    



double LatticeCalculation::TracePlaquette(int n1 ,int n2, int n3, int n4)
{
  //cout<<"n1 n2 n3 n4"<<n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<n4<<endl;
  TComplex U1,U2,U3,U4,TraceU;

  
  //cout<<"Coordinate of Link n1 = "<<LinksCoordinate[n1][0]<<LinksCoordinate[n1][1]<<LinksCoordinate[n1][2]<<LinksCoordinate[n1][3]<<LinksDirection[n1]<<"\n";
  //cout<<"LinksContent = "<<LinksContent[n1][0]<<"\t"<<LinksContent[n1][1]<<"\t"<<LinksContent[n1][2]<<"\t"<<LinksContent[n1][3]<<"\t Modulus = "<<pow(LinksContent[n1][0],2) + pow(LinksContent[n1][1],2) + pow(LinksContent[n1][2],2) + pow(LinksContent[n1][3],2)<<endl;

  //cout<<"Coordinate of Link n2 = "<<LinksCoordinate[n2][0]<<LinksCoordinate[n2][1]<<LinksCoordinate[n2][2]<<LinksCoordinate[n2][3]<<LinksDirection[n2]<<"\n";
  //cout<<"LinksContent = "<<LinksContent[n2][0]<<"\t"<<LinksContent[n2][1]<<"\t"<<LinksContent[n2][2]<<"\t"<<LinksContent[n2][3]<<"\t Modulus = "<<pow(LinksContent[n2][0],2) + pow(LinksContent[n2][1],2) + pow(LinksContent[n2][2],2) + pow(LinksContent[n2][3],2)<<endl;

  //cout<<"Coordinate of Link n1 = "<<LinksCoordinate[n3][0]<<LinksCoordinate[n3][1]<<LinksCoordinate[n3][2]<<LinksCoordinate[n3][3]<<LinksDirection[n3]<<"\n";
  //cout<<"LinksContent = "<<LinksContent[n3][0]<<"\t"<<LinksContent[n3][1]<<"\t"<<LinksContent[n3][2]<<"\t"<<LinksContent[n3][3]<<"\t Modulus = "<<pow(LinksContent[n3][0],2) + pow(LinksContent[n3][1],2) + pow(LinksContent[n3][2],2) + pow(LinksContent[n3][3],2)<<endl;

  //cout<<"Coordinate of Link n1 = "<<LinksCoordinate[n4][0]<<LinksCoordinate[n4][1]<<LinksCoordinate[n4][2]<<LinksCoordinate[n4][3]<<LinksDirection[n4]<<"\n";
  //cout<<"LinksContent = "<<LinksContent[n4][0]<<"\t"<<LinksContent[n4][1]<<"\t"<<LinksContent[n4][2]<<"\t"<<LinksContent[n4][3]<<"\t Modulus = "<<pow(LinksContent[n4][0],2) + pow(LinksContent[n4][1],2) + pow(LinksContent[n4][2],2) + pow(LinksContent[n4][3],2)<<endl;
  
  
  
  U1 = TComplex(0.0,0.0); U2 = TComplex(0.0,0.0);
  U3 = TComplex(0.0,0.0); U4 = TComplex(0.0,0.0);TraceU = TComplex(0.0,0.0);

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

double LatticeCalculation::TraceFieldStrengthSquare(int n1,int n2,int n3, int n4)
{
  
  //cout<<"n1 n2 n3 n4 =  "<<n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<n4<<endl; 
  //cout<<"Total links = "<<4*Nt*Nx*Ny*Nz<<endl;

  TComplex PlaquetteMatrix[2][2],FieldStrengthMatrix[2][2],U1,U2,U3,U4,TraceFSquare;

  U1 = TComplex(0.0,0.0); U2 = TComplex(0.0,0.0);
  U3 = TComplex(0.0,0.0); U4 = TComplex(0.0,0.0);TraceFSquare = TComplex(0.0,0.0);

  PlaquetteMatrix[0][0] = TComplex(0.0,0.0);PlaquetteMatrix[0][1] = TComplex(0.0,0.0);
  PlaquetteMatrix[1][0] = TComplex(0.0,0.0);PlaquetteMatrix[1][1] = TComplex(0.0,0.0);

  FieldStrengthMatrix[0][0]= TComplex(0.0,0.0);FieldStrengthMatrix[0][1]= TComplex(0.0,0.0);
  FieldStrengthMatrix[1][0]= TComplex(0.0,0.0);FieldStrengthMatrix[1][1]= TComplex(0.0,0.0);
  

  //Calculating PlaquetteMatrix, Pmunu = U1*U2*U3*U4
  //cout<<"Plaquette Content\n"<<endl;
  //cout<<LinksContent[n1][0]<<"\t"<<LinksContent[n1][1]<<"\t"<<LinksContent[n1][2]<<"\t"<<LinksContent[n1][3]<<endl;
  //cout<<LinksContent[n2][0]<<"\t"<<LinksContent[n2][1]<<"\t"<<LinksContent[n2][2]<<"\t"<<LinksContent[n2][3]<<endl;
  //cout<<LinksContent[n3][0]<<"\t"<<LinksContent[n3][1]<<"\t"<<LinksContent[n3][2]<<"\t"<<LinksContent[n3][3]<<endl;
  //cout<<LinksContent[n4][0]<<"\t"<<LinksContent[n4][1]<<"\t"<<LinksContent[n4][2]<<"\t"<<LinksContent[n4][3]<<endl;
  //cout<<"Absurd LinNumber = "<<n4<<endl;

  
  
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
        {
	  for(int m=0; m<2; m++)
	    {
	      for(int n=0; n<2; n++)
		{
		  for(int p=0; p<2; p++)
		    {
		      U1 = LinksContent[n1][0]*ID[i][m] + IYOTA*(LinksContent[n1][1]*Pauli[0][i][m] + LinksContent[n1][2]*Pauli[1][i][m] + LinksContent[n1][3]*Pauli[2][i][m]);
		      
		      U2 = LinksContent[n2][0]*ID[m][n] + IYOTA*(LinksContent[n2][1]*Pauli[0][m][n] + LinksContent[n2][2]*Pauli[1][m][n] + LinksContent[n2][3]*Pauli[2][m][n]);
		      
		      U3 = LinksContent[n3][0]*ID[p][n] + CIYOTA*(LinksContent[n3][1]*CPauli[0][p][n] + LinksContent[n3][2]*CPauli[1][p][n] + LinksContent[n3][3]*CPauli[2][p][n]);
		      
		      U4 = LinksContent[n4][0]*ID[j][p] + CIYOTA*(LinksContent[n4][1]*CPauli[0][j][p] + LinksContent[n4][2]*CPauli[1][j][p] + LinksContent[n4][3]*CPauli[2][j][p]);
		      
		      PlaquetteMatrix[i][j] = PlaquetteMatrix[i][j] + U1*U2*U3*U4;
		      //cout<<"PlaquettMatrix[i][j]"<<PlaquetteMatrix[i][j]<<endl;
		    }
		}
	    }
	} 
    }   
  
  //Calculating Field-StrengthMatrix
  
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
        {
	  
	  FieldStrengthMatrix[i][j] = (PlaquetteMatrix[i][j] - (TComplex::Conjugate(PlaquetteMatrix[j][i])))/(2.0*IYOTA*BareCouplingConstant*LatticeSpacing*LatticeSpacing);
	  //cout<<"IYOTA,g,a"<<IYOTA<<"\t"<<BareCouplingConstant<<"\t"<<LatticeSpacing<<endl;
	  //cout<<"Fmunu["<<i<<"]["<<j<<"]"<<FieldStrengthMatrix[i][j]<<endl;;
	  //cout<<"PlaquettMatrix["<<i<<"]["<<j<<"]"<<PlaquetteMatrix[i][j]<<endl;
	  //cout<<"Fmunu["<<i<<"]["<<j<<"]"<<FieldStrengthMatrix[i][j]<<endl;
	  
	}
    }
  //cout<<"2*i*g*a^2 = "<<2.0*IYOTA*BareCouplingConstant*LatticeSpacing*LatticeSpacing<<endl;

  //Calculating Trace of SquaredFieldStrengthMatrix
  for(int i=0; i<2; i++)
    {
      for(int m=0; m<2; m++)	
        {
	  TraceFSquare = TraceFSquare + FieldStrengthMatrix[i][m]*FieldStrengthMatrix[m][i];
	}
    }
  
  //cout<<"TraceFSquare is "<<TraceFSquare<<endl;
  //cout<<"TraceFSquare*(-4*g^2*a^4) is "<<TraceFSquare*(-4.0*pow(BareCouplingConstant,2.0)*pow(LatticeSpacing,4.0))<<endl;
  return (TraceFSquare.Re());
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

void AverageCorrelatorF3iF3iMinusF4iF4iSU2(double Beta)
{
  int StartTime = time(NULL);

  //const double PI = 3.1415926; //Value of pi
  //const int NIteration = 10;
  //srand(time(NULL));
  LatticeCalculation LC;
 
  
  //..............Input Parameter.....Lattice.....
  LC.PI = 3.1415926;
  LC.Nt = 3; //Discretizing each dimension in N-points 
  LC.Nx = 12;
  LC.Ny = 12;
  LC.Nz = 12;
  LC.LatticePoints = LC.Nt*LC.Nx*LC.Ny*LC.Nz; //                            
  LC.TotalLinks = 4*LC.Nt*LC.Nx*LC.Ny*LC.Nz; // 4 space-time dimension      
  LC.NLinksX = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in X-direction        
  LC.NLinksY = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Y-direction        
  LC.NLinksZ = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Z-direction        
  LC.NLinksT = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in T-direction        
  LC.Run = 40;                //Number of iteration to get stable gauge field configuration                                     
  LC.TotalConfiguration = 40; //Number of gauge field configuartion over which average of (F3iF3i-F4iF4i) is computed         
  LC.TotalTemperatureRun = 1; //Number of different values of Temperature to see the dependence of correlation function         
  LC.lambdaLattice = 5.3; //(in MeV) A parameter that goes in LatticeSpacing formular using renormalization group equation 
  
  LC.DefineDynamicArrays();
  LC.DefineComplexMatrix();
  //..............End of Input Parameter........Lattice
  cout<<"Program started"<<endl;
  cout<<"Dimension of Lattice used is  \t" <<LC.Nt<<"*"<<LC.Nx<<"*"<<LC.Ny<<"*"<<LC.Nz<<", So Total Links are "<<LC.TotalLinks<<endl;
  cout<<"\t"<< LC.LatticeSize[0]<<","<<LC.LatticeSize[1]<<","<<LC.LatticeSize[2]<<","<<LC.LatticeSize[3]<<endl;
  
  
  //ofstream f1;
  //f1.open("Coordinate.txt",ios::out);
  //f1.open("Average.txt",ios::out);
  //f1.open("Plaquette_wrt_beta_N6_I20.txt",ios::app);
  //Char_t Var[70]; int INTBeta = Beta*1000;
  //sprintf(Var,"10Run/Correlator_F3iF3i_F4iF4i__wrt_beta_%d_Nt12_Ns12_I5000.txt",INTBeta);
  //f1.open(Var,ios::app);    

//Labeling each link with its 4 space-time coordinate and a direction
  LC.LabelLinksCoordinate();
  
  //cout<<"Line 807 This is working"<<endl;
  
     
  int TemperatureLoop =0;
  int CountPlane41,CountPlane42, CountPlane31,CountPlane32;
  double SumFieldStrengthCorrelation41,SumFieldStrengthCorrelation42,SumFieldStrengthCorrelation31,SumFieldStrengthCorrelation32,SumFieldStrengthCorrelationPerPlaquette;
  double *AverageFieldStrengthCorrelation = new double[LC.TotalTemperatureRun];
  double *ScaledAverageFieldStrengthCorrelation = new double[LC.TotalTemperatureRun];
  double *TEMPERATURE = new double[LC.TotalTemperatureRun];
  double *BETA = new double[LC.TotalTemperatureRun];
  double *BetaAverageCorrelator = new double[LC.TotalTemperatureRun];
  double *BetaAvgCorrelatorError = new double[LC.TotalTemperatureRun];

  double *TraceF31F31PerPlaquette = new double[LC.TotalConfiguration];
  double *TraceF32F32PerPlaquette = new double[LC.TotalConfiguration];
  double *TraceF41F41PerPlaquette = new double[LC.TotalConfiguration];
  double *TraceF42F42PerPlaquette = new double[LC.TotalConfiguration];
  
  double AverageTraceF31F31PerPlaquette;
  double AverageTraceF32F32PerPlaquette;
  double AverageTraceF41F41PerPlaquette;
  double AverageTraceF42F42PerPlaquette;

  double SumDeviationSquareF31F31;
  double SumDeviationSquareF32F32;
  double SumDeviationSquareF41F41;
  double SumDeviationSquareF42F42;
  
  double ErrorF31F31, ErrorF32F32, ErrorF41F41, ErrorF42F42;

  double *ScaledAvgFSCorrelationError = new double[LC.TotalTemperatureRun];
  
  
  for(TemperatureLoop = 0; TemperatureLoop<LC.TotalTemperatureRun; TemperatureLoop++)  
    {      
      LC.beta = Beta;//1.80 + 0.1*TemperatureLoop;
      LC.BareCouplingConstant = TMath::Sqrt(4.0/LC.beta);
      LC.LatticeSpacing = TMath::Power(11.0*TMath::Power(LC.BareCouplingConstant/LC.PI,2.0)/24.0,-51.0/121.0)*exp(-12.0*TMath::Power(LC.PI/LC.BareCouplingConstant,2.0)/11.0)/(LC.lambdaLattice);
      LC.temperature = 1/(LC.Nt*LC.LatticeSpacing); 
      TEMPERATURE[TemperatureLoop] = LC.temperature;
      BETA[TemperatureLoop] = LC.beta;

      cout<<"Coupling is "<<LC.BareCouplingConstant<<endl;
      cout<<"Beta (4/g^2)= "<< LC.beta<<endl; 
      cout<<"LatticeSpacing "<<LC.LatticeSpacing<<endl;
      cout<<"Temperature is "<<LC.temperature<<endl;
      cout<<"(nt)^4 "<<pow(LC.Nt,4.0)<<endl;
      cout<<"1/T^4 "<<pow(1/LC.temperature,4.0)<<endl;;
      cout<<"1/(4*g^2*a^4) "<<1/pow(2*LC.BareCouplingConstant*pow(LC.LatticeSpacing,2),2.0)<<endl; 
      cout<<"1/(4*g^2*a^4*T^4) "<<1/pow(2*LC.BareCouplingConstant*pow(LC.LatticeSpacing*LC.temperature,2.0),2.0)<<endl;

      SumFieldStrengthCorrelationPerPlaquette = 0.0;
      SumDeviationSquareF31F31 = 0.0;
      SumDeviationSquareF32F32 = 0.0;
      SumDeviationSquareF41F41 = 0.0;
      SumDeviationSquareF42F42 = 0.0;
      
      
      for(int i=0; i<LC.TotalConfiguration; i++)
	{ TraceF31F31PerPlaquette[i] = 0.0;
	  TraceF32F32PerPlaquette[i] = 0.0; 
	  TraceF41F41PerPlaquette[i] = 0.0;
	  TraceF42F42PerPlaquette[i] = 0.0;
	}

      AverageTraceF31F31PerPlaquette = 0.0;
      AverageTraceF32F32PerPlaquette = 0.0;
      AverageTraceF41F41PerPlaquette = 0.0;
      AverageTraceF42F42PerPlaquette = 0.0;

 


      //Setting all link variables to Unity
      for(Int_t i = 0; i<LC.TotalLinks; i++)
	{ 
	  LC.LinksContent[i][0]=1.0;      //a0
	  LC.LinksContent[i][1]=0.0;      //a1
	  LC.LinksContent[i][2]=0.0;      //a2
	  LC.LinksContent[i][3]=0.0;      //a3
	}

      
      //Constructing first stable gauge field configuration
      
      for(Int_t Iter1=0; Iter1<LC.Run; Iter1++)
	{
	  LC.Iteration = Iter1;    
	  cout<<"Performing Iteration "<<Iter1<<endl;  
	  //Touching each internal link with heat bath 
	  for(Int_t ic = 0 ; ic< LC.TotalLinks ; ic++)
	    {
	      //cout<<"Processing Link"<<ic<<endl;
	      LC.UpdateLink(ic); //Apply heat bath algorithm to each internal link and apply periodic boundary condition
	      //cout<<"Content of Link "<<LC.LinksContent[ic][0]<<"\t"<<LC.LinksContent[ic][1]<<"\t"<<LC.LinksContent[ic][2]<<"\t"<<LC.LinksContent[ic][3]<<endl;   
	    }
	}
      

      //Now We have stable gauge field configuration
      for(Int_t Iter=0; Iter<LC.TotalConfiguration; Iter++)
        {
          LC.Iteration = LC.Iteration + Iter;
          cout<<"\nWorking on Configuration "<<Iter<<endl;

          SumFieldStrengthCorrelation41 = 0.0;SumFieldStrengthCorrelation42 =0.0;
	  SumFieldStrengthCorrelation31 = 0.0;SumFieldStrengthCorrelation32 =0.0;
	  CountPlane41 = 0;
	  CountPlane42 = 0;
	  CountPlane31 = 0;
	  CountPlane32 = 0;

	  
	  //Touching each internal link with heat bath
          for(Int_t ic = 0 ; ic< LC.TotalLinks ; ic++)
            {
	      //cout<<"Processing Link"<<ic<<endl;                                                                                                            
              LC.UpdateLink(ic); //Apply heat bath algorithm to each internal link and apply periodic boundary condition
	    }
	  
	  //Evaluate Trace of Field-Strength-Field-Strength Matrix
	  
	  Int_t xt[4] ={0,0,0,0};
	  Int_t d1=0,d2=0,p1=0,p2=0;
	  for(Int_t id = 0; id <4; id++)
	    {

	      if(id==0){d1=0; d2=1; p1=2; p2=3;}
	      if(id==1){d1=0; d2=2; p1=1; p2=3;}
	      if(id==2){d1=1; d2=3; p1=0; p2=2;}
	      if(id==3){d1=2; d2=3; p1=0; p2=1;}

	      for(Int_t i=0; i<LC.LatticeSize[p1]; i++)  //axis perpendicular to V-axis and H-axis 
		{
		  for(Int_t j=0; j<LC.LatticeSize[p2]; j++)//axis perpendicular to V-axis and H-axis
		    {
		      for(Int_t k=0; k<LC.LatticeSize[d1]; k++)  //Vertical-axis                                             
			{
			  for(Int_t m=0; m<LC.LatticeSize[d2]; m++) //Horizontal axis                                             
			    {
			      //if(i==1 && j==1 && 1 && m ==1 && k==1)
				{
			      xt[d1]=k;
			      xt[d2]=m;
			      xt[p1]=i;
			      xt[p2]=j;
			    
			      //cout<<"Going to calculate Fmunu using Plaquette \n\n";
			      
			      //if(xt[p1]==LC.LatticeSize[p1]-1 && xt[p2]==LC.LatticeSize[p2]-1)
			      //{xt[p1]=0;xt[p2]=0;}
			      //if(xt[p1]==LC.LatticeSize[p1]-1){xt[p1]=0;}
			      //if(xt[p2]==LC.LatticeSize[p2]-1){xt[p2]=0;}
			      
			      
			      Int_t Link4 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d1];
			      Int_t Link1 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d2];

			      xt[d2]=xt[d2]+1;
			      
			      if(xt[d2] == LC.LatticeSize[d2])
				{ xt[d2]=0;
				}			      				  
			      Int_t Link2 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d1];
			      
			      xt[d2] = m;
			      
			      xt[d1] = xt[d1]+1;
			      if( xt[d1] == LC.LatticeSize[d1])
                                {xt[d1]=0;
                                }
			      Int_t Link3 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d2]; 
			      //cout<<"We have formed the plaquette square \n\n";
			      //cout<<LC.LinksContent[Link1][0]<<"\t"<<LC.LinksContent[Link1][1]<<"\t"<<LC.LinksContent[Link1][2]<<"\t"<<LC.LinksContent[Link1][3]<<endl;			     
			      
			      if(id==0 )
				{				  
				  CountPlane41 = CountPlane41 + 1;
				  SumFieldStrengthCorrelation41 = SumFieldStrengthCorrelation41 + LC.TraceFieldStrengthSquare(Link1,Link2,Link3,Link4);
				  //cout<<"In Plane (4,1):  summation TraceUmunu^2 = "<<SumFieldStrengthCorrelation41*(-4.0)*pow(LC.BareCouplingConstant*LC.LatticeSpacing*LC.LatticeSpacing,2.0)/(CountPlane41)<<endl;
				}

			      if(id==1 )
                                {
                                  CountPlane42 = CountPlane42 + 1;
                                  SumFieldStrengthCorrelation42 = SumFieldStrengthCorrelation42 + LC.TraceFieldStrengthSquare(Link1,Link2,Link3,Link4);
                                  //cout<<"In Plane (4,2):  summation TraceUmunu^2 = "<<SumFieldStrengthCorrelation42*(-4.0)*pow(LC.BareCouplingConstant*LC.LatticeSpacing*LC.LatticeSpacing,2.0)/(CountPlane42)<<endl;
                                }

			      if(id==2 )
                                {
				  CountPlane31 = CountPlane31 + 1;
				  SumFieldStrengthCorrelation31 = SumFieldStrengthCorrelation31 + LC.TraceFieldStrengthSquare(Link1,Link2,Link3,Link4);
				  //cout<<"In Plane (3,1):  summation TraceUmunu^2 = "<<SumFieldStrengthCorrelation31*(4.0)*pow(LC.BareCouplingConstant*LC.LatticeSpacing*LC.LatticeSpacing,2.0)/(CountPlane31)<<endl;  
				}
			      
			      if(id==3 )
                                {
                                  CountPlane32 = CountPlane32 + 1;
                                  SumFieldStrengthCorrelation32 = SumFieldStrengthCorrelation32 + LC.TraceFieldStrengthSquare(Link1,Link2,Link3,Link4);
                                  //cout<<"In Plane (3,2):  summation TraceUmunu^2 = "<<SumFieldStrengthCorrelation32*(4.0)*pow(LC.BareCouplingConstant*LC.LatticeSpacing*LC.LatticeSpacing,2.0)/(CountPlane32)<<endl;
                                }
			      //cout<<"Summation (FmunuFmnu)"<<SumFieldStrengthCorrelation<<endl;
			      
				}
			    }
			}
		    }
		}
	    }

	  
	  TraceF31F31PerPlaquette[Iter] = SumFieldStrengthCorrelation31/CountPlane31;
	  TraceF32F32PerPlaquette[Iter] = SumFieldStrengthCorrelation32/CountPlane32;
	  TraceF41F41PerPlaquette[Iter] = SumFieldStrengthCorrelation41/CountPlane41;
	  TraceF42F42PerPlaquette[Iter] = SumFieldStrengthCorrelation42/CountPlane42;

	  SumFieldStrengthCorrelationPerPlaquette = SumFieldStrengthCorrelationPerPlaquette + (TraceF31F31PerPlaquette[Iter] + TraceF32F32PerPlaquette[Iter] - TraceF41F41PerPlaquette[Iter] - TraceF42F42PerPlaquette[Iter]);
          


	  cout<<"Trace(F(4,1)^2)/(N41*T^4)  && (4,2) "<<TraceF41F41PerPlaquette[Iter]/pow(LC.temperature,4.0)<<"\t"<<TraceF42F42PerPlaquette[Iter]/pow(LC.temperature,4.0)<<endl;
	  cout<<"Trace(F(3,1)^2)/(N31*T^4)  && (3,2) "<<TraceF31F31PerPlaquette[Iter]/pow(LC.temperature,4.0)<<"\t"<<TraceF32F32PerPlaquette[Iter]/pow(LC.temperature,4.0)<<endl;
	  cout<<"N41,N42,N31,N32 = "<<CountPlane41<<"\t"<<CountPlane42<<"\t"<<CountPlane31<<"\t"<<CountPlane32<<endl;
	  cout<<"Tr(F3i^2)/(2*N31*T^4) - Tr(F4i^2)/(2*N41*T^4)=  "<<(TraceF31F31PerPlaquette[Iter]+TraceF32F32PerPlaquette[Iter]-TraceF41F41PerPlaquette[Iter]-TraceF42F42PerPlaquette[Iter])/(2.0*pow(LC.temperature,4.0))<<endl;
         cout<<"ScaledAverageFmunu^2 So far = "<<SumFieldStrengthCorrelationPerPlaquette/(2*pow(LC.temperature,4.0)*(Iter+1))<<endl;
	 cout<<"So far a^4*Sigma(F3iF3i-F4iF4i)/2 = "<<pow(LC.LatticeSpacing,4.0)*SumFieldStrengthCorrelationPerPlaquette/(2.0*(Iter+1))<<endl;
	  
	}
      
      for(int i=0; i<LC.TotalConfiguration; i++)
	{
	  AverageTraceF31F31PerPlaquette = AverageTraceF31F31PerPlaquette + (TraceF31F31PerPlaquette[i]/LC.TotalConfiguration);
	  AverageTraceF32F32PerPlaquette = AverageTraceF32F32PerPlaquette + (TraceF32F32PerPlaquette[i]/LC.TotalConfiguration);
	  AverageTraceF41F41PerPlaquette = AverageTraceF41F41PerPlaquette + (TraceF41F41PerPlaquette[i]/LC.TotalConfiguration);
	  AverageTraceF42F42PerPlaquette = AverageTraceF42F42PerPlaquette + (TraceF42F42PerPlaquette[i]/LC.TotalConfiguration);
	}
      
      for(int i=0; i<LC.TotalConfiguration; i++)
        {
          SumDeviationSquareF31F31 =   SumDeviationSquareF31F31 + pow(AverageTraceF31F31PerPlaquette - TraceF31F31PerPlaquette[i],2.0); 
	  SumDeviationSquareF32F32 =   SumDeviationSquareF32F32 + pow(AverageTraceF32F32PerPlaquette - TraceF32F32PerPlaquette[i],2.0);
	  SumDeviationSquareF41F41 =   SumDeviationSquareF41F41 + pow(AverageTraceF41F41PerPlaquette - TraceF41F41PerPlaquette[i],2.0);
          SumDeviationSquareF42F42 =   SumDeviationSquareF42F42 + pow(AverageTraceF42F42PerPlaquette - TraceF42F42PerPlaquette[i],2.0);
	}
      
      ErrorF31F31 = TMath::Sqrt(SumDeviationSquareF31F31/(LC.TotalConfiguration-1.0))/TMath::Sqrt(LC.TotalConfiguration);
      ErrorF32F32 = TMath::Sqrt(SumDeviationSquareF32F32/(LC.TotalConfiguration-1.0))/TMath::Sqrt(LC.TotalConfiguration);
      ErrorF41F41 = TMath::Sqrt(SumDeviationSquareF41F41/(LC.TotalConfiguration-1.0))/TMath::Sqrt(LC.TotalConfiguration);
      ErrorF42F42 = TMath::Sqrt(SumDeviationSquareF42F42/(LC.TotalConfiguration-1.0))/TMath::Sqrt(LC.TotalConfiguration);
            
      ScaledAverageFieldStrengthCorrelation[TemperatureLoop] = (AverageTraceF31F31PerPlaquette  + AverageTraceF32F32PerPlaquette - AverageTraceF41F41PerPlaquette - AverageTraceF42F42PerPlaquette  )/(2*pow(LC.temperature,4.0));

      BetaAverageCorrelator[TemperatureLoop] = pow(LC.LatticeSpacing,4.0)*(AverageTraceF31F31PerPlaquette  + AverageTraceF32F32PerPlaquette - AverageTraceF41F41PerPlaquette - AverageTraceF42F42PerPlaquette)/2.0;
      
      ScaledAvgFSCorrelationError[TemperatureLoop] =  TMath::Sqrt( pow(ErrorF31F31,2.0) + pow(ErrorF32F32,2.0) + pow(ErrorF41F41,2.0) + pow(ErrorF42F42,2.0) )/(2*pow(LC.temperature,4.0));
      
      BetaAvgCorrelatorError[TemperatureLoop] = pow(LC.LatticeSpacing,4.0)*TMath::Sqrt( pow(ErrorF31F31,2.0) + pow(ErrorF32F32,2.0) + pow(ErrorF41F41,2.0) + pow(ErrorF42F42,2.0) )/2.0;
      //cout<<"Number of Plaquette in (3,1) plane = "<<CountPlane31<<endl;
      //cout<<"Number of Plaquette in (3,2) plane = "<<CountPlane32<<endl;
      //cout<<"Number of Plaquette in (4,1) plane = "<<CountPlane41<<endl;
      //cout<<"Number of Plaquette in (4,2) plane = "<<CountPlane42<<endl;
      //cout<<LC.temperature<<"\t"<<ScaledAverageFieldStrengthCorrelation[TemperatureLoop]<<endl;
      //f1<<TEMPERATURE[TemperatureLoop]<<"\t"<<ScaledAverageFieldStrengthCorrelation[TemperatureLoop]<<"\t"<<ScaledAvgFSCorrelationError[TemperatureLoop]<<"\t"<<ScaledAvgFSCorrelationError[TemperatureLoop]*100.0/ScaledAverageFieldStrengthCorrelation[TemperatureLoop];
      //f1<<"\t"<<LC.beta<<"\t"<<BetaAverageCorrelator[TemperatureLoop]<<"\t"<<BetaAvgCorrelatorError[TemperatureLoop]<<"\t"<<BetaAvgCorrelatorError[TemperatureLoop]*100.0/BetaAverageCorrelator[TemperatureLoop]<<endl;
      
    }

  /*        
 for(Int_t id= 0; id<LC.TotalTemperatureRun; id++)
   {
     cout<<"Scaled Correlation function at temperature "<<TEMPERATURE[id]<<"\t is \t"<<ScaledAverageFieldStrengthCorrelation[id]<<"\t+-"<<ScaledAvgFSCorrelationError[id]<<"\t"<<ScaledAvgFSCorrelationError[id]*100.0/ScaledAverageFieldStrengthCorrelation[id]<<endl;
     cout<<"Beta Plot value "<<LC.beta<<"\t"<<BetaAverageCorrelator[id]<<"\t"<<BetaAvgCorrelatorError[id]<<"\t"<<BetaAvgCorrelatorError[id]*100.0/BetaAverageCorrelator[id]<<endl;
     //f1<<id<<"\t"<<AveragePlaquette[id]<<endl;
     }
  */
   cout<<"Ending program"<<endl;
  
   /*    
//Plotting the ScaledAverageFieldStrengthCorrelation vs Temperature
    
   TCanvas *c1 = new TCanvas("C1","c1");    
   TGraph *GAvg = new TGraph(LC.TotalTemperatureRun,TEMPERATURE,ScaledAverageFieldStrengthCorrelation);

 GAvg->Draw("A*");
 GAvg->GetXaxis()->SetRangeUser(-1.,LC.Run);
 GAvg->GetYaxis()->SetRangeUser(0.0,3.0);
 GAvg->SetTitle("FieldStrength Correlation vsTemperature");
 GAvg->GetYaxis()->SetTitle("Scaled Average FieldStrength Correlation");
 GAvg->GetXaxis()->SetTitle("Temperature (MeV)");
 GAvg->GetYaxis()->CenterTitle();
 GAvg->GetXaxis()->CenterTitle();
 c1->SaveAs("F3iF3iT4_3_12_12_12_500.root");

 TCanvas *c2 = new TCanvas("c2","c2");
 TGraph *G2 = new TGraph(LC.TotalTemperatureRun,BETA,BetaAverageCorrelator);
 G2->Draw("A*");
 G2->GetYaxis()->SetRangeUser(0.,0.1);
 G2->GetXaxis()->SetRangeUser(1.0,3.0);
 G2->SetTitle("A^4*FieldStrength Correlation vs Beta");
 G2->GetYaxis()->SetTitle("Average a^4*FieldStrength Correlation");
 G2->GetXaxis()->SetTitle("Beta");
 G2->GetYaxis()->CenterTitle();
 G2->GetXaxis()->CenterTitle();
 c2->SaveAs("a4F3i4i_3_12_12_12_500.root");
    //   f1<<LC.beta<<"\t"<<AveragePlaquette[LC.Run-1]<<endl;
   //f1.close();  
   */  int EndTime = time(NULL);
  int Hour = (EndTime-StartTime)/3600;
  int Minute = ((EndTime-StartTime)/60)-Hour*60;
  int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
  cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
   
}
