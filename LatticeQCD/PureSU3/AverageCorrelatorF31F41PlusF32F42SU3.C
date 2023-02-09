//This program calculates the average plaquette for SU(3) gauge field for N*N*N*N Lattice with Periodic BC's
//  Following is lattice structure for 3*4 (Notice the missing link at boundary;they are to incorporate the periodic BC )
//
//     |_|_|_|_
//     |_|_|_|_
//     |_|_|_|_
//
//We use cpp function "rand()" to generate random number
//We use 5 variables to store the coordinate of the link (first four for t,x,y,z and last-one for its direction; We use 0 to point in t-dir, 1 (x-dir).. ) 
//We use Heat-bath algorithm to update each link
//To compile this code (for beta =2.3) in terminal type:  root -l "AveragePlaquetteSU3.C(2.3)" 
//Contact Amit Kumar for questions (amitkr2410@gmail.com or kumar.amit@wayne.edu)

#include"iostream"
#include"fstream"
#include"TMath.h"
//#include"stdlib.h" //for srand and rand
#include"time.h"  //for time() function
#include"TCanvas.h"
#include"TGraphErrors.h"
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
  int TotalConfiguration; //Number of gauge field configuartion over which average of (F3iF3i-F4iF4i) is computed 
  int TotalTemperatureRun; //Number of different values of Temperature to see the dependence of correlation function  
  int* LatticeSize;//LatticeSize[4] = {Nt,Nx,Ny,Nz}
  double lambdaLattice; //(in MeV) A parameter that goes in LatticeSpacing formula using renormalization group equation    

  double beta, temperature, LatticeSpacing, BareCouplingConstant;

  MTRand ObjectMTRand;//MTRand Class is defined in "MTwister.h" which has rand() function defined
  
  int CurrentLink, Iteration;

  TComplex*** LinksContent;//To define LinksContent[TotalLinks][3][3]; 
  int** LinksCoordinate; //To define LinksCoordinate[TotalLinks][4];
  int* LinksDirection;  //To define LinksDirection[TotalLinks];

  int***** LinksNumber; //LinkNumber[t,x,y,z,d]  

  TComplex ID[2][2],Pauli[3][2][2],CPauli[3][2][2],IYOTA,CIYOTA; //Complex Matrix
  TComplex SumProduct3Us[3][3];  //Complex matrix Staple U2*U3*U4 (It is staple)

  void DefineDynamicArrays();
  void LabelLinksCoordinate();
  Int_t OutSideBoundary(Int_t,Int_t ,Int_t, Int_t,Int_t);
  void UpdateLink(Int_t);
  //Int_t SearchLinkNumber(Int_t, Int_t, Int_t, Int_t, Int_t);
  void DefineComplexMatrix();
  double TracePlaquette(int,int,int,int);
  double TraceFieldStrengthSquare(int,int,int,int);
  double TraceFieldStrengthCrossedTerm(int* ,int* ,int* ,int*);
};

void LatticeCalculation::DefineDynamicArrays()
{

  LinksDirection = new int [TotalLinks];
  LinksCoordinate = new int* [TotalLinks];
  LinksContent = new TComplex** [TotalLinks];

  LatticeSize = new int [4];

  for(int i=0; i<TotalLinks; i++)
    {
      LinksCoordinate[i] = new int [4];
      LinksContent[i] = new TComplex* [3];
    }
  
  for(int i=0; i<TotalLinks; i++)
    {
      for(int j=0; j<3; j++)
	{
	  LinksContent[i][j] = new TComplex [3];
	}
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
		    
		    LinksNumber[a][b][c][d][3] = CurrentLink;

		    //cout<<"\n txyzd \t"<<a<<","<<b<<","<<c<<","<<d<<","<<3<<endl;
		    //cout<<"Links Number = "<<  LinksNumber[a][b][c][d][3]<<"\n";
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

		    //cout<<" txyzd \t"<<a<<","<<b<<","<<c<<","<<d<<","<<2<<endl;
                    //cout<<"Links Number = "<<  LinksNumber[a][b][c][d][2]<<"\n";
		    // std::cout<<LinksCoordinate[CurrentLink+NLinks][0]<<LinksCoordinate[CurrentLink][1]<<LinksCoordinate[CurrentLink][2]<<LinksCoordinate[CurrentLink][3]<<endl;
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
	      LinkNumberStaple[count][i2] = LinksNumber[StapleCoordinate[p][ra][0]][ StapleCoordinate[p][ra][1]][ StapleCoordinate[p][ra][2]][ StapleCoordinate[p][ra][3]][ StapleCoordinate[p][ra][4] ];
	      //cout<<"LinkNumberFound "<<LinkNumberStaple[count][i2]<<endl;
	      // cout<<"Content is a0 a1 a2 a3 "<<LinksContent[LinkNumberStaple[count][i2]][0]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][1]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][2]<<"\t"<<LinksContent[LinkNumberStaple[count][i2]][3]<<endl;
	}
    }
  
  //cout<<"Six Loops identified"<<endl;
  //Lets do matrix multiplication and summation ( summation U2*U3*U4), first initialize SumProduct3Us to zero
  //cout<<"Starting matrixmultiplication and summation"<<endl;

  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  SumProduct3Us[i][j] = TComplex(0.0,0.0);
	}
    }
  
  
  

  TComplex A,B,C; 
  A=TComplex(0.0,0.0);
  B=TComplex(0.0,0.0);
  C=TComplex(0.0,0.0);

  for( Int_t p = 0; p<6; p++)
    {
      //cout<<"LinksContent ="<<LinksContent[LinkNumberStaple[p][2]]<<endl;
      for(Int_t i = 0; i<3; i++)
	{
	  for(Int_t j = 0; j<3; j++)
	    {
	      for(Int_t m = 0; m<3; m++)
		{
		  for(Int_t n = 0; n<3; n++)
		    {
		      if(p==0 || p==2 || p==4)
			{
			  
			  A = LinksContent[LinkNumberStaple[p][2]][i][m];

			  B = TComplex::Conjugate(LinksContent[LinkNumberStaple[p][1]][n][m]); //(Dagger)
												 
			  C = TComplex::Conjugate(LinksContent[LinkNumberStaple[p][0]][j][n]); //(Dagger)
			}
		      
		      if(p==1 || p==3 || p==5)
			{
			  A = TComplex::Conjugate(LinksContent[LinkNumberStaple[p][2]][m][i]);//(Dagger)
 
			  B = TComplex::Conjugate(LinksContent[LinkNumberStaple[p][1]][n][m]); // (Dagger)

			  C = LinksContent[LinkNumberStaple[p][0]][n][j];
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
  
  //Finding Block_matrix R
  // R = U(0)*SumProduct3Us;  or U(0)*Staple[3][3]; where U(0)= U matrix at Current Link; SumProduct= Staple
  // R is 3*3 Complex matrix
  TComplex UNew[3][3]={0.0}, UOld[3][3]={0.0};
  TComplex ROld[3][3]={0.0}, RNew[3][3]={0.0};
  TComplex AA[3][3]={0.0};

  //cout<<"U(0)*Staple Matrix at _v=0_ = \n";

  for(int v = 0; v<3; v++)
    {

      if(v==0)
	{
	  for(int i = 0; i<3; i++)
	    {
	      for(int j = 0; j<3; j++)
		{
		  for(int m = 0; m<3; m++)
		    {	     
		      RNew[i][j] = RNew[i][j] + LinksContent[CurrentLink][i][m]*SumProduct3Us[m][j];

		    }
		  ROld[i][j] = RNew[i][j];
		  UNew[i][j] = LinksContent[CurrentLink][i][j];
		  UOld[i][j] = UNew[i][j];

		  //cout<<ROld[i][j]<<"\t";
		}
	      //cout<<"\n"<<endl;
	    }
	}


      
  //Generating A(1) using R matrix
  //Generate Block matric corresponding to R-matrix and call it r(2,2)

  Double_t r[4]={0.0}, rVec=0.0;

  if(v==0)
    {
      r[0] = (RNew[0][0].Re() + RNew[1][1].Re())/2.0;
      r[3] = (RNew[0][0].Im() - RNew[1][1].Im())/2.0 ;
      r[1] = (RNew[0][1].Im() + RNew[1][0].Im())/2.0;
      r[2] = (RNew[0][1].Re() - RNew[1][0].Re())/2.0 ;
    }

  if(v==1)
    {
      r[0] = (RNew[1][1].Re()  +  RNew[2][2].Re())/2.0;
      r[3] = (RNew[1][1].Im() - RNew[2][2].Im())/2.0 ;
      r[1] = (RNew[1][2].Im() + RNew[2][1].Im())/2.0;
      r[2] = (RNew[1][2].Re() - RNew[2][1].Re())/2.0 ;
    }

  if(v==2)
    {
      r[0] = (RNew[0][0].Re() + RNew[2][2].Re())/2.0;
      r[3] = (RNew[0][0].Im() - RNew[2][2].Im())/2.0;
      r[1] = (RNew[0][2].Im() + RNew[2][0].Im())/2.0;
      r[2] = (RNew[0][2].Re() - RNew[2][0].Re())/2.0 ;
    }


  rVec = TMath::Sqrt( TMath::Power(r[1],2.0) + TMath::Power(r[2],2.0) + TMath::Power(r[3],2.0) );
  Double_t chi = TMath::Sqrt( TMath::Power(rVec,2.0) + TMath::Power(r[0],2.0) );

	  //Create pdf and generate x0 x1 x2 x3 with apporpiate weight 
	  
          Double_t x0, x1, x2, x3,Theta, Phi;
	  
	  /////////////////////............Differential Sampling for a0.......////////////////////////
	  
	  Double_t x = 1.0;
	  Double_t Zmin, Zmax, Z, dpdz, dpdzmax=0.0, dz;
	  Double_t r1 , r2;
	  Int_t Success =0;
	  Zmin = TMath::Exp(-x*beta*chi*2.0/3.0);
	  Zmax = TMath::Exp(x*beta*chi*2.0/3.0);

	  //dz = (Zmax-Zmin)/(100000.0-1.0);                                                                                                                
	  //for(Int_t i=0; i< 100000; i++)                                                                                                                  
	  //  {Z = Zmin + (i*dz);dpdz = sqrt(1 - (pow(beta*K,-2)*pow(log(Z),2)) );                                                                          
	  //if( dpdzmax < dpdz ){dpdzmax = dpdz;}}                                                                                                          
	  dpdzmax =1.0;
	  Int_t SearchCount = 0;//cout<<"Search for x0 started "<<endl;                                                                                     
	  while(Success==0)
	    {
	      Double_t Dummy1 = ObjectMTRand.rand();
	      Double_t Dummy2 = ObjectMTRand.rand();
	      r1 = Dummy1;r2 = Dummy2;
	      Z = Zmin + (Zmax - Zmin)*r1;
	      dpdz = TMath::Sqrt(1 - (TMath::Power(3.0/(2.0*beta*chi),2.0)*TMath::Power(log(Z),2)) );

	      SearchCount = SearchCount + 1;
	      //cout<<"Search count ="<<SearchCount<<endl;                                                                                                  
	      if( dpdz > (dpdzmax*r2) )
		{
		  Success =1;
		}
	      else
		{Success =0;}
	    }
	  
	  x0 = log(Z)*3.0/(2.0*beta*chi);

	    /////////////.......End of differential sampling algorithm for a0...///////////
	    
              
	  ///////////.........Algorithm as defined in QCD Lattice Book............//////////////
	  /*
	  Double_t r1,r2,r3,r4; Double_t lambdaTrialS=0.0;
	  
          Int_t Success=0;
          while(Success==0)                                                                                                                                  
            {       
	      r1 = 1.0 - ObjectMTRand.rand();               
	      r2 = 1.0 - ObjectMTRand.rand();
	      r3 = 1.0 - ObjectMTRand.rand(); 
	      lambdaTrialS = - ( ( log(r1) +( TMath::Power(cos(2.0*PI*r2),2.0)*log(r3)) )/(4.0*beta*r[0]/3.0 ));
	      r4 = ObjectMTRand.rand();				
	      
	      if(r4 <= (TMath::Sqrt(1- lambdaTrialS)*TMath::Exp(-4.0*rVec*beta*TMath::Sqrt(lambdaTrialS)*(TMath::Sqrt(1- lambdaTrialS))/3.0 ) ) )
	      {                                                                                                                     
		Success = 1;   a0= (1.0 - 2.0*lambdaTrialS);
	      }	      
              else{ Success = 0;}
	      
            }                      
	 */
	  ///////////.........End of the gaussian algorithm.......///////
	  //{cout<<"\n random number x0 = "<<x0<<endl;}
	     
	    Double_t Dummy1 = (ObjectMTRand.rand())*2.0;
	    Theta =  TMath::ACos( Dummy1 - 1.0);
	    Double_t Dummy2 = ObjectMTRand.rand();                                                                                                               
	               
	    Phi = Dummy2*2.0*PI;
	    Double_t xVec = TMath::Sqrt(1.0- x0*x0 );
	    
	    x1 = xVec*cos(Theta);
	    x2 = xVec*sin(Theta)*cos(Phi);
	    x3 = xVec*sin(Theta)*sin(Phi);
	    //cout<<"Theta =\t "<<Theta*180.0/PI;
	    //cout<<"\t Modulus of 4-vector x = "<<(x0*x0 + x1*x1 + x2*x2 + x3*x3)<<endl;
	    
	    Double_t a0, a1, a2, a3;
	    for(int i = 0; i<3; i++)
	      {for(int j = 0; j<3; j++)
		  { 
		    ROld[i][j] = RNew[i][j];
		    RNew[i][j] = 0.0;

		    UOld[i][j] = UNew[i][j];
		    UNew[i][j] = 0.0;
		  }
	      }
	    
	    if(v==0)
	      {
		
		a0 = ( x0*r[0] + x1*r[1] + x2*r[2] + x3*r[3])/chi;
		a3 = (-x0*r[3] + x3*r[0] - x2*r[1] + x1*r[2])/chi;
		a1 = (-x0*r[1] - x3*r[2] + x2*r[3] + x1*r[0])/chi;
		a2 = (-x0*r[2] + x3*r[1] + x2*r[0] - x1*r[3])/chi;

		//cout<<"\t Modulus of 4-vector a = "<<(a0*a0 + a1*a1 + a2*a2 + a3*a3)<<endl;
		//cout<<"RMatrix for _v=1 is \n";
		AA[0][0] = TComplex(a0,a3);
		AA[0][1] = TComplex(a2,a1);
		AA[1][0] = TComplex(-a2,a1);
		AA[1][1] = TComplex(a0,-a3);
		AA[0][2] = TComplex(0.0,0.0);
		AA[1][2] = TComplex(0.0,0.0);
		AA[2][2] = TComplex(1.0,0.0);
		AA[2][0] = TComplex(0.0,0.0);
		AA[2][1] = TComplex(0.0,0.0);

		for(int i = 0; i<3; i++)
		  {
		    for(int j = 0; j<3; j++)
		      { 
			for(int m = 0; m<3; m++)
			  {
			    RNew[i][j] = RNew[i][j] + AA[i][m]*ROld[m][j];
			    UNew[i][j] = UNew[i][j] + AA[i][m]*UOld[m][j];
			  }//cout<<RNew[i][j]<<"\t";
		      }
		    //cout<<"\n";
		  }
	      }
		    
	    if(v==1)
              {
		a0 = ( x0*r[0] + x1*r[1] + x2*r[2] + x3*r[3])/chi;
                a3 = (-x0*r[3] + x3*r[0] - x2*r[1] + x1*r[2])/chi;
		a1 = (-x0*r[1] - x3*r[2] + x2*r[3] + x1*r[0])/chi;
		a2 = (-x0*r[2] + x3*r[1] + x2*r[0] - x1*r[3])/chi;

		//cout<<"RMatrix for _v=2_ is \n";
                AA[0][0] = TComplex(1.0,0.0);
                AA[0][1] = TComplex(0.0,0.0);
                AA[0][2] = TComplex(0.0,0.0);
                AA[1][0] = TComplex(0.0,0.0);
                AA[2][0] = TComplex(0.0,0.0);
                AA[1][1] = TComplex(a0,a3);
                AA[1][2] = TComplex(a2,a1);
                AA[2][1] = TComplex(-a2,a1);
                AA[2][2] = TComplex(a0,-a3);

                for(int i = 0; i<3; i++)
                  {
                    for(int j = 0; j<3; j++)
                      {
                        for(int m = 0; m<3; m++)
                          {
                            RNew[i][j] = RNew[i][j] + AA[i][m]*ROld[m][j];
			    UNew[i][j] = UNew[i][j] + AA[i][m]*UOld[m][j];
                          }//cout<<RNew[i][j]<<"\t";
                      }
		    //cout<<"\n";
                  }

	      }
	    
	    if(v==2)
              {
		a0 = ( x0*r[0] + x1*r[1] + x2*r[2] + x3*r[3])/chi;
                a3 = (-x0*r[3] + x3*r[0] - x2*r[1] + x1*r[2])/chi;
		a1 = (-x0*r[1] - x3*r[2] + x2*r[3] + x1*r[0])/chi;
		a2 = (-x0*r[2] + x3*r[1] + x2*r[0] - x1*r[3])/chi;

                AA[1][1] = TComplex(1.0,0.0);
                AA[0][1] = TComplex(0.0,0.0);
                AA[1][0] = TComplex(0.0,0.0);
                AA[1][2] = TComplex(0.0,0.0);
                AA[2][1] = TComplex(0.0,0.0);
                AA[0][0] = TComplex(a0,a3);
                AA[0][2] = TComplex(a2,a1);
                AA[2][0] = TComplex(-a2,a1);
                AA[2][2] = TComplex(a0,-a3);

                for(int i = 0; i<3; i++)
                  {
                    for(int j = 0; j<3; j++)
                      {
                        for(int m = 0; m<3; m++)
                          {
                            RNew[i][j] = RNew[i][j] + AA[i][m]*ROld[m][j];
			    UNew[i][j] = UNew[i][j] + AA[i][m]*UOld[m][j];
                          }
                      }

                  }
	      }

	    //Determinant of RNew
	    //cout<<"\n Determinant of UNew is "<<UNew[0][0]*( UNew[1][1]*UNew[2][2] - UNew[1][2]*UNew[2][1] ) - UNew[0][1]*( UNew[1][0]*UNew[2][2] - UNew[1][2]*UNew[2][0] ) + UNew[0][2]*( UNew[1][0]*UNew[2][1] - UNew[1][1]*UNew[2][0] )<<endl;


    

    }

  //cout<<"Checking Unitarity of matrix \n";

  TComplex UCheck[3][3]={0.0};
 
  for(int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
	{
	  for(int m = 0; m<3; m++)
	    {
	      UCheck[i][j] = UCheck[i][j] + UNew[i][m]*TComplex::Conjugate(UNew[j][m]);
	    }//cout<<UCheck[i][j]<<"\t";
	}//cout<<"\n";
    }  

  //cout<<"\n Determinant of UNew is "<<UNew[0][0]*( UNew[1][1]*UNew[2][2] - UNew[1][2]*UNew[2][1] ) - UNew[0][1]*( UNew[1][0]*UNew[2][2] - UNew[1][2]*UNew[2][0] ) + UNew[0][2]*( UNew[1][0]*UNew[2][1] - UNew[1][1]*UNew[2][0] )<<endl;   

  //cout<<"\n New Link matrix is \n";

  for(Int_t i=0; i<3; i++)
    {
      for(Int_t j =0; j<3; j++)
	{
	  LinksContent[CurrentLink][i][j] = UNew[i][j] ;		 
	  //cout<<UNew[i][j]<<"\t";
	}//cout<<"\n";
    }
  

  //cout<<"Update is done \n";
              
	      
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
	
	  for(Int_t i=0; i<3; i++)
	    {
	      for(Int_t m=0; m<3; m++)
		{
		  for(Int_t n=0; n<3; n++)
		    {
		      for(Int_t p=0; p<3; p++)
			{
			  U1 = LinksContent[n1][i][m];

			  U2 = LinksContent[n2][m][n];

			  U3 = TComplex::Conjugate(LinksContent[n3][p][n]);//(Dagger)

			  U4 = TComplex::Conjugate(LinksContent[n4][i][p]);//(Dagger)

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

  TComplex PlaquetteMatrix[3][3],FieldStrengthMatrix[3][3],U1,U2,U3,U4,TraceFSquare;

  U1 = TComplex(0.0,0.0); U2 = TComplex(0.0,0.0);
  U3 = TComplex(0.0,0.0); U4 = TComplex(0.0,0.0);TraceFSquare = TComplex(0.0,0.0);

  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
        {
	  PlaquetteMatrix[i][j] = TComplex(0.0,0.0);
	  FieldStrengthMatrix[i][j] = TComplex(0.0,0.0);
	}
    }
  //PlaquetteMatrix[0][0] = TComplex(0.0,0.0);PlaquetteMatrix[0][1] = TComplex(0.0,0.0);
  //PlaquetteMatrix[1][0] = TComplex(0.0,0.0);PlaquetteMatrix[1][1] = TComplex(0.0,0.0);

  //FieldStrengthMatrix[0][0]= TComplex(0.0,0.0);FieldStrengthMatrix[0][1]= TComplex(0.0,0.0);
  //FieldStrengthMatrix[1][0]= TComplex(0.0,0.0);FieldStrengthMatrix[1][1]= TComplex(0.0,0.0);
  

  //Calculating PlaquetteMatrix, Pmunu = U1*U2*U3*U4
  //cout<<"Plaquette Content\n"<<endl;
  //cout<<LinksContent[n1][0]<<"\t"<<LinksContent[n1][1]<<"\t"<<LinksContent[n1][2]<<"\t"<<LinksContent[n1][3]<<endl;
  //cout<<LinksContent[n2][0]<<"\t"<<LinksContent[n2][1]<<"\t"<<LinksContent[n2][2]<<"\t"<<LinksContent[n2][3]<<endl;
  //cout<<LinksContent[n3][0]<<"\t"<<LinksContent[n3][1]<<"\t"<<LinksContent[n3][2]<<"\t"<<LinksContent[n3][3]<<endl;
  //cout<<LinksContent[n4][0]<<"\t"<<LinksContent[n4][1]<<"\t"<<LinksContent[n4][2]<<"\t"<<LinksContent[n4][3]<<endl;
  //cout<<"Absurd LinNumber = "<<n4<<endl;

  //cout<<"Plaquette Matrix is = \n";
  
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
        {
	  for(int m=0; m<3; m++)
	    {
	      for(int n=0; n<3; n++)
		{
		  for(int p=0; p<3; p++)
		    {
		      U1 = LinksContent[n1][i][m];
		      
		      U2 = LinksContent[n2][m][n];
		      
		      U3 = TComplex::Conjugate(LinksContent[n3][p][n]);
		      
		      U4 = TComplex::Conjugate(LinksContent[n4][j][p]);
		      
		      PlaquetteMatrix[i][j] = PlaquetteMatrix[i][j] + U1*U2*U3*U4;
		      //cout<<"PlaquettMatrix[i][j]"<<PlaquetteMatrix[i][j]<<endl;
		    }
		}
	    }//cout<<PlaquetteMatrix[i][j]<<"\t";
	} //cout<<"\n";
    }   
  
  //Calculating Field-StrengthMatrix
  
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
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
  for(int i=0; i<3; i++)
    {
      for(int m=0; m<3; m++)	
        {
	  TraceFSquare = TraceFSquare + FieldStrengthMatrix[i][m]*FieldStrengthMatrix[m][i];
	}
    }
  
  //cout<<"TraceFSquare is "<<TraceFSquare<<endl;
  //cout<<"TraceFSquare*(-4*g^2*a^4) is "<<TraceFSquare*(-4.0*pow(BareCouplingConstant,2.0)*pow(LatticeSpacing,4.0))<<endl;
  return (TraceFSquare.Re());
}


double LatticeCalculation::TraceFieldStrengthCrossedTerm(int* N1,int* N2,int* N3, int* N4)
{
  //F31==0
  //F41==1
  //F32==2
  //F42==3
  //       a2
  //        _
  //     a1|_|a3
  //       a4
  //
  //
  //N1={a1,a2,a3,a4},  N2={b1,b2,b3,b4},  N3={c1..c4}, N4={d1...d4}
  int n1[4],n2[4],n3[4],n4[4];
  TComplex PlaquetteMatrix[4][3][3],FieldStrengthMatrix[4][3][3],U1,U2,U3,U4,TraceF31F41,TraceF32F42;

  U1=TComplex(0.0,0.0); U2=TComplex(0.0,0.0);
  U3=TComplex(0.0,0.0); U4=TComplex(0.0,0.0);TraceF31F41=TComplex(0.0,0.0);TraceF32F42=TComplex(0.0,0.0);

  for(int i=0; i<4; i++)
    {for(int j=0; j<3; j++)
        {for(int k=0; k<3; k++)
            {PlaquetteMatrix[i][j][k]=TComplex(0.0,0.0); FieldStrengthMatrix[i][j][k]=TComplex(0.0,0.0);}
        }
    }
  for(int k=0; k<4; k++)
    {
      n1[k]=N1[k];n2[k]=N2[k];n3[k]=N3[k];n4[k]=N4[k];
      //cout<<"Plaquette = "<<n1[k]<<"\t"<<n2[k]<<"\t"<<n3[k]<<"\t"<<n4[k]<<endl;
    }


  for(int k=0; k<4; k++)
    {
      for(int i=0; i<3; i++)
        {
          for(int j=0; j<3; j++)
            {
              for(int m=0; m<3; m++)
                {
                  for(int n=0; n<3; n++)
                    {
                      for(int p=0; p<3; p++)
                        {
                          U1 = LinksContent[n1[k]][i][m];

                          U2 = LinksContent[n2[k]][m][n];

                          U3 = TComplex::Conjugate(LinksContent[n3[k]][p][n]);//conjugate

			  U4 = TComplex::Conjugate(LinksContent[n4[k]][j][p]);//Conjugate

                          PlaquetteMatrix[k][i][j] = PlaquetteMatrix[k][i][j] + U1*U2*U3*U4;


			}
                    }
                }
            }
        }
    }

  //Calculating Field-StrengthMatrix
  for(int k=0; k<4; k++)
    {
      for(int i=0; i<3; i++)
        {
          for(int j=0; j<3; j++)
            {

              FieldStrengthMatrix[k][i][j] = (PlaquetteMatrix[k][i][j] - (TComplex::Conjugate(PlaquetteMatrix[k][j][i])))/(2.0*IYOTA*BareCouplingConstant*LatticeSpacing*LatticeSpacing);


	      //cout<<"IYOTA,g,a"<<IYOTA<<"\t"<<BareCouplingConstant<<"\t"<<LatticeSpacing<<endl;
              //cout<<"Fmunu["<<i<<"]["<<j<<"]"<<FieldStrengthMatrix[i][j]<<endl;;
              //cout<<"PlaquettMatrix["<<i<<"]["<<j<<"]"<<PlaquetteMatrix[i][j]<<endl;
              //cout<<"Fmunu["<<i<<"]["<<j<<"]"<<FieldStrengthMatrix[i][j]<<endl;

            }
        }
    }
  //cout<<"2*i*g*a^2 = "<<2.0*IYOTA*BareCouplingConstant*LatticeSpacing*LatticeSpacing<<endl;
  
  
  //Calculating Trace of SquaredFieldStrengthMatrix
  for(int i=0; i<3; i++)
    {
      for(int m=0; m<3; m++)
        {
          TraceF31F41 = TraceF31F41 + FieldStrengthMatrix[0][i][m]*FieldStrengthMatrix[1][m][i];
        }
    }

  for(int i=0; i<3; i++)
    {
      for(int m=0; m<3; m++)
        {
          TraceF32F42 = TraceF32F42 + FieldStrengthMatrix[2][i][m]*FieldStrengthMatrix[3][m][i];
        }
    }

  //cout<<"TraceFSquare is "<<TraceF31F41+TraceF32F42<<endl;
  //cout<<"TraceFSquare*(-4*g^2*a^4) is "<<TraceFSquare*(-4.0*pow(BareCouplingConstant,2.0)*pow(LatticeSpacing,4.0))<<endl;
  return (TraceF31F41.Re() + TraceF32F42.Re());
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

void AverageCorrelatorF31F41PlusF32F42SU3(int nt,int nx,double Beta,int TC)
{ int StartTime = time(NULL);

  //const double PI = 3.1415926; //Value of pi
  //const int NIteration = 10;
  //srand(time(NULL));
  LatticeCalculation LC;
 
  
  //..............Input Parameter.....Lattice.....
  LC.PI = 3.1415926;
  LC.Nt = nt; //Discretizing each dimension in N-points 
  LC.Nx = nx;
  LC.Ny = nx;
  LC.Nz = nx;
  LC.LatticePoints = LC.Nt*LC.Nx*LC.Ny*LC.Nz; //                            
  LC.TotalLinks = 4*LC.Nt*LC.Nx*LC.Ny*LC.Nz; // 4 space-time dimension      
  LC.NLinksX = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in X-direction        
  LC.NLinksY = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Y-direction        
  LC.NLinksZ = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Z-direction        
  LC.NLinksT = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in T-direction        
  LC.Run = 50;                //Number of iteration to get stable gauge field configuration                                     
  LC.TotalConfiguration = TC; //Number of gauge field configuartion over which average of (F3iF3i-F4iF4i) is computed         
  LC.TotalTemperatureRun = 1; //Number of different values of Temperature to see the dependence of correlation function         
  LC.lambdaLattice = 5.3; //(in MeV) A parameter that goes in LatticeSpacing formular using renormalization group equation 
  
  LC.DefineDynamicArrays();
  LC.DefineComplexMatrix();
  Int_t Vacuum =0;
  Int_t Thermal =1;
  //..............End of Input Parameter........Lattice
  cout<<"Program started"<<endl;
  cout<<"Dimension of Lattice used is  \t" <<LC.Nt<<"*"<<LC.Nx<<"*"<<LC.Ny<<"*"<<LC.Nz<<", So Total Links are "<<LC.TotalLinks<<endl;
  cout<<"\t"<< LC.LatticeSize[0]<<","<<LC.LatticeSize[1]<<","<<LC.LatticeSize[2]<<","<<LC.LatticeSize[3]<<endl;
  
  Double_t Normalization = 2*3*3.0;
  ofstream f1;
  //f1.open("Coordinate.txt",ios::out);
  //f1.open("Average.txt",ios::out);
  //f1.open("Plaquette_wrt_beta_N6_I20.txt",ios::app);
  Char_t Var[200]; int INTBeta = Beta*100;
  sprintf(Var,"10RunSU3/SU3Correlator_F31F41_PlusF32F42__wrt_betaTimes100_%i_Nt%i_Ns%i_I%i.txt",INTBeta,LC.Nt,LC.Nx,LC.TotalConfiguration);
    f1.open(Var,ios::app);    

//Labeling each link with its 4 space-time coordinate and a direction
  LC.LabelLinksCoordinate();
  
  //cout<<"Line 807 This is working"<<endl;
  
     
  int TemperatureLoop =0;
  int CountPlane31;
  double SumFieldStrengthF3iF4i;//Some Over all plaquette for given configuration
  double SumScaledF3iF4iPerPlaquette;
  double *ScaledAverageFieldStrengthCorrelation = new double[LC.TotalTemperatureRun];
  double *TEMPERATURE = new double[LC.TotalTemperatureRun];
  double *BETA = new double[LC.TotalTemperatureRun];
 
  double *TraceScaledF3iF4iPerPlaquette = new double[LC.TotalConfiguration];
  double AverageScaledTraceF3iF4i;
  double SumDeviationSquareF3iF4i;

  double ErrorF3iF4i;
  double *ScaledAvgFSCorrelationError = new double[LC.TotalTemperatureRun];

  //double BETAINT[11]={1.87,1.9,2.0,2.1,2.3,2.4,2.5,2.6,2.645,2.7,2.2};
  //double ThermalTemperatureNt3[11]={102.439,110.312,141.297,181.176,298.708,384.022,494.064, 636.067, 712.799, 819.397, 232.5 };


  for(TemperatureLoop = 0; TemperatureLoop<LC.TotalTemperatureRun; TemperatureLoop++)
    {
      LC.beta = Beta;//BETAINT[INTBeta];//1.80 + 0.1*TemperatureLoop;
      LC.BareCouplingConstant =1.0;// TMath::Sqrt(6.0/LC.beta);
      LC.LatticeSpacing =1.0;// pow(11.0*TMath::Power(LC.BareCouplingConstant/LC.PI,2.0)/24.0,-51.0/121.0)*TMath::Exp(-12.0*TMath::Power(LC.PI/LC.BareCouplingConstant,2.0)/11.0)/(LC.lambdaLattice);
      LC.temperature =1.0;// 1/(LC.Nt*LC.LatticeSpacing);
      TEMPERATURE[TemperatureLoop] = LC.temperature;
      BETA[TemperatureLoop] = LC.beta;

      cout<<"Coupling is "<<LC.BareCouplingConstant<<endl;
      cout<<"Beta (4/g^2)= "<< LC.beta<<endl;
      cout<<"LatticeSpacing "<<LC.LatticeSpacing<<endl;
      cout<<"Temperature is "<<LC.temperature<<endl;
      cout<<"(nt)^4 "<<TMath::Power(LC.Nt,4.0)<<endl;
      cout<<"1/T^4 "<<TMath::Power(1/LC.temperature,4.0)<<endl;;
      cout<<"1/(4*g^2*a^4) "<<1/TMath::Power(2*LC.BareCouplingConstant*TMath::Power(LC.LatticeSpacing,2),2.0)<<endl;
      cout<<"1/(4*g^2*a^4*T^4) "<<1/TMath::Power(2*LC.BareCouplingConstant*TMath::Power(LC.LatticeSpacing*LC.temperature,2.0),2.0)<<endl;

      SumFieldStrengthF3iF4i  = 0.0;
      SumScaledF3iF4iPerPlaquette=0.0;
      AverageScaledTraceF3iF4i=0.0;
      SumDeviationSquareF3iF4i = 0.0;
      ScaledAverageFieldStrengthCorrelation[TemperatureLoop]=0.0;
      ScaledAvgFSCorrelationError[TemperatureLoop]=0.0;  


      //Setting all link variables to Unity
      for(Int_t i = 0; i<LC.TotalLinks; i++)
	{ 
	  for(Int_t m = 0; m<3; m++)
	    {
	      for(Int_t n = 0; n<3; n++)
		{      
		  LC.LinksContent[i][m][n]=TComplex(0.0,0.0);
		}
	    }
      
	  for(Int_t m = 0; m<3; m++)
	    {
	      LC.LinksContent[i][m][m]=TComplex(1.0,0.0);      //
	    }
	  
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

          SumFieldStrengthF3iF4i = 0.0;
          CountPlane31 = 0;
          TraceScaledF3iF4iPerPlaquette[Iter]=0.0;


          //Touching each internal link with heat bath
          for(Int_t ic = 0 ; ic< LC.TotalLinks ; ic++)
            {
              //cout<<"Processing Link"<<ic<<endl;
              LC.UpdateLink(ic); //Apply heat bath algorithm to each internal link and apply periodic boundary condition
            }

          //cout<<"Going to evaluate the Field Strength matric for generated configuration"<<endl;
          //Evaluate Trace of Field-Strength-Field-Strength Matrix

          Int_t xt[4] ={0,0,0,0};
          Int_t N1[4]={0},N2[4]={0},N3[4]={0},N4[4]={0};
          Int_t d1=0,d2=0,p1=0,p2=0;


          for(Int_t i=0; i<LC.LatticeSize[0]; i++)  //time
            {
              for(Int_t j=0; j<LC.LatticeSize[1]; j++)//x-axis
                {
                  for(Int_t k=0; k<LC.LatticeSize[2]; k++)  //y-axis
                    {
                      for(Int_t m=0; m<LC.LatticeSize[3]; m++) //z-axis
			{
			  //cout<<"Going to calculate Fmunu using Plaquette \n\n";
			
			  for(Int_t id = 0; id <4; id++)
			    {
			      
			      if(id==0){d1=1; d2=3; p1=0; p2=2;}
			      if(id==1){d1=1; d2=0; p1=2; p2=3;}
			      if(id==2){d1=2; d2=3; p1=0; p2=1;}
			      if(id==3){d1=2; d2=0; p1=1; p2=3;}
			      
			      xt[0]=i;
			      xt[1]=j;
			      xt[2]=k;
			      xt[3]=m;
			      //cout<<"t x y z  d2(v) d1(h) = "<<i<<"\t"<<j<<"\t"<<k<<"\t"<<m<<"\t"<<d2<<"\t"<<d1<<endl;
			      
			      
			      Int_t Link4 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d1];
			      
			      Int_t Link1 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d2];
			      
			      Int_t temp1 = xt[d2];
			      xt[d2]=xt[d2]+1;
			      
			      if(xt[d2] == LC.LatticeSize[d2])
				{ xt[d2]=0;
				}
			      Int_t Link2 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d1];
			      
			      xt[d2] = temp1;
			      xt[d1] = xt[d1]+1;
			      
			      if( xt[d1] == LC.LatticeSize[d1])
				{xt[d1]=0;
				}
			      Int_t Link3 = LC.LinksNumber[xt[0]][xt[1]][xt[2]][xt[3]][d2];
			      //cout<<"We have formed the plaquette square \n\n";
			      
			      N1[id]=Link1;
			      N2[id]=Link2;
			      N3[id]=Link3;
			      N4[id]=Link4;
			    }
			  
			  
			  CountPlane31 = CountPlane31 + 1;
			  SumFieldStrengthF3iF4i = SumFieldStrengthF3iF4i + LC.TraceFieldStrengthCrossedTerm(N1,N2,N3,N4);
			  //cout<<"In Plane (3,2):  summation TraceUmunu^2 = "<<SumFieldStrengthCorrelation32*(4.0)*pow(LC.BareCouplingConstant*LC.LatticeSpacing*LC.LatticeSpacing,2.0)/(CountPlane32)<<endl;
			  
			  //cout<<"Summation (FmunuFmnu)"<<SumFieldStrengthCorrelation<<endl;
			}
		    }
		}
	    }


  if(Thermal ==1)
    {
      TraceScaledF3iF4iPerPlaquette[Iter] = SumFieldStrengthF3iF4i/(TMath::Power(LC.temperature,4.0)*CountPlane31*Normalization);
    }
  if(Vacuum ==1)
    {
      //TraceScaledF3iF4iPerPlaquette[Iter] = SumFieldStrengthF3iF4i/(TMath::Power(ThermalTemperatureNt3[INTBeta],4.0)*CountPlane31);
    }

  SumScaledF3iF4iPerPlaquette = SumScaledF3iF4iPerPlaquette + (TraceScaledF3iF4iPerPlaquette[Iter]);

  cout<<" Trace(F3iF4i)/N31 = "<<SumFieldStrengthF3iF4i/CountPlane31<<endl;
  cout<<" Trace(F3iF4i)/(N31*T^4*Normalization)  "<<TraceScaledF3iF4iPerPlaquette[Iter]<<endl;
  cout<<"N31 Plane= "<<CountPlane31<<endl;

  cout<<"ScaledAverageF3iF4i So far = "<<SumScaledF3iF4iPerPlaquette/(Iter+1)<<endl;

	}

AverageScaledTraceF3iF4i = SumScaledF3iF4iPerPlaquette/LC.TotalConfiguration;

for(int i=0; i<LC.TotalConfiguration; i++)
  {
    SumDeviationSquareF3iF4i =   SumDeviationSquareF3iF4i + TMath::Power(AverageScaledTraceF3iF4i - TraceScaledF3iF4iPerPlaquette[i],2.0);
  }

ErrorF3iF4i = TMath::Sqrt(SumDeviationSquareF3iF4i/(LC.TotalConfiguration-1.0))/TMath::Sqrt(LC.TotalConfiguration);

ScaledAverageFieldStrengthCorrelation[TemperatureLoop] = AverageScaledTraceF3iF4i;

ScaledAvgFSCorrelationError[TemperatureLoop] =  ErrorF3iF4i;


//cout<<"Number of Plaquette in (3,1) plane = "<<CountPlane31<<endl;
//cout<<LC.temperature<<"\t"<<ScaledAverageFieldStrengthCorrelation[TemperatureLoop]<<endl;
if(Thermal == 1)
  {
    f1<<TEMPERATURE[TemperatureLoop]<<"\t"<<LC.beta<<"\t"<<ScaledAverageFieldStrengthCorrelation[TemperatureLoop]<<"\t"<<ScaledAvgFSCorrelationError[TemperatureLoop]<<"\t"<<ScaledAvgFSCorrelationError[TemperatureLoop]*100.0/ScaledAverageFieldStrengthCorrelation[TemperatureLoop]<<endl;
  }

 if(Vacuum == 1)
   {
          //f1<<ThermalTemperatureNt3[INTBeta]<<"\t"<<ScaledAverageFieldStrengthCorrelation[TemperatureLoop]<<"\t"<<ScaledAvgFSCorrelationError[TemperatureLoop]<<"\t"<<ScaledAvgFSCorrelationError[TemperatureLoop]*100.0/ScaledAverageFieldStrengthCorrelation[TemperatureLoop]<<"\t"<<LC.beta<<"\t"<<TEMPERATURE[TemperatureLoop]<<endl;
   }
}        
  f1.close();  
  int EndTime = time(NULL);
  int Hour = (EndTime-StartTime)/3600;
  int Minute = ((EndTime-StartTime)/60)-Hour*60;
  int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
  cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
      
}
