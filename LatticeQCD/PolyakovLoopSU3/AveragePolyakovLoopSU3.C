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
  int* LatticeSize;  //  LatticeSize[4]={Nt,Nx,Ny,Nz};
  double beta;

  MTRand ObjectMTRand;//Define a object to access rand() function defined inside class MTRand in file MersenneTwister.h

  TComplex*** LinksContent;   //LinksContent[TotalLinks][3][3];
  int** LinksCoordinate;   //LinksCoordinate[TotalLinks][4];
  int* LinksDirection;     //LinksDirection[TotalLinks];
  int***** LinksNumber; //LinkNumber[t,x,y,z,d]

  int CurrentLink, Iteration;
   
  TComplex ID[2][2],Pauli[3][2][2],CPauli[3][2][2],IYOTA,CIYOTA; //Complex Matrix
  
  TComplex SumProduct3Us[3][3];  //Complex matrix Staple U2*U3*U4 (It is staple)


  void DefineDynamicArrays();
  void LabelLinksCoordinate();
  Int_t OutSideBoundary(Int_t,Int_t ,Int_t, Int_t,Int_t);
  void UpdateLink(Int_t);
  Int_t SearchLinkNumber(Int_t, Int_t, Int_t, Int_t, Int_t);
  void DefineComplexMatrix();
  double TracePlaquette(int,int,int,int);
  TComplex TracePolyakov(int,int,int);
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

void LatticeCalculation::UpdateLink(int  GG)
{
  int StapleCoordinate[3][6][5]={0};
  int CoordinateLinkI[5]={0};
  int LinkNumberStaple[6][3]={0}; //6 loop; each with three links

  CurrentLink = GG;
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
	      LinkNumberStaple[count][i2] =  SearchLinkNumber(StapleCoordinate[p][ra][0], StapleCoordinate[p][ra][1], StapleCoordinate[p][ra][2],  StapleCoordinate[p][ra][3],  StapleCoordinate[p][ra][4] );
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
  TComplex UNew[3][3], UOld[3][3];
  TComplex ROld[3][3], RNew[3][3];
  TComplex AA[3][3];
  for(int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
	{
	  UNew[i][j]=TComplex(0.0,0.0);UOld[i][j]=TComplex(0.0,0.0);
	  RNew[i][j]=TComplex(0.0,0.0);ROld[i][j]=TComplex(0.0,0.0);
	  AA[i][j]=TComplex(0.0,0.0);
	}
    }

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
      r[0] = ((RNew[0][0].Re()) + (RNew[1][1].Re()))/2.0;
      r[3] = ((RNew[0][0].Im()) - (RNew[1][1].Im()))/2.0 ;
      r[1] = ((RNew[0][1].Im()) + (RNew[1][0].Im()))/2.0;
      r[2] = ((RNew[0][1].Re()) - (RNew[1][0].Re()))/2.0 ;
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

  //cout<<"In beta: The factor chi, and net factor = "<<chi<<","<<beta*chi*2.0/3.0<<endl;

	  //Create pdf and generate x0 x1 x2 x3 with apporpiate weight 
	  
          Double_t x0, x1, x2, x3,Theta, Phi;
	  
	  /////////////////////............Differential Sampling for a0.......////////////////////////
	  
	  Double_t x = 1.0;
	  Double_t Zmin, Zmax, Z, dpdz, dpdzmax=0.0, dz;
	  Double_t r1 , r2;
	  Int_t Success =0;
	  Zmin = TMath::Exp(-x*beta*chi*2.0/3.0);
	  Zmax = TMath::Exp(x*beta*chi*2.0/3.0);
	  /*
	  dz = (Zmax-Zmin)/(10000000.0-1.0);                                                                                                                
	  for(Int_t i=0; i< 10000000; i++)                                                                                                                  
	    {Z = Zmin + (i*dz);dpdz = TMath::Sqrt(1.0 - (TMath::Power(3.0/(2.0*beta*chi),2.0)*TMath::Power(log(Z),2.0)) );
	      if( dpdzmax < dpdz ){dpdzmax = dpdz;}} cout<<"dpdzmax is "<<dpdzmax<<",chi is "<<chi<<endl;
	  */
	  dpdzmax =1.0;
	  Int_t SearchCount = 0;//cout<<"Search for x0 started "<<endl;                                                                                     
	  while(Success==0)
	    {
	      Double_t Dummy1 = ObjectMTRand.rand();
	      Double_t Dummy2 = ObjectMTRand.rand();
	      r1 = Dummy1;r2 = Dummy2;
	      Z = Zmin + (Zmax - Zmin)*r1;
	      dpdz = TMath::Sqrt(1.0 - (TMath::Power(3.0/(2.0*beta*chi),2.0)*TMath::Power(log(Z),2.0)) );

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


int LatticeCalculation::SearchLinkNumber(int t, int x, int y, int z, int d)
{
  //cout<<"\t\t  ..... Search Function started"<<endl;
  //cout<<"\t \t .... txyzd = "<<t<<","<<x<<","<<y<<","<<z<<","<<d<<endl; 
  Int_t found =0; //1 means found
  Int_t Link = 0;
  Int_t ra = 0;
  //If Lattice point is outside the grid, perform a periodic transformation 
  if( t<0 || t > Nt-1 )
    {
      //cout<<"Alert Link is outside the grid, using periodic transformation to bring inside"<<endl; 
      Int_t Temp=t; t = Temp%Nt;
    }
  if( x<0 || x > Nx-1 )
    {
      //cout<<"Alert Link is outside the grid, using periodic transformation to bring inside"<<endl;   
      Int_t Temp=x; x =Temp%Nx;
    }
  if( y<0 || y > Ny-1 )
    {
      //cout<<"Alert Link is outside the grid, using periodic transformation to bring inside"<<endl;     
      Int_t Temp=y; y =Temp%Ny;
    }
  if( z<0 || z > Nz-1 )
    {
      //cout<<"Alert Link is outside the grid, using periodic transformation to bring inside"<<endl;    
      Int_t Temp=z; z =Temp%Nz;
    }

  /*
  while(ra< TotalLinks && found == 0) 
    {              
      if(d == LinksDirection[ra])       {                      
	if(t == LinksCoordinate[ra][0] && x == LinksCoordinate[ra][1] && y == LinksCoordinate[ra][2] && z == LinksCoordinate[ra][3])          
	  {                                                                                                                                   
	    found = 1;              Link = ra;
	  }     else {ra = ra+1;}                                                                                                             
      }      else {ra = ra + 1;}                                                                                                              
    } 
  */
  Link = LinksNumber[t][x][y][z][d];
  if(Link > -1 && Link < TotalLinks)
    {
      //cout<<"\t\t .... End using searchFunction with Link number = "<<Link<<endl;
      return Link;
    }
  else{
    cout<<"\t Alert Link not found"<<endl;
    return -1;
  }

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

	  cout<<"Trace U = "<<TraceU<<endl;
  //if(TraceU.Re() < 0.0){TraceU=TComplex(0.0,0.0);}
  return (TraceU.Re());
}

TComplex LatticeCalculation::TracePolyakov(int x, int y, int z)
{

  TComplex BNew[3][3], BOld[3][3],TraceU;
  Int_t LNum=0, t1 =0; 
  //cout<<"t = "<<x<<endl;
  TraceU  = TComplex(0.0,0.0);
  
  LNum = SearchLinkNumber(0,x,y,z,0);
  //cout<<"LNUm at t0 = "<<LNum<<endl;

  for(int i=0 ; i<3;i++)
    {
      for(int j=0; j<3; j++)
        {
	  BOld[i][j] = TComplex(0.0,0.0);
	  BNew[i][j] = LinksContent[LNum][i][j];
	  //cout<<BNew[i][j]<<"\t";
	}
      //cout<<"\n";
    }

  for(int m=1;m<Nt; m++)
    {
      
      for(int i1=0 ; i1<3;i1++)
	{
	  for(int j1=0; j1<3; j1++)
	    {
	      BOld[i1][j1] = BNew[i1][j1];
	      BNew[i1][j1] = TComplex(0.0,0.0);
	    }
	}

      t1 = 0 + m;
      LNum = SearchLinkNumber(t1, x, y, z, 0);
      //cout<<" LNum for Next t \t"<<LNum<<endl;
      for(int i=0 ; i<3;i++)
	{
	  for(int j=0; j<3; j++)
	    {
	      for(int k=0; k<3; k++)
		{
		  BNew[i][j] = BNew[i][j] + BOld[i][k]*LinksContent[LNum][k][j];
		}
	      //cout<<BNew[i][j]<<"\t";
	    }//cout<<"\n";
	}
    }
  
  //Calculating Trace
  
  for(int i=0; i<3; i++)
    {
	  TraceU = TraceU + BNew[i][i];

    }
  //TracePolyakovLine = TraceU;
  //cout<<"PolyakovLoop Trace is "<<TraceU<<endl;
  //cout<<"Abs PolyakovLoop Trace is "<<TMath::Sqrt(TMath::Power(TraceU.Re(),2.0) + TMath::Power(TraceU.Im(),2.0))<<endl;
  //return TraceU.Re();
  return TraceU;
  //return ( TMath::Sqrt(TMath::Power(TraceU.Re(),2.0) + TMath::Power(TraceU.Im(),2.0)) );

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
int main(int argc, char* argv[])
//void AveragePolyakovLoopSU3(int nt,int nx,double Beta)
{
  int StartTime = time(NULL);
  int nt=atoi(argv[1]), nx=atoi(argv[2]), totalrun=atoi(argv[4]); double Beta=atof(argv[3])/1.0;
  //const double PI = 3.1415926; //Value of pi
  //const int NIteration = 10;
  //srand(time(NULL));
  LatticeCalculation LC;

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
  LC.Run = totalrun;                //Number of iteration to get stable gauge field configuration  
  

  LC.DefineDynamicArrays();
  LC.DefineComplexMatrix();
  //.............End of Input Parameter.....Lattice

  LC.beta = Beta;
  //LC.beta = 0.4 + Beta*0.4;
  cout<<"Program started"<<endl;
  cout<<"Dimension of Lattice used is  \t" <<LC.Nt<<"*"<<LC.Nx<<"*"<<LC.Ny<<"*"<<LC.Nz<<", So Total Links are "<<LC.TotalLinks<<endl;
  cout<<"Total Plaquette should be 6*Nt*Nx*Ny*Nz "<<6*LC.Nt*LC.Nx*LC.Ny*LC.Nz<<endl;
  cout<<"Total Polyakov loop should be Nx*Ny*Nz "<<LC.Nx*LC.Ny*LC.Nz<<endl;
  Char_t FileName1[5000],FileName2[5000];
  sprintf(FileName2,"/wsu/home/fy/fy41/fy4125/Lattice/PolyakovLoopSU3/OutPutFile/APolyakovEachConfiguration_wrt_I_Nt_%i_Ns_%i_betaTimes100_%i.txt",LC.Nt,LC.Nx,int(Beta*100));
  sprintf(FileName1,"/wsu/home/fy/fy41/fy4125/Lattice/PolyakovLoopSU3/OutPutFile/APolyakov_wrt_betaTimes100_%i_Nt_%i_Ns_%i_I500.txt",int(Beta*100),LC.Nt,LC.Nx);
  ofstream f1,f2;
  f1.open(FileName1,ios::out);
  f2.open(FileName2,ios::out);
  //f1.open("AveragePlaquetteForBeta3_2.txt",ios::out);
  //f1.open("SU3/DataFilesAndPlots/Plaquette_wrt_beta_N6_I40.txt",ios::app);
  //f1.open("Plaquette_wrt_beta_N3_12_12_12_I30.txt",ios::app);
  f1<<"#Beta = "<<Beta<<"\t AveragePolyakovLoop = "<<"\tErrorAPL "<<endl;
  f2<<"#Iteration \t"<<"AveragePolyakovLoop \t"<<"\t ErrorAvgPolyakov"<<endl;

//OBLabeling each link with its 4 space-time coordinate and a direction
  LC.LabelLinksCoordinate();     
    
//Setting all link variables to Unity
  for(Int_t i = 0; i<LC.TotalLinks ; i++)
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

  //Skip first 40 iteration
  //Touching each internal link with heat bath
  for(int i = 0; i<1; i++)
    {
      for(Int_t ic = 0 ; ic< LC.TotalLinks ; ic++)
	{
	  //cout<<"Processing Link"<<ic<<endl; 
	  LC.UpdateLink(ic); //Apply heat bath algorithm to each internal link and apply periodic boundary condition
      //cout<<"Content of Link "<<LC.LinksContent[ic][0]<<"\t"<<LC.LinksContent[ic][1]<<"\t"<<LC.LinksContent[ic][2]<<"\t"<<LC.LinksContent[ic][3]<<endl;
	}
    }

  
  int TotalPlaquette = 6*LC.Nt*LC.Nx*LC.Ny*LC.Nz;
  int TotalPolyakovLoop = (LC.Nx)*(LC.Ny)*(LC.Nz);
  //int TotalPolyakovLoop = (LC.Nx+1)*(LC.Ny+1)*(LC.Nz+1);
  TComplex *ValueEachPolyakov = new TComplex[TotalPolyakovLoop];  
  double *AveragePolyakov = new double[LC.Run];

  TComplex SumPolyakov;
  double *RUN = new double[LC.Run];
  double *RUNError = new double[LC.Run];

  for(Int_t in=0; in<LC.Run; in++)
    { AveragePolyakov[in] = 0.0;RUN[in]=in;RUNError[in]=in;
    }
  
  for(Int_t Iter=0; Iter<LC.Run; Iter++)
    {
      
      LC.Iteration = Iter;    
      cout<<"Performing Iteration "<<Iter<<endl;  

      //Touching each internal link with heat bath 
      for(Int_t ii=0;ii<1;ii++)
	{
	  for(Int_t ic = 0 ; ic< LC.TotalLinks ; ic++)
	    {
	      //cout<<"Processing Link"<<ic<<endl;
	      LC.UpdateLink(ic); //Apply heat bath algorithm to each internal link and apply periodic boundary condition
	      //cout<<"Content of Link "<<LC.LinksContent[ic][0]<<"\t"<<LC.LinksContent[ic][1]<<"\t"<<LC.LinksContent[ic][2]<<"\t"<<LC.LinksContent[ic][3]<<endl;   
	    }
	}
    
   
  //Evaluating plaquette
   Int_t xt[4] ={0,0,0,0};
   Int_t x1=0,x2=0,x3=0;
   Int_t PolyakovLoopCount=0;
   SumPolyakov = TComplex(0.0,0.0);

   PolyakovLoopCount =0;

	      
	  for(Int_t i=0; i< LC.Nx   ; i++)  
	    {
	      for(Int_t j=0; j< LC.Ny  ; j++)
		{
		  for(Int_t k=0; k< LC.Nz  ; k++) 
		    {
			      PolyakovLoopCount = PolyakovLoopCount + 1;

			      x1 = i; x2=j; x3=k;
			      
			      if(x1 == LC.Nx )
				{
				  x1=0;
				} 
			      if(x2 == LC.Ny )
                                {
                                  x2=0;
                                }
			      if(x3 == LC.Nz )
                                {
                                  x3=0;
                                } 
			      //cout<<"i j k = "<<i<<","<<j<<","<<k<<endl;
			      //cout<<"x y z = "<<x1<<","<<x2<<","<<x3<<endl;

			      ValueEachPolyakov[PolyakovLoopCount-1] = (LC.TracePolyakov(x1,x2,x3));			     
			      SumPolyakov = SumPolyakov + ValueEachPolyakov[PolyakovLoopCount-1];
			      //cout<<"Value of current plaquette is "<<ValueEachPolyakov[PolyakovLoopCount-1]<<endl;
			      
		    }
			  
		}
	    }
	    
	
    
	  //cout<<"Sum = "<<TMath::Sqrt( TMath::Power(TracePL.Re(),2.0) + TMath::Power(TracePL.Im(),2.0) )/PolyakovLoopCount<<endl;
	  //cout<<"Total number of Polyakov loop = "<<PolyakovLoopCount<<" , "<<TotalPolyakovLoop<<endl;
	  //cout<<"SumPolyakov for current configuration is "<<TMath::Sqrt( TMath::Power(SumPolyakov.Re(),2.0) + TMath::Power(SumPolyakov.Im(),2.0) )/PolyakovLoopCount;
	  AveragePolyakov[Iter] = TMath::Sqrt( TMath::Power(SumPolyakov.Re(),2.0) + TMath::Power(SumPolyakov.Im(),2.0) )/PolyakovLoopCount;
    
        
    cout<<Iter<<"\t"<<AveragePolyakov[Iter]<<endl;
    f2<<Iter<<"\t"<<AveragePolyakov[Iter]<<endl;
    }
  
   cout<<"Ending program"<<endl;
   
/* 
   //Plotting the AveragePlaquette vs Iteration
  
  TCanvas *c1 = new TCanvas("C1","c1",500,450);    
  TGraph *GAvg = new TGraph(LC.Run,RUN,AveragePlaquette);
   
   GAvg->Draw("AP");
   GAvg->SetMarkerStyle(6);
   GAvg->SetMarkerSize(1.0);
   GAvg->GetXaxis()->SetRangeUser(-1.,LC.Run);
   GAvg->GetYaxis()->SetRangeUser(0.1,1.0);
   Char_t var1[100];sprintf(var1,"Average Plaquette (Lattice Size= %d ^4, #beta= %.3f, Identity as I.C.) ",LC.Nt,LC.beta);
   GAvg->SetTitle(var1);
   GAvg->GetYaxis()->SetTitle("Average Plaquette");
   GAvg->GetXaxis()->SetTitle("Iteration");
   GAvg->GetYaxis()->CenterTitle();
   GAvg->GetXaxis()->CenterTitle();
   Char_t var[50];
   sprintf(var,"Plaquette_N%d_I%d_Beta%d.png",LC.Nt,LC.Run,int(LC.beta*1000));
   c1->SetGrid();

   //c1->SaveAs("Plaquette_N3_I30_Beta2_3.png");
   //c1->SaveAs(var);
*/
   //Calculating Average Plaquette by averaging for runs after skipping first 10 iteration
   //AveragePlaquette error using standard deviation method
   
   Double_t AP=0.0 , APError=0.0, SDSquare=0.0;
   Int_t SkipN =50;
   for(int i=SkipN; i<LC.Run; i++)
     {
       AP = AP + ( AveragePolyakov[i])/(LC.Run-SkipN);
     }
   for(int i=SkipN; i<LC.Run; i++)
     {
       SDSquare = SDSquare + TMath::Power(AveragePolyakov[i] - AP ,2.0);
     }

   APError = TMath::Sqrt(SDSquare/(LC.Run - SkipN - 1))/TMath::Sqrt(LC.Run - SkipN);
   
   cout<<"Average Over configurations: The average Polyakov =  "<<AP<<" +- "<<APError<<endl;
   f1<<LC.beta<<"\t"<<AP<<"\t"<<APError<<endl;
   f1.close();
   f2.close();

   int EndTime = time(NULL);
   int Hour = (EndTime-StartTime)/3600;
   int Minute = ((EndTime-StartTime)/60)-Hour*60;
   int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
   cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
   return 0;
}
