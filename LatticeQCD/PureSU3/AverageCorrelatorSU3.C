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
  int***** LinksNumber;//To store the link number of link at location [t,x,y,z,d]

  TComplex ID[2][2],Pauli[3][2][2],CPauli[3][2][2],IYOTA,CIYOTA, IDENTITY[3][3]; //Complex Matrix
  TComplex SumProduct3Us[3][3];  //Complex matrix Staple U2*U3*U4 (It is staple)

  void DefineDynamicArrays();
  void LabelLinksCoordinate();
  Int_t OutSideBoundary(Int_t,Int_t ,Int_t, Int_t,Int_t);
  void UpdateLink(Int_t);
  Int_t SearchLinkNumber(Int_t, Int_t, Int_t, Int_t, Int_t); //Replaced by LinksNumber[Nt][Nx][Ny][Nz][4]
  void DefineComplexMatrix();
  int Modulo(int , int );

  //Compute thermal expectation value of all operators
  //void ComputeTraceOfOperator(int t, int x, int y, int z);
  void ComputeFmunuUsingSinglePlaquette();
  void ComputeFmunuUsingFourPlaquette();
  TComplex TraceF3iF3iMinusF4iF4i, TraceF3iF4iPlusF4iF3i, TraceF3iDzF3iMinusF4iDzF4i, TraceF3iDzF4iPlusF4iDzF3i;
  TComplex TraceF3iF3iMinusF4iF4iTraceless, TraceF3iF4iPlusF4iF3iTraceless, TraceF3iDzF3iMinusF4iDzF4iTraceless, TraceF3iDzF4iPlusF4iDzF3iTraceless;

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

  for(int i=0; i<Nt; i++) {for(int j=0; j<Nx; j++) {for(int k=0; k<Ny; k++) {for(int l=0; l<Nz; l++) { LinksNumber[i][j][k][l] = new int [4];}}   }}
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
                  //t-direction                                                           
                  {
                    LinksCoordinate[CurrentLink][0]=a;
                    LinksCoordinate[CurrentLink][1]=b;
                    LinksCoordinate[CurrentLink][2]=c;
                    LinksCoordinate[CurrentLink][3]=d;
                    LinksDirection[CurrentLink]=0;

		    LinksNumber[a][b][c][d][0] = CurrentLink;
		    //std::cout<<LinksCoordinate[CurrentLink][0]<<LinksCoordinate[CurrentLink][1]<<LinksCoordinate[CurrentLink][2]<<LinksCoordinate[CurrentLink][3]<<endl;
                  }
		  
		  //x-dir
		  {
		    LinksCoordinate[CurrentLink + NLinksT][0]=a;
                    LinksCoordinate[CurrentLink + NLinksT][1]=b;
                    LinksCoordinate[CurrentLink + NLinksT][2]=c;
                    LinksCoordinate[CurrentLink + NLinksT][3]=d;
                    LinksDirection[CurrentLink + NLinksT]=1;

                    LinksNumber[a][b][c][d][1] = CurrentLink +  NLinksT ;
		    // std::cout<<LinksCoordinate[CurrentLink+NLinks][0]<<LinksCoordinate[CurrentLink][1]<<LinksCoordinate[CurrentLink][2]<<LinksCoordinate[CurrentLink][3]<<endl;
                  }

                  //y-direction                                                          
                  {
		    LinksCoordinate[CurrentLink + NLinksT + NLinksX][0]=a;
                    LinksCoordinate[CurrentLink + NLinksT + NLinksX][1]=b;
                    LinksCoordinate[CurrentLink + NLinksT + NLinksX][2]=c;
                    LinksCoordinate[CurrentLink + NLinksT + NLinksX][3]=d;
                    LinksDirection[CurrentLink + NLinksT + NLinksX]=2;

                    LinksNumber[a][b][c][d][2] = CurrentLink +  NLinksT + NLinksX;
                  }

                  //z-direction                                                           
                  {
		    LinksCoordinate[CurrentLink + NLinksT + NLinksX + NLinksY][0]=a;
                    LinksCoordinate[CurrentLink + NLinksT + NLinksX + NLinksY][1]=b;
                    LinksCoordinate[CurrentLink + NLinksT + NLinksX + NLinksY][2]=c;
                    LinksCoordinate[CurrentLink + NLinksT + NLinksX + NLinksY][3]=d;
                    LinksDirection[CurrentLink + NLinksT + NLinksX + NLinksY]=3;

                    LinksNumber[a][b][c][d][3] = CurrentLink +  NLinksT + NLinksX + NLinksY;
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
	      LinkNumberStaple[count][i2] = SearchLinkNumber(StapleCoordinate[p][ra][0], StapleCoordinate[p][ra][1], StapleCoordinate[p][ra][2], StapleCoordinate[p][ra][3], StapleCoordinate[p][ra][4]);
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
  //
  for(int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
	{UNew[i][j] = TComplex(0.0,0.0);  UOld[i][j] = TComplex(0.0,0.0);
	  RNew[i][j] = TComplex(0.0,0.0); ROld[i][j] = TComplex(0.0,0.0); AA[i][j] = TComplex(0.0,0.0);
	}}
  
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

  TComplex UCheck[3][3]={0.0}; for(int i = 0; i<3; i++) { for(int j = 0; j<3; j++) {UCheck[i][j]=TComplex(0.0,0.0);}}
 
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

int LatticeCalculation::Modulo(int a, int b)
{
  int m = a % b;
  if (m < 0) {
    // m += (b < 0) ? -b : b; // avoid this form: it is UB when b == INT_MIN                                                                  
    m = (b < 0) ? m - b : m + b;
  }
  return m;
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
      Int_t Temp=t; t = Modulo(Temp, Nt);
    }
  if( x<0 || x > Nx-1 )
    {
      //cout<<"Alert Link is outside the grid, using periodic transformation to bring inside"<<endl;
      Int_t Temp=x; x =  Modulo(Temp, Nx);
    }
  if( y<0 || y > Ny-1 )
    {
      //cout<<"Alert Link is outside the grid, using periodic transformation to bring inside"<<endl;
      Int_t Temp=y; y = Modulo(Temp, Ny);
    }
  if( z<0 || z > Nz-1 )
    {
      //cout<<"Alert Link is outside the grid, using periodic transformation to bring inside"<<endl;
      Int_t Temp=z; z = Modulo(Temp, Nz); //cout<<"z is "<<z<<endl;
    }
  /*
  while(ra< TotalLinks && found == 0)
    {
      if(d == LinksDirection[ra])	{
	  if(t == LinksCoordinate[ra][0] && x == LinksCoordinate[ra][1] && y == LinksCoordinate[ra][2] && z == LinksCoordinate[ra][3])
	    {
	      found = 1;	      Link = ra;
	      
	    }	  else {ra = ra+1;}
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


void LatticeCalculation::ComputeFmunuUsingSinglePlaquette()
{
  //cout<<"\n Begining of a the function"<<endl;
  // [0][3][3] is Fmunu matrix at (t,x,y,z)
  // [1][3][3] is Fmunu matrix at (t,x,y,z+1) 
  //Variable without using Traceless Fmunu
  TraceF3iF3iMinusF4iF4i = TComplex(0.0,0.0);
  TraceF3iF4iPlusF4iF3i  = TComplex(0.0,0.0);
  TraceF3iDzF3iMinusF4iDzF4i = TComplex(0.0,0.0);
  TraceF3iDzF4iPlusF4iDzF3i  = TComplex(0.0,0.0);
  //Variable Using Traceless Fmunu
  TraceF3iF3iMinusF4iF4iTraceless = TComplex(0.0,0.0);
  TraceF3iF4iPlusF4iF3iTraceless  = TComplex(0.0,0.0);
  TraceF3iDzF3iMinusF4iDzF4iTraceless = TComplex(0.0,0.0);
  TraceF3iDzF4iPlusF4iDzF3iTraceless  = TComplex(0.0,0.0);

  TComplex F31[2][3][3], F32[2][3][3], F41[2][3][3], F42[2][3][3]; 
  TComplex F31Traceless[2][3][3], F32Traceless[2][3][3], F41Traceless[2][3][3], F42Traceless[2][3][3];
  TComplex F31Trace[2], F32Trace[2], F41Trace[2], F42Trace[2];
  // DzF31 is Covariant derivative of Fmunu at (t,x,y,z)
  TComplex DzF31[3][3], DzF32[3][3], DzF41[3][3], DzF42[3][3];
  TComplex DzF31Traceless[3][3], DzF32Traceless[3][3], DzF41Traceless[3][3], DzF42Traceless[3][3];
  //LinkMatrix at (t,x,y,z) along z-direction
  TComplex LinkZMatrix[3][3];
  //Plaquette Matrix
  TComplex PlaquetteMatrix[3][3];
  int t=0, x=0, y=0, z=0;int zPlusStep=0; int NPoints=0;

  for(Int_t xt=0; xt<LatticeSize[0]; xt++)  //axis perpendicular to V-axis and H-axis                                                       
    {
      for(Int_t xx=0; xx<LatticeSize[1]; xx++)//axis perpendicular to V-axis and H-axis                                                     
	{
	  for(Int_t xy=0; xy<LatticeSize[2]; xy++)  //Vertical-axis                                                                       
	    {
	      for(Int_t xz=0; xz<LatticeSize[3]; xz++) //Horizontal axis                                                                    
		{
		  t=xt; x=xx; y=xy; z=xz;
		  if(t==3 && x==4 && y==4 && z==4)
		    { NPoints=NPoints+1;
		  //Matrix to compute stuff
		  Int_t Link1=0, Link2=0, Link3=0, Link4=0;
		  for(int i=0; i<3; i++)
		    {  for(int j=0; j<3; j++)
			{ F31[0][i][j]= TComplex(0.0,0.0); F31[1][i][j]= TComplex(0.0,0.0);
			  F32[0][i][j]= TComplex(0.0,0.0); F32[1][i][j]= TComplex(0.0,0.0);
			  F41[0][i][j]= TComplex(0.0,0.0); F41[1][i][j]= TComplex(0.0,0.0);
			  F42[0][i][j]= TComplex(0.0,0.0); F42[1][i][j]= TComplex(0.0,0.0);
			  F31Traceless[0][i][j]= TComplex(0.0,0.0); F31Traceless[1][i][j]= TComplex(0.0,0.0);
			  F32Traceless[0][i][j]= TComplex(0.0,0.0); F32Traceless[1][i][j]= TComplex(0.0,0.0);
                          F41Traceless[0][i][j]= TComplex(0.0,0.0); F41Traceless[1][i][j]= TComplex(0.0,0.0);
			  F42Traceless[0][i][j]= TComplex(0.0,0.0); F42Traceless[1][i][j]= TComplex(0.0,0.0);

			  F31Trace[0]= TComplex(0.0,0.0);   F31Trace[1]= TComplex(0.0,0.0);
			  F32Trace[0]= TComplex(0.0,0.0);   F32Trace[1]= TComplex(0.0,0.0);
                          F41Trace[0]= TComplex(0.0,0.0);   F41Trace[1]= TComplex(0.0,0.0);
			  F42Trace[0]= TComplex(0.0,0.0);   F42Trace[1]= TComplex(0.0,0.0);
			  
			  DzF31[i][j]=TComplex(0.0,0.0); DzF32[i][j]=TComplex(0.0,0.0); 
			  DzF41[i][j]=TComplex(0.0,0.0); DzF42[i][j]=TComplex(0.0,0.0);
			  DzF31Traceless[i][j]=TComplex(0.0,0.0); DzF32Traceless[i][j]=TComplex(0.0,0.0);
                          DzF41Traceless[i][j]=TComplex(0.0,0.0); DzF42Traceless[i][j]=TComplex(0.0,0.0);

			  LinkZMatrix[i][j]=TComplex(0.0,0.0);  PlaquetteMatrix[i][j]=  TComplex(0.0,0.0);
			}}
		  
		  //PlaneType=0 : (3,1) at (t,x,y,z) and at (t,x,y,z+1)
		  //PlaneType=1 : (3,2) at (t,x,y,z) and at (t,x,y,z+1)
		  //PlaneType=2 : (4,1) at (t,x,y,z) and at (t,x,y,z+1)
		  //PlaneType=3 : (4,2) at (t,x,y,z) and at (t,x,y,z+1)
		  for(Int_t PlaneType = 0; PlaneType < 4; PlaneType++)
		    {
		      for(Int_t Step =0; Step <2 ; Step++)
			{ zPlusStep=z+Step;
			  //PlaneType=0 (3,1)
			  if(PlaneType==0)
			    {
			      //cout<<"Plane (3,1) at step="<<Step<<endl;
			      Link1 = SearchLinkNumber(t,  x,    y,  zPlusStep,   3);
			      Link2 = SearchLinkNumber(t,  x,    y,  zPlusStep+1, 1);
			      Link3 = SearchLinkNumber(t,  x+1,  y,  zPlusStep,   3);
			      Link4 = SearchLinkNumber(t,  x,    y,  zPlusStep,   1);
			    }
	  
			  //PlaneType=1 (3,2)
			  if(PlaneType==1)
			    {
			      //cout<<"Plane (3,2) at step="<<Step<<endl;
			      Link1 = SearchLinkNumber(t, x, y,   zPlusStep,    3);
			      Link2 = SearchLinkNumber(t, x, y,   zPlusStep+1,  2);
			      Link3 = SearchLinkNumber(t, x, y+1, zPlusStep,    3);
			      Link4 = SearchLinkNumber(t, x, y,   zPlusStep,    2);
			    }
	  
			  //PlaneType=2 (4,1)
			  if(PlaneType==2)
			    {
			      //cout<<"Plane (4,1) at step="<<Step<<endl;
			      Link1 = SearchLinkNumber(t,   x,   y, zPlusStep, 0);
			      Link2 = SearchLinkNumber(t+1, x,   y, zPlusStep, 1);
			      Link3 = SearchLinkNumber(t,   x+1, y, zPlusStep, 0);
			      Link4 = SearchLinkNumber(t,   x,   y, zPlusStep, 1);
			    }
	    
			  //PlaneType=3 (4,2)
			  if(PlaneType==3)
			    {
			      //cout<<"Plane (4,2) at step="<<Step<<endl;
			      Link1 = SearchLinkNumber(t,   x, y,   zPlusStep, 0);
			      Link2 = SearchLinkNumber(t+1, x, y,   zPlusStep, 2);
			      Link3 = SearchLinkNumber(t,   x, y+1, zPlusStep, 0);
			      Link4 = SearchLinkNumber(t,   x, y,   zPlusStep, 2);
			    }
	    
			  //Reset Plaquette Matrix to zero
			  for(Int_t i=0; i<3; i++){ for(Int_t j=0; j<3; j++) {PlaquetteMatrix[i][j]=TComplex(0.0,0.0); }}
	    
			  //Computing Plaquette Matrix
			  for(Int_t i=0; i<3;i++)
			    {
			      for(Int_t j=0; j<3;j++)
				{
				  for(Int_t m=0; m<3;m++)
				    {
				      for(Int_t n=0; n<3;n++)
					{
					  for(Int_t k=0; k<3;k++)
					    {
					      PlaquetteMatrix[i][j] = PlaquetteMatrix[i][j] + (LinksContent[Link1][i][k]*LinksContent[Link2][k][m]*TComplex::Conjugate(LinksContent[Link3][n][m])*TComplex::Conjugate(LinksContent[Link4][j][n]));  
					    }
					}
				    }
				  //cout<<PlaquetteMatrix[i][j]<<"\t";
				}
			      //cout<<endl;
			    }
			  //cout<<endl;
			  //Computing Fmunu   
			  for(Int_t i=0; i<3; i++)
			    { 
			      for(Int_t j=0; j<3; j++) 
				{
				  if(PlaneType==0)
				    {
				      F31[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA); 
				    }
				  
				  if(PlaneType==1)
				    {
				      F32[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA);
				    }
				  
				  if(PlaneType==2)
				    {
				      F41[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA);
				    }
				  
				  if(PlaneType==3)
				    {
				      F42[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA);
				    }
				}
			    }
			  //cout<<"Step="<<Step<<"\t is done"<<endl;
			} //end of for loop "Step=0,1"      
		    } //end of for loop "PlaneType=0,1,2,3"
		
		  //Link Matrix along z-direction at (t,x,y,z)
		  for(Int_t i=0; i<3; i++)
		    {
		      for(Int_t j=0; j<3; j++)
			{
			  LinkZMatrix[i][j] = LinksContent[ SearchLinkNumber(t, x, y,  z,  3)] [i] [j] ;
			}
		    }
		  
		  //Also compute trace of Fmunu's
		  for(Int_t m=0; m<2; m++)
                    {
		      for(Int_t i=0; i<3; i++)
			{
			  F31Trace[m]= F31Trace[m] + F31[m][i][i];
			  F32Trace[m]= F32Trace[m] + F32[m][i][i];
			  F41Trace[m]= F41Trace[m] + F41[m][i][i];
			  F42Trace[m]= F42Trace[m] + F42[m][i][i];			  
			}
		    }
		  //Now compute Fmunu's as traceless variable
		  for(Int_t m=0; m<2; m++)
                    {
		      for(Int_t i=0; i<3; i++)
			{
			  for(Int_t j=0; j<3; j++)
			    {
			      F31Traceless[m][i][j] = F31[m][i][j] - (F31Trace[m]*IDENTITY[i][j]/3.0);
			      F32Traceless[m][i][j] = F32[m][i][j] - (F32Trace[m]*IDENTITY[i][j]/3.0);
			      F41Traceless[m][i][j] = F41[m][i][j] - (F41Trace[m]*IDENTITY[i][j]/3.0);
			      F42Traceless[m][i][j] = F42[m][i][j] - (F42Trace[m]*IDENTITY[i][j]/3.0);
			    }}}

		  //Covariant-derivative of Fmunu w.r.t z-direction: DzF21[3][3] at (t,x,y,z)
		  for(Int_t i=0; i<3; i++)
		    {
		      for(Int_t j=0; j<3; j++)
			{
			  for(Int_t m=0; m<3; m++)
			    {	      
			      DzF31[i][j] = (DzF31[i][j] - (LinkZMatrix[i][m]*F31[0][m][j]));
			      DzF32[i][j] = (DzF32[i][j] - (LinkZMatrix[i][m]*F32[0][m][j]));
			      DzF41[i][j] = (DzF41[i][j] - (LinkZMatrix[i][m]*F41[0][m][j]));
			      DzF42[i][j] = (DzF42[i][j] - (LinkZMatrix[i][m]*F42[0][m][j]));

			      DzF31Traceless[i][j] = (DzF31Traceless[i][j] - (LinkZMatrix[i][m]*F31Traceless[0][m][j]));
                              DzF32Traceless[i][j] = (DzF32Traceless[i][j] - (LinkZMatrix[i][m]*F32Traceless[0][m][j]));
                              DzF41Traceless[i][j] = (DzF41Traceless[i][j] - (LinkZMatrix[i][m]*F41Traceless[0][m][j]));
                              DzF42Traceless[i][j] = (DzF42Traceless[i][j] - (LinkZMatrix[i][m]*F42Traceless[0][m][j]));
			    }
			  DzF31[i][j] = DzF31[i][j] + F31[1][i][j];
			  DzF32[i][j] = DzF32[i][j] + F32[1][i][j];
			  DzF41[i][j] = DzF41[i][j] + F41[1][i][j];
			  DzF42[i][j] = DzF42[i][j] + F42[1][i][j];
			  
			  DzF31Traceless[i][j] = DzF31Traceless[i][j] + F31Traceless[1][i][j];
			  DzF32Traceless[i][j] = DzF32Traceless[i][j] + F32Traceless[1][i][j];
			  DzF41Traceless[i][j] = DzF41Traceless[i][j] + F41Traceless[1][i][j];
			  DzF42Traceless[i][j] = DzF42Traceless[i][j] + F42Traceless[1][i][j];
			}
		    }
		  
		  //Saving Trace of FFCorrelator and other operators
		  //TraceF3iF3iMinusF4iF4i, TraceF3iF4iPlusF4iF3i, TraceF3iDzF3iMinusF4iDzF4i, TraceF3iDzF4iPlusF4iDzF3i 
		  
		  for(Int_t i=0; i<3; i++)
		    {
		      for(Int_t m=0; m<3; m++)
			{
			  TraceF3iF3iMinusF4iF4i = TraceF3iF3iMinusF4iF4i + (F31[0][i][m]*F31[0][m][i] +  F32[0][i][m]*F32[0][m][i] -  F41[0][i][m]*F41[0][m][i] -  F42[0][i][m]*F42[0][m][i] );
			  TraceF3iF4iPlusF4iF3i = TraceF3iF4iPlusF4iF3i + (F31[0][i][m]*F41[0][m][i] +  F32[0][i][m]*F42[0][m][i]);
			  
			  TraceF3iDzF3iMinusF4iDzF4i  =  TraceF3iDzF3iMinusF4iDzF4i + ( F31[0][i][m]*DzF31[m][i] +  F32[0][i][m]*DzF32[m][i] -  F41[0][i][m]*DzF41[m][i] -  F42[0][i][m]*DzF42[m][i] );
			  TraceF3iDzF4iPlusF4iDzF3i   =  TraceF3iDzF4iPlusF4iDzF3i + ( F31[0][i][m]*DzF41[m][i] + F32[0][i][m]*DzF42[m][i] + F41[0][i][m]*DzF31[m][i] + F42[0][i][m]*DzF32[m][i] );


			  cout<<"SP(i,m)=("<<i<<","<<m<<") F31, F41="<<F31Traceless[0][i][m]/(1.0*NPoints)<<","<<F41Traceless[0][m][i]/(1.0*NPoints)<<", F32, F42="<<F32Traceless[0][i][m]/(1.0*NPoints)<<","<<F42Traceless[0][m][i]/(1.0*NPoints)<<endl;
		          cout<<"SP(i,m)=("<<i<<","<<m<<") F31F41="<<(F31Traceless[0][i][m]*F41Traceless[0][m][i])/(1.0*NPoints)<<" F32F42="<<(F32Traceless[0][i][m]*F42Traceless[0][m][i])/(1.0*NPoints)<<endl;

			  TraceF3iF3iMinusF4iF4iTraceless = TraceF3iF3iMinusF4iF4iTraceless + (F31Traceless[0][i][m]*F31Traceless[0][m][i] +  F32Traceless[0][i][m]*F32Traceless[0][m][i] -  F41Traceless[0][i][m]*F41Traceless[0][m][i] -  F42Traceless[0][i][m]*F42Traceless[0][m][i] );
                          TraceF3iF4iPlusF4iF3iTraceless = TraceF3iF4iPlusF4iF3iTraceless +(F31Traceless[0][i][m]*F41Traceless[0][m][i] +  F32Traceless[0][i][m]*F42Traceless[0][m][i]);

                          TraceF3iDzF3iMinusF4iDzF4iTraceless  =  TraceF3iDzF3iMinusF4iDzF4iTraceless + ( F31Traceless[0][i][m]*DzF31Traceless[m][i] +  F32Traceless[0][i][m]*DzF32Traceless[m][i] - F41Traceless[0][i][m]*DzF41Traceless[m][i] -  F42Traceless[0][i][m]*DzF42Traceless[m][i] );
                          TraceF3iDzF4iPlusF4iDzF3iTraceless   =  TraceF3iDzF4iPlusF4iDzF3iTraceless + ( F31Traceless[0][i][m]*DzF41Traceless[m][i] + F32Traceless[0][i][m]*DzF42Traceless[m][i] + F41Traceless[0][i][m]*DzF31Traceless[m][i] + F42Traceless[0][i][m]*DzF32Traceless[m][i] );
			}
		    }
		  cout<<"SP TraceF3iF4iPlusF4iF3iTraceless"<<TraceF3iF4iPlusF4iF3iTraceless/(1.0*NPoints)<<"\n\n";
		} //if condition test	  
		}}}} //Loops for nt, nx, ny, nz 

  cout<<"LatticePoints="<<LatticePoints<<", and NPoints="<<NPoints<<endl;
  TraceF3iF3iMinusF4iF4i    = TraceF3iF3iMinusF4iF4i/(NPoints*1.0);
  TraceF3iF4iPlusF4iF3i     = 2.0*TraceF3iF4iPlusF4iF3i/(NPoints*1.0);
  TraceF3iDzF3iMinusF4iDzF4i= TraceF3iDzF3iMinusF4iDzF4i/(NPoints*1.0);
  TraceF3iDzF4iPlusF4iDzF3i = TraceF3iDzF4iPlusF4iDzF3i/(NPoints*1.0);

  TraceF3iF3iMinusF4iF4iTraceless     = TraceF3iF3iMinusF4iF4iTraceless/(NPoints*1.0);
  TraceF3iF4iPlusF4iF3iTraceless      = 2.0*TraceF3iF4iPlusF4iF3iTraceless/(NPoints*1.0);
  TraceF3iDzF3iMinusF4iDzF4iTraceless = TraceF3iDzF3iMinusF4iDzF4iTraceless/(NPoints*1.0);
  TraceF3iDzF4iPlusF4iDzF3iTraceless  = TraceF3iDzF4iPlusF4iDzF3iTraceless/(NPoints*1.0);

  
  //cout<<"End of the function "<<endl;
  //return (TraceFinalFmunuDzFmunu.Re());
  
}


void LatticeCalculation::ComputeFmunuUsingFourPlaquette()
{
  //cout<<"\n Begining of a the function"<<endl;
  // [0][3][3] is Fmunu matrix at (t,x,y,z)
  // [1][3][3] is Fmunu matrix at (t,x,y,z+1) 
  //Variable without using Traceless Fmunu
  TraceF3iF3iMinusF4iF4i = TComplex(0.0,0.0);
  TraceF3iF4iPlusF4iF3i  = TComplex(0.0,0.0);
  TraceF3iDzF3iMinusF4iDzF4i = TComplex(0.0,0.0);
  TraceF3iDzF4iPlusF4iDzF3i  = TComplex(0.0,0.0);
  //Variable Using Traceless Fmunu
  TraceF3iF3iMinusF4iF4iTraceless = TComplex(0.0,0.0);
  TraceF3iF4iPlusF4iF3iTraceless  = TComplex(0.0,0.0);
  TraceF3iDzF3iMinusF4iDzF4iTraceless = TComplex(0.0,0.0);
  TraceF3iDzF4iPlusF4iDzF3iTraceless  = TComplex(0.0,0.0);
  
  TComplex F31[2][3][3], F32[2][3][3], F41[2][3][3], F42[2][3][3]; 
  TComplex F31Traceless[2][3][3], F32Traceless[2][3][3], F41Traceless[2][3][3], F42Traceless[2][3][3];
  TComplex F31Trace[2], F32Trace[2], F41Trace[2], F42Trace[2];
  // DzF31 is Covariant derivative of Fmunu at (t,x,y,z)
  TComplex DzF31[3][3], DzF32[3][3], DzF41[3][3], DzF42[3][3];
  TComplex DzF31Traceless[3][3], DzF32Traceless[3][3], DzF41Traceless[3][3], DzF42Traceless[3][3];
  //LinkMatrix at (t,x,y,z) along z-direction
  TComplex LinkZMatrix[3][3];
  //Plaquette Matrix
  TComplex PlaquetteMatrix[3][3];
  int Link1[4], Link2[4], Link3[4], Link4[4];
  int t=0, x=0, y=0, z=0; int zPlusStep=0; int NPoints=0;

  for(Int_t xt=0; xt<LatticeSize[0]; xt++)  //axis perpendicular to V-axis and H-axis                                                       
    {
      for(Int_t xx=0; xx<LatticeSize[1]; xx++)//axis perpendicular to V-axis and H-axis                                                     
	{
	  for(Int_t xy=0; xy<LatticeSize[2]; xy++)  //Vertical-axis                                                                       
	    {
	      for(Int_t xz=0; xz<LatticeSize[3]; xz++) //Horizontal axis                                                                    
		{
		  t=xt; x=xx; y=xy; z=xz;
		  if(t==3 && x==4 && y==4 && z==4)
		    { NPoints=NPoints+1;
		  //Initialize matrix to zero before computing stuff
		  for(int i=0; i<3; i++)
		    {  for(int j=0; j<3; j++)
			{ F31[0][i][j]= TComplex(0.0,0.0); F31[1][i][j]= TComplex(0.0,0.0);
			  F32[0][i][j]= TComplex(0.0,0.0); F32[1][i][j]= TComplex(0.0,0.0);
			  F41[0][i][j]= TComplex(0.0,0.0); F41[1][i][j]= TComplex(0.0,0.0);
			  F42[0][i][j]= TComplex(0.0,0.0); F42[1][i][j]= TComplex(0.0,0.0);
			  F31Traceless[0][i][j]= TComplex(0.0,0.0); F31Traceless[1][i][j]= TComplex(0.0,0.0);
			  F32Traceless[0][i][j]= TComplex(0.0,0.0); F32Traceless[1][i][j]= TComplex(0.0,0.0);
			  F41Traceless[0][i][j]= TComplex(0.0,0.0); F41Traceless[1][i][j]= TComplex(0.0,0.0);
			  F42Traceless[0][i][j]= TComplex(0.0,0.0); F42Traceless[1][i][j]= TComplex(0.0,0.0);

			  F31Trace[0]= TComplex(0.0,0.0);   F31Trace[1]= TComplex(0.0,0.0);
			  F32Trace[0]= TComplex(0.0,0.0);   F32Trace[1]= TComplex(0.0,0.0);
			  F41Trace[0]= TComplex(0.0,0.0);   F41Trace[1]= TComplex(0.0,0.0);
			  F42Trace[0]= TComplex(0.0,0.0);   F42Trace[1]= TComplex(0.0,0.0);
			  
			  DzF31[i][j]=TComplex(0.0,0.0); DzF32[i][j]=TComplex(0.0,0.0);
			  DzF41[i][j]=TComplex(0.0,0.0); DzF42[i][j]=TComplex(0.0,0.0);
			  DzF31Traceless[i][j]=TComplex(0.0,0.0); DzF32Traceless[i][j]=TComplex(0.0,0.0);
			  DzF41Traceless[i][j]=TComplex(0.0,0.0); DzF42Traceless[i][j]=TComplex(0.0,0.0);
			  LinkZMatrix[i][j]=TComplex(0.0,0.0);  PlaquetteMatrix[i][j]=  TComplex(0.0,0.0);
			  Link1[i]=0, Link2[i]=0, Link3[i]=0, Link4[i]=0;
			  Link1[i+1]=0, Link2[i+1]=0, Link3[i+1]=0, Link4[i+1]=0;
			}}
		
		  //PlaneType=0 : (3,1) at (t,x,y,z) and at (t,x,y,z+1)
		  //PlaneType=1 : (3,2) at (t,x,y,z) and at (t,x,y,z+1)
		  //PlaneType=2 : (4,1) at (t,x,y,z) and at (t,x,y,z+1)
		  //PlaneType=3 : (4,2) at (t,x,y,z) and at (t,x,y,z+1)
		  for(Int_t PlaneType = 0; PlaneType < 4; PlaneType++)
		    {
		      for(Int_t Step =0; Step <2 ; Step++)
			{ zPlusStep = z+Step;
			  //PlaneType=0 (3,1)
			  if(PlaneType==0)
			    {
			      //First quadrant: (+3,+1) dir
			      Link1[0] = SearchLinkNumber(t,  x,    y,  zPlusStep,   3);
			      Link2[0] = SearchLinkNumber(t,  x,    y,  zPlusStep+1, 1);
			      Link3[0] = SearchLinkNumber(t,  x+1,  y,  zPlusStep,   3);
			      Link4[0] = SearchLinkNumber(t,  x,    y,  zPlusStep,   1);
			      //Second quadrant: (-3,+1) dir				
			      Link1[1] = SearchLinkNumber(t,  x,    y,  zPlusStep  -1,   3);
			      Link2[1] = SearchLinkNumber(t,  x,    y,  zPlusStep+1-1,   1);
			      Link3[1] = SearchLinkNumber(t,  x+1,  y,  zPlusStep  -1,   3);
			      Link4[1] = SearchLinkNumber(t,  x,    y,  zPlusStep  -1,   1);
                              //Third quadrant:  (-3,-1) dir
			      Link1[2] = SearchLinkNumber(t,  x-1,    y,  zPlusStep-1,   3);
			      Link2[2] = SearchLinkNumber(t,  x-1,    y,  zPlusStep,   1);
			      Link3[2] = SearchLinkNumber(t,  x,      y,  zPlusStep-1,   3);
			      Link4[2] = SearchLinkNumber(t,  x-1,    y,  zPlusStep-1,   1);
			      //Fourth quadrant: (+3,-1) dir
			      Link1[3] = SearchLinkNumber(t,  x-1,    y,  zPlusStep,   3);
			      Link2[3] = SearchLinkNumber(t,  x-1,    y,  zPlusStep+1, 1);
			      Link3[3] = SearchLinkNumber(t,  x,      y,  zPlusStep,   3);
			      Link4[3] = SearchLinkNumber(t,  x-1,    y,  zPlusStep,   1);
			    }
	                  
			  //PlaneType=1 (3,2)
			  if(PlaneType==1)
			    {
			      //cout<<"Plane (3,2) at step="<<Step<<endl;
                              //First quadrant: (+3,+2) dir
			      Link1[0] = SearchLinkNumber(t, x, y,   zPlusStep,    3);
			      Link2[0] = SearchLinkNumber(t, x, y,   zPlusStep+1,  2);
			      Link3[0] = SearchLinkNumber(t, x, y+1, zPlusStep,    3);
			      Link4[0] = SearchLinkNumber(t, x, y,   zPlusStep,    2);
			      //Second quadrant:(-3,+2) dir
			      Link1[1] = SearchLinkNumber(t, x, y,   zPlusStep-1,    3);
			      Link2[1] = SearchLinkNumber(t, x, y,   zPlusStep,  2);
			      Link3[1] = SearchLinkNumber(t, x, y+1, zPlusStep-1,    3);
			      Link4[1] = SearchLinkNumber(t, x, y,   zPlusStep-1,    2);
			      //Third quadrant:(-3,-2) dir
			      Link1[2] = SearchLinkNumber(t, x, y-1,   zPlusStep-1,    3);
			      Link2[2] = SearchLinkNumber(t, x, y-1,   zPlusStep,  2);
			      Link3[2] = SearchLinkNumber(t, x, y,     zPlusStep-1,    3);
			      Link4[2] = SearchLinkNumber(t, x, y-1,   zPlusStep-1,    2);
                              //Fourth quadrant: (+3,-2) dir
			      Link1[3] = SearchLinkNumber(t, x, y-1,   zPlusStep,    3);
			      Link2[3] = SearchLinkNumber(t, x, y-1,   zPlusStep+1,  2);
			      Link3[3] = SearchLinkNumber(t, x, y,     zPlusStep,    3);
			      Link4[3] = SearchLinkNumber(t, x, y-1,   zPlusStep,    2);
			    }
	  
			  //PlaneType=2 (4,1)
			  if(PlaneType==2)
			    {
			      //cout<<"Plane (4,1) at step="<<Step<<endl;
                              //First quadrant: (+4,+1) dir
			      Link1[0] = SearchLinkNumber(t,   x,   y, zPlusStep, 0);
			      Link2[0] = SearchLinkNumber(t+1, x,   y, zPlusStep, 1);
			      Link3[0] = SearchLinkNumber(t,   x+1, y, zPlusStep, 0);
			      Link4[0] = SearchLinkNumber(t,   x,   y, zPlusStep, 1);
                              //Second quadrant: (-4,+1) dir
			      Link1[1] = SearchLinkNumber(t-1,   x,   y, zPlusStep, 0);
			      Link2[1] = SearchLinkNumber(t,     x,   y, zPlusStep, 1);
			      Link3[1] = SearchLinkNumber(t-1,   x+1, y, zPlusStep, 0);
			      Link4[1] = SearchLinkNumber(t-1,   x,   y, zPlusStep, 1);
                              //Third quadrant: (-4,-1) dir
			      Link1[2] = SearchLinkNumber(t-1,   x-1,   y, zPlusStep, 0);
			      Link2[2] = SearchLinkNumber(t,     x-1,   y, zPlusStep, 1);
			      Link3[2] = SearchLinkNumber(t-1,   x,     y, zPlusStep, 0);
			      Link4[2] = SearchLinkNumber(t-1,   x-1,   y, zPlusStep, 1);
                              //Fourth quadrant: (+4,-1) dir
			      Link1[3] = SearchLinkNumber(t,   x-1,   y, zPlusStep, 0);
			      Link2[3] = SearchLinkNumber(t+1, x-1,   y, zPlusStep, 1);
			      Link3[3] = SearchLinkNumber(t,   x,     y, zPlusStep, 0);
			      Link4[3] = SearchLinkNumber(t,   x-1,   y, zPlusStep, 1);
			    }
	    
			  //PlaneType=3 (4,2)
			  if(PlaneType==3)
			    {
			      //cout<<"Plane (4,2) at step="<<Step<<endl;
                              //First quadrant: (+4,+2) dir
			      Link1[0] = SearchLinkNumber(t,   x, y,   zPlusStep, 0);
			      Link2[0] = SearchLinkNumber(t+1, x, y,   zPlusStep, 2);
			      Link3[0] = SearchLinkNumber(t,   x, y+1, zPlusStep, 0);
			      Link4[0] = SearchLinkNumber(t,   x, y,   zPlusStep, 2);
                              //Second quadrant: (-4,+2) dir
			      Link1[1] = SearchLinkNumber(t-1,   x, y,   zPlusStep, 0);
			      Link2[1] = SearchLinkNumber(t,     x, y,   zPlusStep, 2);
			      Link3[1] = SearchLinkNumber(t-1,   x, y+1, zPlusStep, 0);
			      Link4[1] = SearchLinkNumber(t-1,   x, y,   zPlusStep, 2);
                              //Third quadrant: (-4,-2) dir
			      Link1[2] = SearchLinkNumber(t-1,   x, y-1,   zPlusStep, 0);
			      Link2[2] = SearchLinkNumber(t,     x, y-1,   zPlusStep, 2);
			      Link3[2] = SearchLinkNumber(t-1,   x, y,     zPlusStep, 0);
			      Link4[2] = SearchLinkNumber(t-1,   x, y-1,   zPlusStep, 2);
                              //Fourth quadrant: (+4,-2) dir
			      Link1[3] = SearchLinkNumber(t,   x, y-1,   zPlusStep, 0);
			      Link2[3] = SearchLinkNumber(t+1, x, y-1,   zPlusStep, 2);
			      Link3[3] = SearchLinkNumber(t,   x, y,     zPlusStep, 0);
			      Link4[3] = SearchLinkNumber(t,   x, y-1,   zPlusStep, 2);
			    }
	                    
			  //Reset Plaquette Matrix to zero
			  for(Int_t i=0; i<3; i++){ for(Int_t j=0; j<3; j++) {PlaquetteMatrix[i][j]=TComplex(0.0,0.0); }}
	                  
			  //Computing Plaquette Matrix
			  for(Int_t i=0; i<3;i++)
			    {
			      for(Int_t j=0; j<3;j++)
				{
				  for(Int_t m=0; m<3;m++)
				    {
				      for(Int_t n=0; n<3;n++)
					{
					  for(Int_t k=0; k<3;k++)
					    {
						 PlaquetteMatrix[i][j] = PlaquetteMatrix[i][j] + (LinksContent[Link1[0]][i][k]*LinksContent[Link2[0]][k][m]*TComplex::Conjugate(LinksContent[Link3[0]][n][m])*TComplex::Conjugate(LinksContent[Link4[0]][j][n]));
						 PlaquetteMatrix[i][j] = PlaquetteMatrix[i][j] + (LinksContent[Link2[1]][i][k]*TComplex::Conjugate(LinksContent[Link3[1]][m][k])*TComplex::Conjugate(LinksContent[Link4[1]][n][m])*LinksContent[Link1[1]][n][j]);
						 PlaquetteMatrix[i][j] = PlaquetteMatrix[i][j] + (TComplex::Conjugate(LinksContent[Link3[2]][k][i])*TComplex::Conjugate(LinksContent[Link4[2]][m][k])*LinksContent[Link1[2]][m][n]*LinksContent[Link2[2]][n][j]);
						 PlaquetteMatrix[i][j] = PlaquetteMatrix[i][j] + (TComplex::Conjugate(LinksContent[Link4[3]][k][i])*LinksContent[Link1[3]][k][m]*LinksContent[Link2[3]][m][n]*TComplex::Conjugate(LinksContent[Link3[3]][j][n]));
						}
					}
				    }
				  //cout<<PlaquetteMatrix[i][j]<<"\t";
				}
			      //cout<<endl;
			    }//For loop int i=0,1,2, i.e. Plaquette matrix calculation is done (Clover SUM)
			  
			  //cout<<endl;
			  //Computing Fmunu   
			  for(Int_t i=0; i<3; i++)
			    { 
			      for(Int_t j=0; j<3; j++) 
				{
				  if(PlaneType==0)
				    {
				      F31[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA); 
				    }
				  
				  if(PlaneType==1)
				    {
				      F32[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA);
				    }
				  
				  if(PlaneType==2)
				    {
				      F41[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA);
				    }
				  
				  if(PlaneType==3)
				    {
				      F42[Step][i][j] = (PlaquetteMatrix[i][j] -  TComplex::Conjugate( PlaquetteMatrix[j][i] ))/(2.0*IYOTA);
				    }
				}
			    }
			  //cout<<"Step="<<Step<<"\t is done"<<endl;
			} //end of for loop "Step=0,1"      
		    }//end of for loop "PlaneType=0,1,2,3
		  
		  //Link Matrix along z-direction at (t,x,y,z)
		  for(Int_t i=0; i<3; i++)
		    {
		      for(Int_t j=0; j<3; j++)
			{
			  LinkZMatrix[i][j] = LinksContent[ SearchLinkNumber(t, x, y,  z,  3)] [i] [j] ;
			}
		    }
		  
		  //Also compute trace of Fmunu's
		  for(Int_t m=0; m<2; m++)
		    {
		      for(Int_t i=0; i<3; i++)
			{
			  F31Trace[m]= F31Trace[m] + F31[m][i][i];
			  F32Trace[m]= F32Trace[m] + F32[m][i][i];
			  F41Trace[m]= F41Trace[m] + F41[m][i][i];
			  F42Trace[m]= F42Trace[m] + F42[m][i][i];
			}
		    }
		  //Now compute Fmunu's as traceless variable
		  for(Int_t m=0; m<2; m++)
		    {
		      for(Int_t i=0; i<3; i++)
			{
			  for(Int_t j=0; j<3; j++)
			    {
			      F31Traceless[m][i][j] = F31[m][i][j] - (F31Trace[m]*IDENTITY[i][j]/3.0);
			      F32Traceless[m][i][j] = F32[m][i][j] - (F32Trace[m]*IDENTITY[i][j]/3.0);
			      F41Traceless[m][i][j] = F41[m][i][j] - (F41Trace[m]*IDENTITY[i][j]/3.0);
			      F42Traceless[m][i][j] = F42[m][i][j] - (F42Trace[m]*IDENTITY[i][j]/3.0);
			    }}}
		
		  //Covariant-derivative of Fmunu w.r.t z-direction: DzF21[3][3] at (t,x,y,z)
		  for(Int_t i=0; i<3; i++)
		    {
		      for(Int_t j=0; j<3; j++)
			{
			  for(Int_t m=0; m<3; m++)
			    {
			      DzF31[i][j] = (DzF31[i][j] - (LinkZMatrix[i][m]*F31[0][m][j]));
			      DzF32[i][j] = (DzF32[i][j] - (LinkZMatrix[i][m]*F32[0][m][j]));
			      DzF41[i][j] = (DzF41[i][j] - (LinkZMatrix[i][m]*F41[0][m][j]));
			      DzF42[i][j] = (DzF42[i][j] - (LinkZMatrix[i][m]*F42[0][m][j]));

			      DzF31Traceless[i][j] = (DzF31Traceless[i][j] - (LinkZMatrix[i][m]*F31Traceless[0][m][j]));
			      DzF32Traceless[i][j] = (DzF32Traceless[i][j] - (LinkZMatrix[i][m]*F32Traceless[0][m][j]));
			      DzF41Traceless[i][j] = (DzF41Traceless[i][j] - (LinkZMatrix[i][m]*F41Traceless[0][m][j]));
			      DzF42Traceless[i][j] = (DzF42Traceless[i][j] - (LinkZMatrix[i][m]*F42Traceless[0][m][j]));
			    }
			  DzF31[i][j] = DzF31[i][j] + F31[1][i][j];
			  DzF32[i][j] = DzF32[i][j] + F32[1][i][j];
			  DzF41[i][j] = DzF41[i][j] + F41[1][i][j];
			  DzF42[i][j] = DzF42[i][j] + F42[1][i][j];

			  DzF31Traceless[i][j] = DzF31Traceless[i][j] + F31Traceless[1][i][j];
			  DzF32Traceless[i][j] = DzF32Traceless[i][j] + F32Traceless[1][i][j];
			  DzF41Traceless[i][j] = DzF41Traceless[i][j] + F41Traceless[1][i][j];
			  DzF42Traceless[i][j] = DzF42Traceless[i][j] + F42Traceless[1][i][j];
			}
		    }
		  
		  //Saving Trace of FFCorrelator and other operators
		  //TraceF3iF3iMinusF4iF4i, TraceF3iF4iPlusF4iF3i, TraceF3iDzF3iMinusF4iDzF4i, TraceF3iDzF4iPlusF4iDzF3i 
		  
		  for(Int_t i=0; i<3; i++)
		    {
		      for(Int_t m=0; m<3; m++)
			{
			  TraceF3iF3iMinusF4iF4i = TraceF3iF3iMinusF4iF4i + (F31[0][i][m]*F31[0][m][i] +  F32[0][i][m]*F32[0][m][i] -  F41[0][i][m]*F41[0][m][i] -  F42[0][i][m]*F42[0][m][i] );
			  TraceF3iF4iPlusF4iF3i = TraceF3iF4iPlusF4iF3i + (F31[0][i][m]*F41[0][m][i] +  F32[0][i][m]*F42[0][m][i]);

			  TraceF3iDzF3iMinusF4iDzF4i  =  TraceF3iDzF3iMinusF4iDzF4i + ( F31[0][i][m]*DzF31[m][i] +  F32[0][i][m]*DzF32[m][i] -  F41[0][i][m]*DzF41[m][i] -  F42[0][i][m]*DzF42[m][i] );
			  TraceF3iDzF4iPlusF4iDzF3i   =  TraceF3iDzF4iPlusF4iDzF3i + ( F31[0][i][m]*DzF41[m][i] + F32[0][i][m]*DzF42[m][i] + F41[0][i][m]*DzF31[m][i] + F42[0][i][m]*DzF32[m][i] );

			  TraceF3iF3iMinusF4iF4iTraceless = TraceF3iF3iMinusF4iF4iTraceless + (F31Traceless[0][i][m]*F31Traceless[0][m][i] +  F32Traceless[0][i][m]*F32Traceless[0][m][i] -  F41Traceless[0][i][m]*F41Traceless[0][m][i] -  F42Traceless[0][i][m]*F42Traceless[0][m][i] );
			  TraceF3iF4iPlusF4iF3iTraceless = TraceF3iF4iPlusF4iF3iTraceless +(F31Traceless[0][i][m]*F41Traceless[0][m][i] +  F32Traceless[0][i][m]*F42Traceless[0][m][i]);

			  cout<<"Clover(i,m)=("<<i<<","<<m<<") F31, F41="<<F31Traceless[0][i][m]/(4.0*NPoints)<<","<<F41Traceless[0][m][i]/(4.0*NPoints)<<", F32, F42="<<F32Traceless[0][i][m]/(4.0*NPoints)<<","<<F42Traceless[0][m][i]/(4.0*NPoints)<<endl;
			  cout<<"Clover(i,m)=("<<i<<","<<m<<") F31F41="<<(F31Traceless[0][i][m]*F41Traceless[0][m][i])/(16.0*NPoints)<<" F32F42="<<(F32Traceless[0][i][m]*F42Traceless[0][m][i])/(16.0*NPoints)<<endl;
			  
			  TraceF3iDzF3iMinusF4iDzF4iTraceless  =  TraceF3iDzF3iMinusF4iDzF4iTraceless + ( F31Traceless[0][i][m]*DzF31Traceless[m][i] +  F32Traceless[0][i][m]*DzF32Traceless[m][i] - F41Traceless[0][i][m]*DzF41Traceless[m][i] -  F42Traceless[0][i][m]*DzF42Traceless[m][i] );
			  TraceF3iDzF4iPlusF4iDzF3iTraceless   =  TraceF3iDzF4iPlusF4iDzF3iTraceless + ( F31Traceless[0][i][m]*DzF41Traceless[m][i] + F32Traceless[0][i][m]*DzF42Traceless[m][i] + F41Traceless[0][i][m]*DzF31Traceless[m][i] + F42Traceless[0][i][m]*DzF32Traceless[m][i] );
			  
			}
		    }

		  //cout<<"Clover TraceF3iF4iPlusF4iF3iTraceless"<<TraceF3iF4iPlusF4iF3iTraceless<<"\n";
		  cout<<"Clover TraceF3iF4iPlusF4iF3iTraceless"<<TraceF3iF4iPlusF4iF3iTraceless/(16.0*NPoints)<<"\n\n";} //if condition test
		}}}} //Loop over Nt, Nx, Ny, Nz

  cout<<"LatticePoints="<<LatticePoints<<", and NPoints="<<NPoints<<endl;
  TraceF3iF3iMinusF4iF4i    = TraceF3iF3iMinusF4iF4i/(NPoints*16.0);
  TraceF3iF4iPlusF4iF3i     = 2.0*TraceF3iF4iPlusF4iF3i/(NPoints*16.0);
  TraceF3iDzF3iMinusF4iDzF4i= TraceF3iDzF3iMinusF4iDzF4i/(NPoints*16.0);
  TraceF3iDzF4iPlusF4iDzF3i = TraceF3iDzF4iPlusF4iDzF3i/(NPoints*16.0);

  TraceF3iF3iMinusF4iF4iTraceless     = TraceF3iF3iMinusF4iF4iTraceless/(NPoints*16.0);
  TraceF3iF4iPlusF4iF3iTraceless      = 2.0*TraceF3iF4iPlusF4iF3iTraceless/(NPoints*16.0);
  TraceF3iDzF3iMinusF4iDzF4iTraceless = TraceF3iDzF3iMinusF4iDzF4iTraceless/(NPoints*16.0);
  TraceF3iDzF4iPlusF4iDzF3iTraceless  = TraceF3iDzF4iPlusF4iDzF3iTraceless/(NPoints*16.0);

  //return (TraceFinalFmunuDzFmunu.Re());
  
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
  //3x3 Identity matrix
  IDENTITY[0][0]=TComplex(1.0,0.0); IDENTITY[0][1]=TComplex(0.0,0.0); IDENTITY[0][2]=TComplex(0.0,0.0);
  IDENTITY[1][0]=TComplex(0.0,0.0); IDENTITY[1][1]=TComplex(1.0,0.0); IDENTITY[1][2]=TComplex(0.0,0.0);
  IDENTITY[2][0]=TComplex(0.0,0.0); IDENTITY[2][1]=TComplex(0.0,0.0); IDENTITY[2][2]=TComplex(1.0,0.0);
  CIYOTA = TComplex(0.0,-1.0);

}
int main(int argc, char* argv[])
//int main(int nt,int nx,double Beta)
{ int StartTime = time(NULL);

  //const double PI = 3.1415926; //Value of pi
  //const int NIteration = 10;
  //srand(time(NULL));
  LatticeCalculation LC; 
  
  //..............Input Parameter.....Lattice.....
  LC.PI = 3.1415926;
  LC.Nt = atoi(argv[1]); //Discretizing each dimension in N-points 
  LC.Nx = atoi(argv[2]);
  LC.Ny = atoi(argv[2]);
  LC.Nz = atoi(argv[2]);
  LC.LatticePoints = LC.Nt*LC.Nx*LC.Ny*LC.Nz; //                            
  LC.TotalLinks = 4*LC.Nt*LC.Nx*LC.Ny*LC.Nz; // 4 space-time dimension      
  LC.NLinksX = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in X-direction        
  LC.NLinksY = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Y-direction        
  LC.NLinksZ = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in Z-direction        
  LC.NLinksT = LC.Nt*LC.Nx*LC.Ny*LC.Nz;      // Links in T-direction        
  LC.Run = 100;                //Number of iteration to get stable gauge field configuration                                     
  LC.TotalConfiguration  = atoi(argv[3]); //Number of gauge field configuartion over which average of (F3iF3i-F4iF4i) is computed         
  LC.TotalTemperatureRun = 1; //Number of different values of Temperature to see the dependence of correlation function         
  LC.lambdaLattice = 5.3; //(in MeV) A parameter that goes in LatticeSpacing formular using renormalization group equation 
  
  LC.DefineDynamicArrays();
  LC.DefineComplexMatrix();
  //Labeling each link with its 4 space-time coordinate and a direction
  LC.LabelLinksCoordinate();
  
  //..............End of Input Parameter........Lattice
  cout<<"Program started"<<endl;
  cout<<"Dimension of Lattice used is  \t" <<LC.Nt<<"*"<<LC.Nx<<"*"<<LC.Ny<<"*"<<LC.Nz<<", So Total Links are "<<LC.TotalLinks<<endl;
  cout<<"\t"<< LC.LatticeSize[0]<<","<<LC.LatticeSize[1]<<","<<LC.LatticeSize[2]<<","<<LC.LatticeSize[3]<<endl;
  
  LC.beta = atof(argv[4]);//1.80 + 0.1*TemperatureLoop;
  LC.BareCouplingConstant = 1.0;
  LC.LatticeSpacing =1.0; //
  LC.temperature = 1.0;//
  
  cout<<"Coupling is "<<LC.BareCouplingConstant<<endl;
  cout<<"Beta (6/g^2)= "<< LC.beta<<endl;
  cout<<"LatticeSpacing "<<LC.LatticeSpacing<<endl;
  cout<<"Temperature is "<<LC.temperature<<endl;
  
  
  ofstream f0, f1,f2,f3, f4, f5, f6, f7, f8;
  char FileName1[25000],FileName2[25000],FileName3[25000], FileName4[25000];
  char FileName5[25000],FileName6[25000],FileName7[25000], FileName8[25000], FileName0[25000];
  sprintf(FileName1,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataLO_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);
  sprintf(FileName2,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataNLO_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);
  sprintf(FileName3,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataLO_Traceless_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);
  sprintf(FileName4,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataNLO_Traceless_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);

  sprintf(FileName5,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataLO_Clover_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);
  sprintf(FileName6,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataNLO_Clover_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);
  sprintf(FileName7,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataLO_Clover_Traceless_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);
  sprintf(FileName8,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Output/FFCorrelatorEachConfigurationDataNLO_Clover_Traceless_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);

  sprintf(FileName0,"/wsu/home/fy/fy41/fy4125/Lattice/SU3/Log/TimeFFCLO_NLO_Clover_Traceless_Nt%i_Ns%i_Beta_%.4f.txt",LC.Nt,LC.Nx,LC.beta);
  f1.open(FileName1,ios::out);  f2.open(FileName2,ios::out);f3.open(FileName3,ios::out);  f4.open(FileName4,ios::out);
  f5.open(FileName5,ios::out);  f6.open(FileName6,ios::out);f7.open(FileName7,ios::out);  f8.open(FileName8,ios::out);
  f0.open(FileName0,ios::out);
  
  f1<<"#Iter \t TraceF3iF3iMinusF4iF4i.Re \t TraceF3iF3iMinusF4iF4i.Im \t CumulativeF3iF3iMinusF4iF4i.Re \t CumulativeF3iF3iMinusF4iF4i.Im \t TraceF3iF4iPlusF4iF3i.Re  \t TraceF3iF4iPlusF4iF3i.Im \t CumulativeF3iF4iPlusF4iF3i.Re \t CumulativeF3iF4iPlusF4iF3i.Im"<<endl;
      f2<<"#Iter \t TraceF3iDzF3iMinusF4iDzF4i.Re \t TraceF3iDzF3iMinusF4iDzF4i.Im \t CumulativeF3iDzF3iMinusF4iDzF4i.Re \t CumulativeF3iDzF3iMinusF4iDzF4i.Im \t TraceF3iDzF4iPlusF4iDzF3i.Re \t TraceF3iDzF4iPlusF4iDzF3i.Im \t CumulativeF3iDzF4iPlusF4iDzF3i.Re  \t CumulativeF3iDzF4iPlusF4iDzF3i.Im "<<endl;

      f3<<"#Iter \t TraceF3iF3iMinusF4iF4iTraceless.Re \t TraceF3iF3iMinusF4iF4iTraceless.Im \t CumulativeF3iF3iMinusF4iF4iTraceless.Re \t CumulativeF3iF3iMinusF4iF4iTraceless.Im \t TraceF3iF4iPlusF4iF3iTraceless.Re \t TraceF3iF4iPlusF4iF3iTraceless.Im \t CumulativeF3iF4iPlusF4iF3iTraceless.Re \t CumulativeF3iF4iPlusF4iF3iTraceless.Im"<<endl;
      f4<<"#Iter \t TraceF3iDzF3iMinusF4iDzF4iTraceless.Re \t TraceF3iDzF3iMinusF4iDzF4iTraceless.Im \t CumulativeF3iDzF3iMinusF4iDzF4iTraceless.Re \t CumulativeF3iDzF3iMinusF4iDzF4iTraceless.Im \t TraceF3iDzF4iPlusF4iDzF3iTraceless.Re \t TraceF3iDzF4iPlusF4iDzF3iTraceless.Im \t CumulativeF3iDzF4iPlusF4iDzF3iTraceless.Re \t CumulativeF3iDzF4iPlusF4iDzF3iTraceless.Im"<<endl;
      
      //Colver based calculation of Fmunu
      f5<<"#Iter \t TraceF3iF3iMinusF4iF4iClover.Re \t TraceF3iF3iMinusF4iF4iClover.Im \t CumulativeF3iF3iMinusF4iF4iClover.Re \t CumulativeF3iF3iMinusF4iF4iClover.Im \t TraceF3iF4iPlusF4iF3iClover.Re \t TraceF3iF4iPlusF4iF3iClover.Im \t CumulativeF3iF4iPlusF4iF3iClover.Re \t CumulativeF3iF4iPlusF4iF3iClover.Im"<<endl;
      f6<<"#Iter \t TraceF3iDzF3iMinusF4iDzF4iClover.Re \t TraceF3iDzF3iMinusF4iDzF4iClover.Im \t CumulativeF3iDzF3iMinusF4iDzF4iClover.Re \t CumulativeF3iDzF3iMinusF4iDzF4iClover.Im \t TraceF3iDzF4iPlusF4iDzF3iClover.Re \t TraceF3iDzF4iPlusF4iDzF3iClover.Im \t CumulativeF3iDzF4iPlusF4iDzF3iClover.Re \t CumulativeF3iDzF4iPlusF4iDzF3iClover.Im"<<endl;

      f7<<"#Iter \t TraceF3iF3iMinusF4iF4iCloverTraceless.Re \t TraceF3iF3iMinusF4iF4iCloverTraceless.Im \t CumulativeF3iF3iMinusF4iF4iCloverTraceless.Re \t CumulativeF3iF3iMinusF4iF4iCloverTraceless.Im \t TraceF3iF4iPlusF4iF3iCloverTraceless.Re \t TraceF3iF4iPlusF4iF3iCloverTraceless.Im \t CumulativeF3iF4iPlusF4iF3iCloverTraceless.Re \t CumulativeF3iF4iPlusF4iF3iCloverTraceless.Im"<<endl;
      f8<<"#Iter \t TraceF3iDzF3iMinusF4iDzF4iCloverTraceless.Re \t TraceF3iDzF3iMinusF4iDzF4iCloverTraceless.Im \t CumulativeF3iDzF3iMinusF4iDzF4iCloverTraceless.Re \t CumulativeF3iDzF3iMinusF4iDzF4iCloverTraceless.Im \t TraceF3iDzF4iPlusF4iDzF3iCloverTraceless.Re \t TraceF3iDzF4iPlusF4iDzF3iCloverTraceless.Im \t CumulativeF3iDzF4iPlusF4iDzF3iCloverTraceless.Re \t CumulativeF3iDzF4iPlusF4iDzF3iCloverTraceless.Im()"<<endl;
  //Variable for Fmunu calculations
  TComplex TraceF3iF3iMinusF4iF4i,     TraceF3iF3iMinusF4iF4iTraceless,     TraceF3iF4iPlusF4iF3i,     TraceF3iF4iPlusF4iF3iTraceless;
  TComplex TraceF3iDzF3iMinusF4iDzF4i, TraceF3iDzF3iMinusF4iDzF4iTraceless, TraceF3iDzF4iPlusF4iDzF3i, TraceF3iDzF4iPlusF4iDzF3iTraceless;
  TComplex SumF3iF3iMinusF4iF4i     = TComplex(0.0,0.0), CumulativeF3iF3iMinusF4iF4i     = TComplex(0.0,0.0);
  TComplex SumF3iF4iPlusF4iF3i      = TComplex(0.0,0.0), CumulativeF3iF4iPlusF4iF3i      = TComplex(0.0,0.0);
  TComplex SumF3iDzF3iMinusF4iDzF4i = TComplex(0.0,0.0), CumulativeF3iDzF3iMinusF4iDzF4i = TComplex(0.0,0.0);
  TComplex SumF3iDzF4iPlusF4iDzF3i  = TComplex(0.0,0.0), CumulativeF3iDzF4iPlusF4iDzF3i  = TComplex(0.0,0.0);

  TComplex SumF3iF3iMinusF4iF4iTraceless     = TComplex(0.0,0.0), CumulativeF3iF3iMinusF4iF4iTraceless     = TComplex(0.0,0.0);
  TComplex SumF3iF4iPlusF4iF3iTraceless      = TComplex(0.0,0.0), CumulativeF3iF4iPlusF4iF3iTraceless      = TComplex(0.0,0.0);
  TComplex SumF3iDzF3iMinusF4iDzF4iTraceless = TComplex(0.0,0.0), CumulativeF3iDzF3iMinusF4iDzF4iTraceless = TComplex(0.0,0.0);
  TComplex SumF3iDzF4iPlusF4iDzF3iTraceless  = TComplex(0.0,0.0), CumulativeF3iDzF4iPlusF4iDzF3iTraceless  = TComplex(0.0,0.0);

  TComplex TraceF3iF3iMinusF4iF4iClover,TraceF3iF3iMinusF4iF4iCloverTraceless;
  TComplex TraceF3iF4iPlusF4iF3iClover, TraceF3iF4iPlusF4iF3iCloverTraceless;
  TComplex TraceF3iDzF3iMinusF4iDzF4iClover, TraceF3iDzF3iMinusF4iDzF4iCloverTraceless;
  TComplex TraceF3iDzF4iPlusF4iDzF3iClover, TraceF3iDzF4iPlusF4iDzF3iCloverTraceless;
  TComplex SumF3iF3iMinusF4iF4iClover     = TComplex(0.0,0.0), CumulativeF3iF3iMinusF4iF4iClover     = TComplex(0.0,0.0);
  TComplex SumF3iF4iPlusF4iF3iClover      = TComplex(0.0,0.0), CumulativeF3iF4iPlusF4iF3iClover      = TComplex(0.0,0.0);
  TComplex SumF3iDzF3iMinusF4iDzF4iClover = TComplex(0.0,0.0), CumulativeF3iDzF3iMinusF4iDzF4iClover = TComplex(0.0,0.0);
  TComplex SumF3iDzF4iPlusF4iDzF3iClover  = TComplex(0.0,0.0), CumulativeF3iDzF4iPlusF4iDzF3iClover  = TComplex(0.0,0.0);

  TComplex SumF3iF3iMinusF4iF4iCloverTraceless     = TComplex(0.0,0.0), CumulativeF3iF3iMinusF4iF4iCloverTraceless     = TComplex(0.0,0.0);
  TComplex SumF3iF4iPlusF4iF3iCloverTraceless      = TComplex(0.0,0.0), CumulativeF3iF4iPlusF4iF3iCloverTraceless      = TComplex(0.0,0.0);
  TComplex SumF3iDzF3iMinusF4iDzF4iCloverTraceless = TComplex(0.0,0.0), CumulativeF3iDzF3iMinusF4iDzF4iCloverTraceless = TComplex(0.0,0.0);
  TComplex SumF3iDzF4iPlusF4iDzF3iCloverTraceless  = TComplex(0.0,0.0), CumulativeF3iDzF4iPlusF4iDzF3iCloverTraceless  = TComplex(0.0,0.0);
  
  f1<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;
  f2<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;
  f3<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;
  f4<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;
  f5<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;
  f6<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;
  f7<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;
  f8<<"# Nt="<<LC.Nt<<"\t Ns="<<LC.Nx<<"\t, SU3 calculation for Beta ="<<LC.beta<<endl;

  
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
	  LC.LinksContent[i][m][m]=TComplex(1.0,0.0); 
	}      
    }

  
  //Constructing first stable gauge field configuration
  LC.Iteration=0.0;
  int MeasurementCount=0.0;
  
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
      LC.Iteration = LC.Iteration + 1;
      MeasurementCount = MeasurementCount +1;
      if(Iter==200)
	{
	  MeasurementCount=1;
	  SumF3iF3iMinusF4iF4i     = TComplex(0.0,0.0), SumF3iF3iMinusF4iF4iTraceless     = TComplex(0.0,0.0);
	  SumF3iF4iPlusF4iF3i      = TComplex(0.0,0.0), SumF3iF4iPlusF4iF3iTraceless      = TComplex(0.0,0.0);
	  SumF3iDzF3iMinusF4iDzF4i = TComplex(0.0,0.0), SumF3iDzF3iMinusF4iDzF4iTraceless = TComplex(0.0,0.0);
	  SumF3iDzF4iPlusF4iDzF3i  = TComplex(0.0,0.0), SumF3iDzF4iPlusF4iDzF3iTraceless  = TComplex(0.0,0.0);

	  SumF3iF3iMinusF4iF4iClover     = TComplex(0.0,0.0), SumF3iF3iMinusF4iF4iCloverTraceless     = TComplex(0.0,0.0);
	  SumF3iF4iPlusF4iF3iClover      = TComplex(0.0,0.0), SumF3iF4iPlusF4iF3iCloverTraceless      = TComplex(0.0,0.0);
	  SumF3iDzF3iMinusF4iDzF4iClover = TComplex(0.0,0.0), SumF3iDzF3iMinusF4iDzF4iCloverTraceless = TComplex(0.0,0.0);
	  SumF3iDzF4iPlusF4iDzF3iClover  = TComplex(0.0,0.0), SumF3iDzF4iPlusF4iDzF3iCloverTraceless  = TComplex(0.0,0.0);
	}
      cout<<"\nComputing operators for Configuration "<<LC.Iteration<<endl;                  
      //Evaluate Trace of FieldStrength-ZdirectedCovariantDerivative-FieldStrength Matrix	       		      
      //cout<<"Going to calculate Fmunu using single Plaquette \n\n";		      
      LC.ComputeFmunuUsingSinglePlaquette();
      TraceF3iF3iMinusF4iF4i     = LC.TraceF3iF3iMinusF4iF4i;
      TraceF3iF4iPlusF4iF3i      = LC.TraceF3iF4iPlusF4iF3i;
      TraceF3iDzF3iMinusF4iDzF4i = LC.TraceF3iDzF3iMinusF4iDzF4i;
      TraceF3iDzF4iPlusF4iDzF3i  = LC.TraceF3iDzF4iPlusF4iDzF3i;

      TraceF3iF3iMinusF4iF4iTraceless     = LC.TraceF3iF3iMinusF4iF4iTraceless;
      TraceF3iF4iPlusF4iF3iTraceless      = LC.TraceF3iF4iPlusF4iF3iTraceless;
      TraceF3iDzF3iMinusF4iDzF4iTraceless = LC.TraceF3iDzF3iMinusF4iDzF4iTraceless;
      TraceF3iDzF4iPlusF4iDzF3iTraceless  = LC.TraceF3iDzF4iPlusF4iDzF3iTraceless;
      
      SumF3iF3iMinusF4iF4i            = SumF3iF3iMinusF4iF4i     + TraceF3iF3iMinusF4iF4i;
      SumF3iF4iPlusF4iF3i             = SumF3iF4iPlusF4iF3i      + TraceF3iF4iPlusF4iF3i;
      SumF3iDzF3iMinusF4iDzF4i        = SumF3iDzF3iMinusF4iDzF4i + TraceF3iDzF3iMinusF4iDzF4i;
      SumF3iDzF4iPlusF4iDzF3i         = SumF3iDzF4iPlusF4iDzF3i  + TraceF3iDzF4iPlusF4iDzF3i;
      
      CumulativeF3iF3iMinusF4iF4i     = SumF3iF3iMinusF4iF4i/(MeasurementCount*1.0);      
      CumulativeF3iF4iPlusF4iF3i      = SumF3iF4iPlusF4iF3i/(MeasurementCount*1.0);      
      CumulativeF3iDzF3iMinusF4iDzF4i = SumF3iDzF3iMinusF4iDzF4i/(MeasurementCount*1.0);      
      CumulativeF3iDzF4iPlusF4iDzF3i  = SumF3iDzF4iPlusF4iDzF3i/(MeasurementCount*1.0);

      SumF3iF3iMinusF4iF4iTraceless            = SumF3iF3iMinusF4iF4iTraceless     + TraceF3iF3iMinusF4iF4iTraceless;
      SumF3iF4iPlusF4iF3iTraceless             = SumF3iF4iPlusF4iF3iTraceless      + TraceF3iF4iPlusF4iF3iTraceless;
      SumF3iDzF3iMinusF4iDzF4iTraceless        = SumF3iDzF3iMinusF4iDzF4iTraceless + TraceF3iDzF3iMinusF4iDzF4iTraceless;
      SumF3iDzF4iPlusF4iDzF3iTraceless         = SumF3iDzF4iPlusF4iDzF3iTraceless  + TraceF3iDzF4iPlusF4iDzF3iTraceless;

      CumulativeF3iF3iMinusF4iF4iTraceless     = SumF3iF3iMinusF4iF4iTraceless/(MeasurementCount*1.0);
      CumulativeF3iF4iPlusF4iF3iTraceless      = SumF3iF4iPlusF4iF3iTraceless/(MeasurementCount*1.0);
      CumulativeF3iDzF3iMinusF4iDzF4iTraceless = SumF3iDzF3iMinusF4iDzF4iTraceless/(MeasurementCount*1.0);
      CumulativeF3iDzF4iPlusF4iDzF3iTraceless  = SumF3iDzF4iPlusF4iDzF3iTraceless/(MeasurementCount*1.0);

	//Computing Fmunu using Clover method
      LC.ComputeFmunuUsingFourPlaquette();
      TraceF3iF3iMinusF4iF4iClover     = LC.TraceF3iF3iMinusF4iF4i;
      TraceF3iF4iPlusF4iF3iClover      = LC.TraceF3iF4iPlusF4iF3i;
      TraceF3iDzF3iMinusF4iDzF4iClover = LC.TraceF3iDzF3iMinusF4iDzF4i;
      TraceF3iDzF4iPlusF4iDzF3iClover  = LC.TraceF3iDzF4iPlusF4iDzF3i;

      TraceF3iF3iMinusF4iF4iCloverTraceless     = LC.TraceF3iF3iMinusF4iF4iTraceless;
      TraceF3iF4iPlusF4iF3iCloverTraceless      = LC.TraceF3iF4iPlusF4iF3iTraceless;
      TraceF3iDzF3iMinusF4iDzF4iCloverTraceless = LC.TraceF3iDzF3iMinusF4iDzF4iTraceless;
      TraceF3iDzF4iPlusF4iDzF3iCloverTraceless  = LC.TraceF3iDzF4iPlusF4iDzF3iTraceless;
      
      SumF3iF3iMinusF4iF4iClover            = SumF3iF3iMinusF4iF4iClover     + TraceF3iF3iMinusF4iF4iClover;
      SumF3iF4iPlusF4iF3iClover             = SumF3iF4iPlusF4iF3iClover      + TraceF3iF4iPlusF4iF3iClover;
      SumF3iDzF3iMinusF4iDzF4iClover        = SumF3iDzF3iMinusF4iDzF4iClover + TraceF3iDzF3iMinusF4iDzF4iClover;
      SumF3iDzF4iPlusF4iDzF3iClover         = SumF3iDzF4iPlusF4iDzF3iClover  + TraceF3iDzF4iPlusF4iDzF3iClover;
      
      CumulativeF3iF3iMinusF4iF4iClover     = SumF3iF3iMinusF4iF4iClover/(MeasurementCount*1.0);      
      CumulativeF3iF4iPlusF4iF3iClover      = SumF3iF4iPlusF4iF3iClover/(MeasurementCount*1.0);      
      CumulativeF3iDzF3iMinusF4iDzF4iClover = SumF3iDzF3iMinusF4iDzF4iClover/(MeasurementCount*1.0);      
      CumulativeF3iDzF4iPlusF4iDzF3iClover  = SumF3iDzF4iPlusF4iDzF3iClover/(MeasurementCount*1.0);

      SumF3iF3iMinusF4iF4iCloverTraceless            = SumF3iF3iMinusF4iF4iCloverTraceless     + TraceF3iF3iMinusF4iF4iCloverTraceless;
      SumF3iF4iPlusF4iF3iCloverTraceless             = SumF3iF4iPlusF4iF3iCloverTraceless      + TraceF3iF4iPlusF4iF3iCloverTraceless;
      SumF3iDzF3iMinusF4iDzF4iCloverTraceless        = SumF3iDzF3iMinusF4iDzF4iCloverTraceless + TraceF3iDzF3iMinusF4iDzF4iCloverTraceless;
      SumF3iDzF4iPlusF4iDzF3iCloverTraceless         = SumF3iDzF4iPlusF4iDzF3iCloverTraceless  + TraceF3iDzF4iPlusF4iDzF3iCloverTraceless;

      CumulativeF3iF3iMinusF4iF4iCloverTraceless     = SumF3iF3iMinusF4iF4iCloverTraceless/(MeasurementCount*1.0);
      CumulativeF3iF4iPlusF4iF3iCloverTraceless      = SumF3iF4iPlusF4iF3iCloverTraceless/(MeasurementCount*1.0);
      CumulativeF3iDzF3iMinusF4iDzF4iCloverTraceless = SumF3iDzF3iMinusF4iDzF4iCloverTraceless/(MeasurementCount*1.0);
      CumulativeF3iDzF4iPlusF4iDzF3iCloverTraceless  = SumF3iDzF4iPlusF4iDzF3iCloverTraceless/(MeasurementCount*1.0);
      cout<<"For current configuration: Lattice points = "<<LC.LatticePoints<<endl;
      cout<<"Tr(F3iF3i-F4iF4i)    = "<<TraceF3iF3iMinusF4iF4i<<"\t Cumulative value ="<<CumulativeF3iF3iMinusF4iF4i<<endl;
      cout<<"Tr(F3iF3i-F4iF4i) TrL= "<<TraceF3iF3iMinusF4iF4iTraceless<<"\t Cumulative value ="<<CumulativeF3iF3iMinusF4iF4iTraceless<<endl;
      cout<<"Tr(F3iF3i-F4iF4i)CLOV= "<<TraceF3iF3iMinusF4iF4iClover<<"\t Cumulative value ="<<CumulativeF3iF3iMinusF4iF4iClover<<endl;
      cout<<"Tr(F3iF3i-F4iF4i) TrL= "<<TraceF3iF3iMinusF4iF4iCloverTraceless<<"\t Cumulative value ="<<CumulativeF3iF3iMinusF4iF4iCloverTraceless<<endl;

      cout<<"Tr(F3iF4i+F4iF3i)    = "<<TraceF3iF4iPlusF4iF3i<<"\t Cumulative value ="<<CumulativeF3iF4iPlusF4iF3i<<endl;
      cout<<"Tr(F3iF4i+F4iF3i) TrL= "<<TraceF3iF4iPlusF4iF3iTraceless<<"\t Cumulative value ="<<CumulativeF3iF4iPlusF4iF3iTraceless<<endl;
      cout<<"Tr(F3iF4i+F4iF3i)CLOV= "<<TraceF3iF4iPlusF4iF3iClover<<"\t Cumulative value ="<<CumulativeF3iF4iPlusF4iF3iClover<<endl;
      cout<<"Tr(F3iF4i+F4iF3i) TrL= "<<TraceF3iF4iPlusF4iF3iCloverTraceless<<"\t Cumulative value ="<<CumulativeF3iF4iPlusF4iF3iCloverTraceless<<endl;

      cout<<"Tr(F3iDzF3i-F4iDzF4i)    = "<<TraceF3iDzF3iMinusF4iDzF4i<<"\t Cumulative value ="<<CumulativeF3iDzF3iMinusF4iDzF4i <<endl;
      cout<<"Tr(F3iDzF3i-F4iDzF4i) TrL= "<<TraceF3iDzF3iMinusF4iDzF4iTraceless<<"\t Cumulative value ="<<CumulativeF3iDzF3iMinusF4iDzF4iTraceless <<endl;
      cout<<"Tr(F3iDzF3i-F4iDzF4i)CLOV= "<<TraceF3iDzF3iMinusF4iDzF4iClover<<"\t Cumulative value ="<<CumulativeF3iDzF3iMinusF4iDzF4iClover <<endl;
      cout<<"Tr(F3iDzF3i-F4iDzF4i) TrL= "<<TraceF3iDzF3iMinusF4iDzF4iCloverTraceless<<"\t Cumulative value ="<<CumulativeF3iDzF3iMinusF4iDzF4iCloverTraceless <<endl;

      cout<<"Tr(F3iDzF4i+F4iDzF3i)    = "<<TraceF3iDzF4iPlusF4iDzF3i<<"\t Cumulative value ="<<CumulativeF3iDzF4iPlusF4iDzF3i<<endl;
      cout<<"Tr(F3iDzF4i+F4iDzF3i) TrL= "<<TraceF3iDzF4iPlusF4iDzF3iTraceless<<"\t Cumulative value ="<<CumulativeF3iDzF4iPlusF4iDzF3iTraceless<<endl;
      cout<<"Tr(F3iDzF4i+F4iDzF3i)CLOV= "<<TraceF3iDzF4iPlusF4iDzF3iClover<<"\t Cumulative value ="<<CumulativeF3iDzF4iPlusF4iDzF3iClover<<endl;
      cout<<"Tr(F3iDzF4i+F4iDzF3i) TrL= "<<TraceF3iDzF4iPlusF4iDzF3iCloverTraceless<<"\t Cumulative value ="<<CumulativeF3iDzF4iPlusF4iDzF3iCloverTraceless<<endl;

      f1<<LC.Iteration<<"\t"<<TraceF3iF3iMinusF4iF4i.Re()<<"\t"<<TraceF3iF3iMinusF4iF4i.Im()<<"\t"<<CumulativeF3iF3iMinusF4iF4i.Re()<<"\t"<<CumulativeF3iF3iMinusF4iF4i.Im()<<"\t"<<TraceF3iF4iPlusF4iF3i.Re()<<"\t"<<TraceF3iF4iPlusF4iF3i.Im()<<"\t"<<CumulativeF3iF4iPlusF4iF3i.Re()<<"\t"<<CumulativeF3iF4iPlusF4iF3i.Im()<<endl;
      f2<<LC.Iteration<<"\t"<<TraceF3iDzF3iMinusF4iDzF4i.Re()<<"\t"<<TraceF3iDzF3iMinusF4iDzF4i.Im()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4i.Re()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4i.Im()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3i.Re()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3i.Im()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3i.Re()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3i.Im()<<endl;

      f3<<LC.Iteration<<"\t"<<TraceF3iF3iMinusF4iF4iTraceless.Re()<<"\t"<<TraceF3iF3iMinusF4iF4iTraceless.Im()<<"\t"<<CumulativeF3iF3iMinusF4iF4iTraceless.Re()<<"\t"<<CumulativeF3iF3iMinusF4iF4iTraceless.Im()<<"\t"<<TraceF3iF4iPlusF4iF3iTraceless.Re()<<"\t"<<TraceF3iF4iPlusF4iF3iTraceless.Im()<<"\t"<<CumulativeF3iF4iPlusF4iF3iTraceless.Re()<<"\t"<<CumulativeF3iF4iPlusF4iF3iTraceless.Im()<<endl;
      f4<<LC.Iteration<<"\t"<<TraceF3iDzF3iMinusF4iDzF4iTraceless.Re()<<"\t"<<TraceF3iDzF3iMinusF4iDzF4iTraceless.Im()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4iTraceless.Re()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4iTraceless.Im()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3iTraceless.Re()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3iTraceless.Im()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3iTraceless.Re()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3iTraceless.Im()<<endl;
      
      //Colver based calculation of Fmunu
      f5<<LC.Iteration<<"\t"<<TraceF3iF3iMinusF4iF4iClover.Re()<<"\t"<<TraceF3iF3iMinusF4iF4iClover.Im()<<"\t"<<CumulativeF3iF3iMinusF4iF4iClover.Re()<<"\t"<<CumulativeF3iF3iMinusF4iF4iClover.Im()<<"\t"<<TraceF3iF4iPlusF4iF3iClover.Re()<<"\t"<<TraceF3iF4iPlusF4iF3iClover.Im()<<"\t"<<CumulativeF3iF4iPlusF4iF3iClover.Re()<<"\t"<<CumulativeF3iF4iPlusF4iF3iClover.Im()<<endl;
      f6<<LC.Iteration<<"\t"<<TraceF3iDzF3iMinusF4iDzF4iClover.Re()<<"\t"<<TraceF3iDzF3iMinusF4iDzF4iClover.Im()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4iClover.Re()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4iClover.Im()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3iClover.Re()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3iClover.Im()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3iClover.Re()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3iClover.Im()<<endl;

      f7<<LC.Iteration<<"\t"<<TraceF3iF3iMinusF4iF4iCloverTraceless.Re()<<"\t"<<TraceF3iF3iMinusF4iF4iCloverTraceless.Im()<<"\t"<<CumulativeF3iF3iMinusF4iF4iCloverTraceless.Re()<<"\t"<<CumulativeF3iF3iMinusF4iF4iCloverTraceless.Im()<<"\t"<<TraceF3iF4iPlusF4iF3iCloverTraceless.Re()<<"\t"<<TraceF3iF4iPlusF4iF3iCloverTraceless.Im()<<"\t"<<CumulativeF3iF4iPlusF4iF3iCloverTraceless.Re()<<"\t"<<CumulativeF3iF4iPlusF4iF3iCloverTraceless.Im()<<endl;
      f8<<LC.Iteration<<"\t"<<TraceF3iDzF3iMinusF4iDzF4iCloverTraceless.Re()<<"\t"<<TraceF3iDzF3iMinusF4iDzF4iCloverTraceless.Im()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4iCloverTraceless.Re()<<"\t"<<CumulativeF3iDzF3iMinusF4iDzF4iCloverTraceless.Im()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3iCloverTraceless.Re()<<"\t"<<TraceF3iDzF4iPlusF4iDzF3iCloverTraceless.Im()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3iCloverTraceless.Re()<<"\t"<<CumulativeF3iDzF4iPlusF4iDzF3iCloverTraceless.Im()<<endl;

      cout<<"\nGenerating Configuration "<<LC.Iteration+1<<endl;
      //Touching each internal link with heat bath
      for(Int_t ic = 0 ; ic< LC.TotalLinks ; ic++)
	{
	  //cout<<"Processing Link"<<ic<<endl;
	  LC.UpdateLink(ic); //Apply heat bath algorithm to each internal link and apply periodic boundary condition
	}
      
    } 

      int EndTime = time(NULL);
      int Hour = (EndTime-StartTime)/3600;
      int Minute = ((EndTime-StartTime)/60)-Hour*60;
      int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
      cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
      f0<<"For Beta = "<<LC.beta<<"\t Programme run time H::M::S = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl; f0.close();
      f1.close(); f2.close(); f3.close(); f4.close(); f5.close(); f6.close(); f7.close(); f8.close();
      return 0;      
}
