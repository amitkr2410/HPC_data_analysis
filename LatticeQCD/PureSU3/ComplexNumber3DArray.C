#include"iostream"
#include"TMath.h"
#include"TComplex.h"
#include"cmath"

using namespace std;

void ComplexNumber3DArray()

{
  int TotalLinks=2;
  int* LinksDirection;
  int** LinksCoordinate;
  TComplex*** A;//A[2][3][3];

  LinksDirection = new int [TotalLinks];
  LinksCoordinate = new int* [TotalLinks];
  A = new TComplex** [2];


for(int i=0; i<TotalLinks; i++)
  { 
    A[i] = new TComplex* [3];
    LinksCoordinate[i] = new int [4];

  }
  
 for(int i=0; i<2; i++)
   {
     for(int j=0; j<3; j++)
       {
	 A[i][j] = new TComplex [3];
       }
   }
 
 
 for(int i=0; i<2; i++)
   {
     
     for(int j=0; j<3; j++)
       {
	 
	 for(int k=0; k<3; k++)
	   {
	     A[i][j][k]= i+j+k;
	     cout<<A[i][j][k]<<"\t";
	   }
	 cout<<"\n";
       }
     cout<<"\n\n";
   }

}
