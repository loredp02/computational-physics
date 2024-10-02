#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>
#include <cmath>
#define _USE_MATH_DEFINES


using namespace std;

int main (int argc, char ** argv)
{
  FILE * SINGLE_DIRECT = fopen("SINGLE_DIRECT.txt", "w");
  FILE * DOUBLE_DIRECT = fopen("DOUBLE_DIRECT.txt", "w");
  FILE * SINGLE_INVERSE = fopen("SINGLE_INVERSE.txt", "w");
  FILE * DOUBLE_INVERSE = fopen("DOUBLE_INVERSE.txt", "w");

  int i = 0;
  int j = 0;
  
  int K = 60000; // ripetizioni
  int h = 200; // incremento 
  int N = 0;

  float fexact_val = (M_PI * M_PI)/6; //valore esatto della sommatoria in precisione singola
  double dexact_val = (M_PI * M_PI)/6; //valore esatto della sommatoria in precisione doppia

  //somma diretta 

      float  single_directsum = 0.; //precisione singola
      double double_directsum = 0.; //precisione doppia

  for(j = 1; j <= K; j = j + h){
      single_directsum = 0.;
      double_directsum = 0.;
      N = j;
      for(int i = 1; i <= N; i++)
    {
      single_directsum += 1./(i*i);
      double_directsum += 1./(i*i);
    } 
    fprintf(SINGLE_DIRECT, "%d %f\n", j, (fexact_val - single_directsum));
    fprintf(DOUBLE_DIRECT, "%d %f\n", j, (dexact_val - double_directsum));

  }
     
  
      cout <<  "somma diretta in precisione singola: " << single_directsum << "\n" << endl;
      cout <<  "somma diretta in doppia precisione: " << double_directsum << "\n" << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //somma inversa 

      float single_inversesum = 0.;
      double double_inversesum = 0.;


  for(j = 1; j <= K; j = j + h){
      single_inversesum = 0.;
      double_inversesum = 0.;
      N = j;
      for(int i = N; i > 0; i--)
	{
	  single_inversesum += 1./(i*i);
	  double_inversesum += 1./(i*i);
	}
    fprintf(SINGLE_INVERSE, "%d %f\n", j, (fexact_val - single_inversesum)); 
    fprintf(DOUBLE_INVERSE, "%d %f\n", j, (dexact_val - double_inversesum));
  }

      cout << "somma inversa in precisione singola: " << single_inversesum << "\n" << endl; 
      cout << "somma inversa in doppia precisione: " << double_inversesum << "\n" << endl;
  


  return 0;
}
