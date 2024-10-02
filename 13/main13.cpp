#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#define FUNC(x) ((*func)(x))


using namespace std;

double rand_range (double min, double max)
{
  return min + (max - min) * rand () / static_cast<double> (RAND_MAX);
}

double funcA(double x)
{
    return asin(x);
}

double funcB(double x)
{
    double z = -1 + 2*x;
    return 2*acos(z);
}



int main ()
{
    int nMAX = 100;
    double yA;
    double yB;
    double sumA;
    double sumB;
    double RisA = 0.57080; // risultato del primo integrale 
    double RisB = M_PI; // risultato del secondo integrale
    FILE * F_dev1 = fopen("devs1.txt", "w");
    FILE * F_dev2 = fopen("devs2.txt", "w");    
    for(int j = 0; j < 10; j++){

        nMAX = nMAX + 1e3;
    
        for(int i = 1; i <= nMAX + 1; i++){
        
        yA = rand_range(0.,1.);
        yB = rand_range(0.,1.);
        sumA = sumA + funcA(yA);
        sumB = sumB + funcB(yB);
    
        }

    double totA = sumA/nMAX;
    double totB = sumB/nMAX;
    sumA = 0;
    sumB = 0;
    fprintf(F_dev1, "%d %.16lf\n", nMAX, abs(RisA-totA));
    fprintf(F_dev2, "%d %.16lf\n", nMAX, abs(RisB-totB));
    //printf("Il primo integrale ha valore: %f, calcolato con l'i.s. viene: %f, con una differenza %f \n", RisA, totA, abs(RisA-totA));
    //printf("Il secondo integrale ha valore: %f, calcolato con l'i.s. viene: %f, con una differenza %f \n", RisB, totB, abs(RisB-totB));

    }
    return 0;
}