// Generare numeri random con distribuzione gaussiana in tutto R secondo la f(x) data. Usare
// il metodo esatto e quello di accept/reject con la funzione g(x) data 
// verifica che la distribuzione sia la medesima

#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#define FUNC(x) ((*func)(x))
#define nMAX 1e4

using namespace std;

double rand_range (double min, double max)
{
  return min + (max - min) * rand () / static_cast<double> (RAND_MAX);
}
////////////////////////////////////////////////////////// METODO ESATTO    
double func(double x) // f(x)
{
    return (1/(sqrt(M_PI)))*exp(-1.*pow(x,2));
}

//Normalizziamo la func a 1

double funcr (double z) // genera r
{   
    double zp = 1 - z;
    return sqrt(-log(zp));
}

double funcx (double r, double a) // per generare x
{
    return r*cos(2*M_PI*a);
}
///////////////////////////////////////////////////// ACCEPT/REJECT

double g(double x){
    double A = (1./sqrt(M_PI));
        if(x < -1.)
    {
        return -A* x * exp(1-pow(x,2));
    }

    else if(-1. <= x && x <= 1.)
    {
         return A;   
    }

        else if(x > 1.){
        return A * x * exp(1-pow(x,2));
    }
}

double h (double x){ //g(x) normalizzata come pdf
    if(x < -1.)
    {
        return -1./3. * x * exp(1-pow(x,2));
    }

    else if(-1. <= x && x <= 1.)
    {
         return 1./3.;   
    }

        else if(x > 1.){
        return 1./3. * x * exp(1-pow(x,2));
    }


}

double winv (double y) { //cumulativa inversa
    if(0 <= y && y < 1./6.)
    {
        return -1.*sqrt(1-log(6*y));
    }

    else if(1./6. <= y && y <= 5./6.)
    {
         return 3*y - 3./2.;   
    }

    else if(5./6. < y && y <= 1){
        return sqrt(1-log(6-6*y));
    }
}


int main ()
{
        // METODO ESATTO
    FILE * importancesampling = fopen("importancesampling.txt", "w");    
    double r = 0.;
    double xs = 0.;
    for(int i = 1; i <= nMAX + 1; i = i + 25    ){
        r = funcr(rand_range(0.,1.));
        xs = funcx(r, rand_range(0., 1.));
        fprintf(importancesampling, "%.16lf\n", xs);
    }

        // METODO A/R

    FILE * accrej3 = fopen("accept3.txt", "w");
    FILE * minimainterazione = fopen("minimint.txt", "w");
    double A = (1./sqrt(M_PI));
    int j = 0;

    for(int i = 0; i < nMAX; i++  ){
        double W = winv(rand_range(0.,1.)); //campioniamo un numero secondo l'inversa della cumulativa di h(x) 
        double U = rand_range(0.,1.);
        //applico il A/R
        if((func(W))/(3*A*h(W)) >= U){ // la costante che mi rende h maggiorante è proprio N. In questo modo 3*A*h(x) > f(x) for all x
            fprintf(accrej3, "%.16lf\n", W);
            fprintf(minimainterazione, "%d\n", i-j);
            j=i; //verifichiamo se il numero di interazioni minime è 1/c, ovvero 1/3A
        }

    }
    return 0;
}

