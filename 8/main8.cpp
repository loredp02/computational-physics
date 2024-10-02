/*Calcolare il volume SM della sfera unitaria in M dimensioni, con la formula analitica data.
Calcolare in funzione di M con 20 punti di integrazione per direzione coordinata,
l’integrale deterministico (usare la formula del midpoint) e l’integrale Monte Carlo
usando e5 punti. Monitorare anche i tempi di calcolo.
*/

#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#include <bits/stdc++.h>
#define nPoints 20  
#define nMCarlo 1e05

using namespace std;

// file per salvare i dati
FILE * DET_TIME = fopen("dettime.txt", "w");
FILE * MC_TIME = fopen("mctime.txt", "w");

// funzione fattoriale 
int factorial(int n)
{
    if (n == 0 || n == 1) {
        return 1;
    }
    
    return n * factorial(n-1);
}

// funzione fattoriale doppio
int doublefactorial(int n)
{
    if (n == 0 || n==1)
      return 1;
    return n * doublefactorial(n-2);
}

// funzione per generare numeri pseudorandomici tra due estremi con pdf uniforme

double rand_range (double min, double max)
{
  return min + (max - min) * rand () / static_cast<double> (RAND_MAX);
}

// funzione per calcolare il valore esatto del volume dell'ipersfera in M Dimensioni
void Sm_Exact (int M){
    double Sm;
    if(M % 2 == 0) {
        Sm = (pow(M_PI, M/2))/(factorial(M/2));
    }

    else {
        Sm = (pow(2, (M+1)/2) * pow(M_PI, (M-1)/2))/(doublefactorial(M));
    }

    printf("%lf     ", Sm);

}

// funzione per calcolare il valore del volume dell'ipersfera in M Dimensioni in modo deterministico
void Sm_Det (int M)
{
    clock_t start, end;
    double tempo;
    start = clock();
    long int npt, dim, n, punto, x_i;
    double norm, volume, h;
    long int n_axis = (nPoints/2);
        


    npt = 1;
    for(int d=1; d <= (M-1); d++) npt = npt * n_axis;
    
    h = 1./(double)n_axis;

    for(int n = 0; n < npt; n++){
        punto = n;

        norm = 0.;
        dim = npt / n_axis;

        for(int d=1; d<=(M-1); d++){
            x_i = punto/dim;
            norm = norm + pow(h*(double)x_i + h/2., 2);

            punto = punto - x_i * dim;
            dim = dim/n_axis;
        }

        if(norm <= 1.) volume = volume + sqrt(1.-norm);
    }

    volume = pow(2,M) * volume/(double)npt;
    end = clock();
    tempo = (double) (end-start)/(double)(CLOCKS_PER_SEC);
    printf("%lf %lf     ", tempo, volume);
    fprintf(DET_TIME, "%.16lf\n", tempo);
    
}




// funzione per calcolare il valore del volume dell'ipersfera in M Dimensioni in modo pseudo-randomico
void Sm_Random (int M){
    clock_t start, end;
    double tempo;
    start = clock();
    double x_i;
    double norm, volume;

    volume = 0;

    for(int n = 0; n < nMCarlo; n++){
        norm = 0.;
        for(int d=1; d<=(M-1); d++){
            x_i = rand_range(0,1);
            norm = norm + pow(x_i,2);
        }

        if(norm <= 1.) volume = volume + sqrt(1.-norm); //f(xi)
    }

    volume = pow(2,M) * volume/(double)nMCarlo;
    end = clock();
    tempo = (double) (end-start)/(double)(CLOCKS_PER_SEC);
    printf("%lf %lf     \n", tempo, volume);
    fprintf(MC_TIME, "%.16lf\n", tempo);


    
}

int main ()
{
    int M = 9;
    int k = 0;
    printf("M       esatto      deterministico          montecarlo\n");
    printf("                    t(s)     val            t(s)     val\n");
    printf("----------------------------------------------------------\n");
    for(int k = 1; k <= M; k++){
    printf("%d      ", k);
    

    Sm_Exact(k);
    Sm_Det(k);
    Sm_Random(k);
    }


    return 0;
}