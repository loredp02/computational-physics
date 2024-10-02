#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#define FUNC(x) ((*func)(x))
#define NMAX 100 // numero dei dati sulla quale medio, che varia.
#define M 1e5 // taglia del pacchetto mediato che NON varia

using namespace std;

FILE * F_N = fopen("N.txt", "w");

//Distrib unif
FILE * F_ux1 = fopen("ux1.txt", "w");
FILE * F_ux2 = fopen("ux2.txt", "w");
FILE * F_ux3 = fopen("ux3.txt", "w");
FILE * F_uDelta1 = fopen("udelta1.txt", "w");
FILE * F_uDelta2 = fopen("udelta2.txt", "w");
//Distrib discr.
FILE * F_dx1 = fopen("dx1.txt", "w");
FILE * F_dx2 = fopen("dx2.txt", "w");
FILE * F_dx3 = fopen("dx3.txt", "w");
FILE * F_dDelta1 = fopen("ddelta1.txt", "w");
FILE * F_dDelta2 = fopen("ddelta2.txt", "w");



// distribuzione uniforme tra (-1,1)
double rand_range (double min, double max)
{
  return min + (max - min) * rand () / static_cast<double> (RAND_MAX);
}

// distribuzione discreta o -1 o 1
double discrete_dist (){
    double u;
    double d;
    u = rand_range(0.,1.);
    if(u <= 0.5)
    {
        d = -1.;
    }

    else{
        d = 1.;
    }

    return d;
}

// calcola <x^n>
double momento (vector<double> campione, int exp, int N)
{
    double somma = 0.;
    for(int i = 0; i < N; i++)
    {
        somma += pow(campione[i],exp); 
    }

    return (somma)/(double)N;
}

int main ()
{   
    for(int N = 1; N <= NMAX + 1; N = N + 5){ //varia il numero di numeri generati randomicamente fino ad un limite NMAX
        vector<double> xi(M); //pacchetto mediato
        vector<double> d_xi(M);
        

        for(int p = 0; p < M; p++){ // ne voglio generare M

            vector<double> r(N); 
            vector<double> d_r(N);
            

            for(int j = 0; j < N; j++) //riempiamo il vettore con N numeri randomici
            {
                r[j] = rand_range(-1.,1.);
                d_r[j] = discrete_dist();
                //printf("r_ %d = %lf\n", j, r[j]);
            }

        // a questo punto ho il vettore ri con N numeri generati randomicamente, mettere la media con N 
            xi[p] = momento(r, 1, N);
            d_xi[p] = momento(d_r, 1, N);
        // possiamo cancellare gli ri tanto è già stata fatta la media
             r.clear();
             d_r.clear();
            //printf("x_ %d = %lf\n", p, xi[p]);

        } // qui finisce il ciclo in p

        double s1 = momento(xi, 1, M);
        double s2 = momento(xi, 2, M);
        double s3 = momento(xi, 3, M);
        double s4 = momento(xi, 4, M);
        double s6 = momento(xi, 6, M);

        double d_s1 = momento(d_xi, 1, M);
        double d_s2 = momento(d_xi, 2, M);
        double d_s3 = momento(d_xi, 3, M);
        double d_s4 = momento(d_xi, 4, M);
        double d_s6 = momento(d_xi, 6, M);

       /* printf("DISTRIB UNIF: Per N = %d, si ha che \n s1 = %.14lf, \n s2 = %.14lf, \n s4 - 3s2^2 = %.16lf,\n s6 - 15s2^3 = %.16lf\n \n", N, s1, s2, (s4 - 3*pow(s2,2)), (s6 - 15*pow(s2,3)));
        printf("DISTRIB DISCRETA: Per N = %d, si ha che \n s1 = %.14lf, \n s2 = %.14lf, \n s4 - 3s2^2 = %.16lf,\n s6 - 15s2^3 = %.16lf\n \n", N, d_s1, d_s2, (d_s4 - 3*pow(d_s2,2)), (d_s6 - 15*pow(d_s2,3)));    
        printf("%lf\n", s4/(pow(s2,2))); */

        // salviamo i numeri nei file txt

        fprintf(F_N, "%d\n", N);
        
        //distrib uniforme
        fprintf(F_ux1, "%.16lf\n", s1);
        fprintf(F_ux2, "%.16lf\n", s2);
        fprintf(F_ux3, "%.16lf\n", s3);
        fprintf(F_uDelta1, "%.16lf\n", (s4 - 3*pow(s2,2)));
        fprintf(F_uDelta2, "%.16lf\n", (s6 - 15*pow(s2,3)));

        //distrib discreta
        fprintf(F_dx1, "%.16lf\n", d_s1);
        fprintf(F_dx2, "%.16lf\n", d_s2);
        fprintf(F_dx3, "%.16lf\n", d_s3);
        fprintf(F_dDelta1, "%.16lf\n", (d_s4 - 3*pow(d_s2,2)));
        fprintf(F_dDelta2, "%.16lf\n", (d_s6 - 15*pow(d_s2,3)));

        
        xi.clear();
    
    }



    return 0;
}