#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
//includiamo una libreria per operare con i numeri complessi
#include <complex>
#define EPS 1E-9  // precisione max
using namespace std;


complex<double> Funzione (complex<double> & z)
{
    //printf("Il numero complesso z passato in Funzione ha Im = (%.lf)\n", z.imag());
    complex<double> p (-1., 0.);
    complex<double> y = pow(z,3) + p;
    //printf("La Funzione restituisce y con Re (%.lf) e Im (%.lf)\n", y.real(), y.imag());
    return y;
}

complex<double> DerivFunzione (complex<double> & z)
{   
    //printf("Il numero complesso z passato in DerivFunzione ha Im = (%.lf)\n", z.imag());
    complex<double> p (3., 0.);
    complex<double> y = p * pow(z,2);
    //printf("La DerivFunzione restituisce y con Re (%.lf) e Im (%.lf)\n", y.real(), y.imag());
    return y;
}


int main ()
{
    double lowExtreme = -2.;
    double upExtreme = 2.;
    double leftExtreme = -2.;
    double rightExtreme = 2.;

    double fraction = 0.25;

    int NRows = (int) ( (double) abs(rightExtreme-leftExtreme))/(fraction);
    int NCol = (int) ( (double) abs(lowExtreme-upExtreme))/(fraction);

    int K = 1E3; // numero massimo di iterazioni

    FILE * file_griglia = fopen("grigliafiles.txt", "w"); //per la griglia
    FILE * file_curva = fopen("curvefiles.txt", "w");  //andiamo a depositare gli step che compie da un punto di partenza (-1.25,0.25) 

    // c'è da iterare con punto di partenza ogni z della griglia, dunque a partire da z = -2 - 2i e andando a coprire tutta la griglia.

    int i = 0; // righe
    int j = 0; // colonne
    
    for(int i = 0; i <= NRows; i++) // ciclo sulle righe
    {
        complex<double> w ((leftExtreme, upExtreme)); 
        double Scaled_imm = upExtreme - fraction*i;
        w.imag(Scaled_imm); //settiamo il valore di w su cui viene applicato NR

        for(int j = 0; j <= NCol; j++)
        { //ciclo sulle colonne

            double Scaled_real = leftExtreme + fraction*j;
            w.real(Scaled_real);
            printf("i = %d, j = %d) w attuale: R(%f) I(%f)\n", i, j, w.real(), w.imag());
            complex<double> c = w; //ricopiamo il numero di partenza 
            for(int q = 0; q < K; q++) //ciclo sulle iterazioni del metodo di NR
            {
                //metodo di Newton
                complex<double> f = Funzione(c);
                //printf("w: R(%f) I(%f) \nf: R(%f) I(%f)\n", w.real(), w.imag(), f.real(), f.imag()); 
                complex<double> df = DerivFunzione(c);
                //printf("df: R(%f) I(%f)\n", df.real(), df.imag()); 
                complex<double> d = (f)/(df);
                //printf("d: R(%f) I(%f)\n", d.real(), d.imag()); 
                complex<double> p; // nuovo valore, ovvero c al passo dopo
                p = c - d;
                //printf("Funzione(w)/DerivFunzione(w): R(%f) I(%f)\n", d.real(), d.imag()); 
                //printf("Allo step %d di NR, w ->  R(%f) I(%f)\n",it, w.real(), w.imag()); 
                
                complex<double> differenza = p - c;


                if(abs(differenza) < EPS)
                {   
                    fprintf(file_griglia, "%f %f %d\n", w.real(), w.imag(), q);
                    // si segna l'iterazione d'arresto rispetto al punto di partenza
                    break;
                }

                if(w.real() == -1.250000 && w.imag() == 0.250000)
                {
                    fprintf(file_curva, "%d %f %f\n", q, c.real(), c.imag());
                }

                    differenza.real(0);
                    differenza.imag(0);
                    c = p;
            } //fine ciclo sulle iterazioni

        }   //fine ciclo colonne
            
    } //fine ciclo sulle righe

return 0;

} // fine main

    

