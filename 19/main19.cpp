#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#define EPS 1e-12  // precisione max
using namespace std;

double func(double x)
{
    return 2*pow(x,2) - 3*x + 1;
}

double derivFunc(double x)
{
    return 4*x - 3;
}

double Legendre(double x)
{
    return (46189*pow(x,10) - 109395*pow(x,8) + 90090*pow(x,6) -30030*pow(x,4) + 3465*pow(x,2) - 63)/(256);
}

double derivLegendre(double x)
{
    return (46189*10*pow(x,9) - 109395*8*pow(x,7) + 90090*6*pow(x,5) -30030*4*pow(x,3) + 3465*2*x)/(256);
}


// metodo bracketing + bisezione
double bracketing_seek(double a, double b, double (*func)(double))
{
    if(func(a) * func(b) > 0){
    printf("BRA: no roots\n");
    return 0;}
    
    double c = a;    

    while((b-a) >= EPS)
    {
        c = (a+b)/2.;
        if (func(c) == 0.) break;

        else if(func(c)*func(a)<0) b=c;
        else a=c;
    }

    //printf("Lo zero trovato (Bracketing)  %f\n", c);
    return c;
}


// metodo Newton-Raphson
double newton_raphson(double x, double (*func)(double), double (*derivFunc)(double))
{
    double delta = func(x) / derivFunc(x);
    while(abs(delta) >= EPS)
    {
        delta = func(x)/derivFunc(x);

        x = x - delta;
    }

   // printf("Lo zero trovato (NR) %f", x);
    return x;

}


using namespace std;

int main ()
{
    char ANSWER;
    char FUNCTION;
    double estremosx;
    double estremodx;
    double a;
    double b;

    //FILE * estremi = fopen("estremiLeg.txt", "w"); 

    printf("Vuoi operare con f(x) o con il polinomio di Legendre? F/L ");
    cin >> FUNCTION;
    //  PER IL POLINOMIO DI LEGENDRE
    if(FUNCTION == 'L')
    {
        printf("Vuoi cercare gli intervalli con zero in un range? Y/N ");
        cin >> ANSWER;
        if(ANSWER == 'Y')
        {
            printf("Scrivi l'estremo sinistro e l'estremo destro ");
            cin >> estremosx;
            cin >> estremodx;
            double x;
            int j = 1;
            double M = 100000; // DEVE ANDARE A VEDERE NELL'INTERVALLO DISCRETIZZATO DOVE LA FUNZIONE CAMBIA SEGNO
            double fraction = abs(estremosx - estremodx)/(M);

            printf("La frazione: %.4lf\n", fraction);

            for(int i = 1; i <= M; i++){
                if(Legendre(estremosx)*Legendre(estremosx + fraction) <= 0){
                printf("LEG %d) Uno zero in [%.16f, %.16f]\n", j, estremosx, (estremosx + fraction));
                //fprintf(estremi, "%d %.16lf %.16lf\n", j, estremosx, (estremosx + fraction));
                j++;
                }
        
            estremosx = estremosx + fraction;
            
            }
        }
        else if(ANSWER == 'N')
        {
            printf("Scrivi prima il valore di a e poi quello di b: ");
            cin >> a;
            cin >> b;
            double xB;
            double xNR;
            xB = bracketing_seek(a, b, Legendre);
            printf("LEG) Il metodo Bracketing ha trovato x = %.16f \n", xB);
            xNR = newton_raphson(a, Legendre, derivLegendre);
            printf("LEG) Il metodo Newton Raphson ha trovato x = %.16f \n", xNR);
        }

        else printf("Errore!"); return 0;
    }

    // PER IL POLINOMIO DI F(X)
    if(FUNCTION == 'F')
    {
        printf("Vuoi cercare gli intervalli con zero in un range? Y/N ");
        cin >> ANSWER;
        if(ANSWER == 'Y')
        {
            printf("Scrivi l'estremo sinistro e poi l'estremo destro: ");
            cin >> estremosx;
            cin >> estremodx;
            double x;
            double M = 1000; // DEVE ANDARE A VEDERE NELL'INTERVALLO DISCRETIZZATO DOVE LA FUNZIONE CAMBIA SEGNO
            double fraction = abs(estremosx - estremodx)/(M);
            for(int i = 1; i <= M; i++){
                if(func(estremosx)*func(estremosx + fraction) < 0){
                printf("f) Uno zero in [%.16f, %.16f]\n", estremosx, (estremosx + fraction));
                }
            estremosx = estremosx + fraction;
            }
        }
        else if(ANSWER == 'N')
        {
            printf("Scrivi prima il valore di a e poi quello di b: ");
            cin >> a;
            cin >> b;
            double xB;
            double xNR;
            xB = bracketing_seek(a, b, func);
            printf("f) Il metodo Bracketing ha trovato x = %.16f \n", xB);
            xNR = newton_raphson(a, func, derivFunc);
            printf("f) Il metodo Newton Raphson ha trovato x = %.16f \n", xNR);
        }

        else printf("Errore!"); return 0;
    }


    else printf("Errore!"); return 0;
    
    
    return 0;
}

