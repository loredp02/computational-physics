/*
Calcolare con Trapezio, Simson e Romberg l'integrale tra 5 e 0 di x^7 e^-x, per TRAP e SIMP plottare le deviazioni in funzione di 1/N con N numero 
di punti. Poi calcolare con GAUSS tramite Legendre (2,4,8 punti) e Laguerre (2,4,8, con alfa = 0)
*/
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#define FUNC(x) ((*func)(x))


using namespace std;

// Metodo del trapezio 

double trapzd (double (*func)(double), double a, double b, int n){
    double x, tnm, sum, del;
    static double s;
    int it,j;

    if(n==1){ return (s = 0.5*(b-a)*(FUNC(a) + FUNC(b)));
    } else {
        for (it=1, j=1; j < n-1; j++) it <<= 1;
        tnm = it;
        del = (b-a)/tnm;
        x = a + (double) (1./2.)*del;
        for (sum=0.0, j=1; j<=it;j++,x+=del) sum += FUNC(x);
        s= (double) (1./2.) *(s+(b-a)*sum/tnm);  
        return s;
    }
}

// Metodo di Simpson

double qsimp(double (*func)(double), double a, double b, int n){
    int j;
    float s,st,ost=0.0,os=0.0;

    for (j=1;j<=n;j++) {
        st=trapzd(func,a,b,j);
        s=(4.0*st-ost)/3.0;
        os=s;
        ost=st;
    }

    return s;
}

// METODO DI ROMBERG

double romberg(double (*func)(double), double a, double b, int N) {
    double h[N+1], r[N+1][N+1];
    for (int i = 1; i < N + 1; ++i) {
        h[i] = (b - a) / pow(2, i - 1);
    }
    r[1][1] = h[1] / 2 * (func(a) + func(b));
    for (int i = 2; i < N + 1; ++i) {
        double coeff = 0;
        for (int k = 1; k <= pow(2, i - 2); ++k) {
            coeff += func(a + (2 * k - 1) * h[i]);
        }
        r[i][1] = 0.5 * (r[i - 1][1] + h[i - 1] * coeff);
    }
    
    for (int i = 2; i < N + 1; ++i) {
        for (int j = 2; j <= i; ++j) {
            r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
        }
    }
    return r[N][N];
}


double gauss_legendre(double (*func) (double), double a, double b, int n){
    // Cambiamo l'intervallo di integrazione della funzione tra [-1,1]
    vector<double> vec_x2 {-5.77350269189625764507e-01, 5.77350269189625764507e-01};
    vector<double> vec_x4 {-8.61136311594052575248e-01, -3.39981043584856264792e-01, 3.39981043584856264792e-01, 8.61136311594052575248e-01};
    vector<double> vec_x8 {-9.60289856497536231661e-01, -7.96666477413626739567e-01, -5.25532409916328985830e-01, -1.83434642495649804936e-01,
                            1.83434642495649804936e-01,  5.25532409916328985830e-01, 7.96666477413626739567e-01, 9.60289856497536231661e-01};

    vector<double> vec_w2 {1, 1};
    vector<double> vec_w4 {3.47854845137453857383e-01, 6.52145154862546142644e-01, 6.52145154862546142644e-01, 3.47854845137453857383e-01};
    vector<double> vec_w8 {1.01228536290376259154e-01, 2.22381034453374470546e-01,  3.13706645877887287338e-01,  3.62683783378361982976e-01,
                           3.62683783378361982976e-01, 3.13706645877887287338e-01,  2.22381034453374470546e-01, 1.01228536290376259154e-01};

    double res = 0.0;
    double coeff = (b-a)*(1./2.);
    double sum = 0.0;

    if(n==2){
        for(int i = 0; i < n; i++){
            sum += vec_w2.at(i)*func(coeff*vec_x2.at(i) + (a+b)*(1./2.));
        }
    }

    else if(n==4){
        for(int i = 0; i < n; i++){
        sum += vec_w4.at(i)*func(coeff*vec_x4.at(i) + (a+b)*(1./2.));
        }
    }

    else if(n==8){
        for(int i = 0; i < n; i++){
        sum += vec_w8.at(i)*func(coeff*vec_x8.at(i) + (a+b)*(1./2.));
        }
    }

    res = coeff*sum;

    return res;
}

double gauss_laguerre(double (*func) (double), int n){

    vector<double> vec_x2 {5.85786437626904951182e-01, 3.41421356237309504876e+00};
    vector<double> vec_x4 {3.22547689619392311802e-01, 1.74576110115834657569e+00, 4.53662029692112798345e+00, 9.39507091230113312950e+00};
    vector<double> vec_x8 {1.70279632305100999786e-01, 9.03701776799379912170e-01, 2.25108662986613068929e+00, 4.26670017028765879378e+00,
                           7.04590540239346569719e+00, 1.07585160101809952241e+01, 1.57406786412780045781e+01, 2.28631317368892641052e+01};

    vector<double> vec_w2 {8.53553390593273762191e-01, 1.46446609406726237796e-01};
    vector<double> vec_w4 {6.03154104341633601660e-01, 3.57418692437799686640e-01, 3.88879085150053842740e-02, 5.39294705561327450102e-04};
    vector<double> vec_w8 {3.69188589341637529929e-01, 4.18786780814342956078e-01, 1.75794986637171805706e-01, 3.33434922612156515224e-02,
                           2.79453623522567252491e-03, 9.07650877335821310457e-05, 8.48574671627253154502e-07, 1.04800117487151038157e-09}; 
                        
    double res = 0.0;
    double sum = 0.0;

    if(n==2){
        for(int i = 0; i < n; i++){
            sum += vec_w2.at(i)*func(vec_x2.at(i));
        }
    }

    else if(n==4){
        for(int i = 0; i < n; i++){
        sum += vec_w4.at(i)*func(vec_x4.at(i));
        }
    }

    else if(n==8){
        for(int i = 0; i < n; i++){
        sum += vec_w8.at(i)*func(vec_x8.at(i));
        }
    }

    return sum;

};






///////////////////////////////////////////////////////////////////////////////////////////////////7

double func(double x){
    double y = 0.0;
    y = pow(x,7) * exp(-x);
    return y;
}

double func_lag1(double x){
    double y = 0.0;
    y = pow(x,7);
    return y;
}

double func_lagshift(double x){
    double y  = 0.0;
    y = pow(x+5, 7);
    return y;
}

int main ()
{   
    double a = 0.;
    double b = 5.;
    double TrueRes = 5040 - (648240/exp(5));
    printf("Il risultato vero: %.16lf \n", TrueRes);

    int nMAX = 20;
    int n = 0;


    FILE * TRAPDEV = fopen("trapezdev.txt", "w");
    FILE * SIMPDEV = fopen("simpdev.txt", "w");
    FILE * ROMBDEV = fopen("rombdev.txt", "w");


    for(int n = 1; n <= nMAX; n += 1){
        ///////////////////////////////////////////////////////////////////////////////////////// TRAPEZIO
        double TrapezioRis = 0.0;
        TrapezioRis = trapzd(func, a, b, n);
        double TrapezioDev = fabs(TrueRes - TrapezioRis);
        fprintf(TRAPDEV, "%f %.16lf\n", 1./n, TrapezioDev);

        printf("Il risultato col trapezio -> %.16lf\n", TrapezioRis);
        ///////////////////////////////////////////////////////////////////////////////////////// SIMPSON
        double SimpsonRis = 0.0;
        SimpsonRis = qsimp(func, a, b, n);
        double SimpsonDev = fabs(TrueRes - SimpsonRis);
        fprintf(SIMPDEV, "%f %.16lf\n", 1./n, SimpsonDev);
        
        printf("Il risultato con Simpson -> %.16lf\n", SimpsonRis);
        ///////////////////////////////////////////////////////////////////////////////////////// ROMBERG
        double RombergRis = 0.0;
        RombergRis = romberg(func, a, b, n);
        double RombDev = fabs(TrueRes - RombergRis);
        fprintf(ROMBDEV, "%f %.16lf\n", 1./n, RombDev);
        printf("Il risultato con Romberg -> %.16lf\n\n", RombergRis);

    }

   double gauss_leg2 = 0.0;
   double gauss_leg4 = 0.0;
   double gauss_leg8 = 0.0;

   gauss_leg2 = gauss_legendre(func, a, b, 2);
   gauss_leg4 = gauss_legendre(func, a, b, 4);
   gauss_leg8 = gauss_legendre(func, a, b, 8);

   printf("Gauss Legendre con 2 punti: %.16lf \n", gauss_leg2);
   printf("Gauss Legendre con 4 punti: %.16lf \n", gauss_leg4);
   printf("Gauss Legendre con 8 punti: %.16lf \n", gauss_leg8);

   printf("------------------------------------------- \n");    

   double gauss_lag2 = 0.0;
   double gauss_lag4 = 0.0;
   double gauss_lag8 = 0.0;
   double diff_coeff = exp(-5);

   gauss_lag2 = gauss_laguerre(func_lag1, 2) - diff_coeff * gauss_laguerre(func_lagshift, 2);
   gauss_lag4 = gauss_laguerre(func_lag1, 4) - diff_coeff * gauss_laguerre(func_lagshift, 4);
   gauss_lag8 = gauss_laguerre(func_lag1, 8) - diff_coeff * gauss_laguerre(func_lagshift, 8);

   printf("Gauss Laguerre con 2 punti: %.16lf \n", gauss_lag2);
   printf("Gauss Laguerre con 4 punti: %.16lf \n", gauss_lag4);
   printf("Gauss Laguerre con 8 punti: %.16lf \n", gauss_lag8);

return 0;
}


