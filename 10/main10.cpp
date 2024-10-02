#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#define FUNC(x) ((*func)(x))
#define nMAX 1e6    


double rand_range (double min, double max)
{
  return min + (max - min) * rand () / static_cast<double> (RAND_MAX);
}



double hit_or_miss(double (*func)(double), double a, double b, double ymax, int N){
    double xrandNum, yrandNum;
    double counter_in = 0.;
    double val;

    for(int i = 0; i < N; i++){
        xrandNum = rand_range(a,b); // trovare a mano il massimo, non in modo algoritmico per il rettangolo.
        yrandNum = rand_range(0, ymax);

        if(yrandNum < func(xrandNum))
        {
            counter_in++;
        }
    }
    
    val = (fabs(b-a)*ymax*counter_in)/N; // R* Nn/N
    
    return val;
}

double montecarlo_sampling(double (*func) (double), double a, double b, int N) {
    double sum = 0.;
    double xr = 0.;
    for (int i = 0; i < N; i++){
        xr = rand_range(a,b);
        sum += func(xr);
    }
    double val = 0.;
    val = (fabs(b-a)*sum)/(N);
    return val;
}


double func3(double x)
{
    double y = 0.0;
    y = pow(x,7) * exp(-x);
    return y;
}

double func4(double x){
    double y = 0.0;
    y = cosh(x);
    return y;
}

double func5(double x){
    double y = 0.0;
    y = pow(x,2) + x*sin(4*x);
    return y;
}

using namespace std;

int main ()
{
    double es3_hit;
    double es3_s;
    double a3 = 0.;
    double b3 = 5.;
    double y3 = 526.4;
    double TrueRes3 = 5040 - (648240/exp(5)); //672.19

    double es4_hit;
    double es4_s;
    double a4 = 3.;
    double b4 = 8.;
    double y4 = 1490.5 ;
    double TrueRes4 = sinh(8) - sinh(3); // 1480.460

    double es5_hit;
    double es5_s;
    double a5 = -1.;
    double b5 = 8.;
    double y5 = 68.45;
    double TrueRes5 = 1./16. * (2736 + sin(4) + sin(32) - 4*cos(4) - 32*cos(32)); // 169.48

    FILE * DEVHM_3 = fopen("DEVHM_3.txt", "w");
    FILE * DEVHM_4 = fopen("DEVHM_4.txt", "w");
    FILE * DEVHM_5 = fopen("DEVHM_5.txt", "w");
    
    FILE * DEVS_3 = fopen("DEVS_3.txt", "w");
    FILE * DEVS_4 = fopen("DEVS_4.txt", "w");
    FILE * DEVS_5 = fopen("DEVS_5.txt", "w");

    FILE * F_N = fopen("N.txt", "w");

    

    for(int N = 1; N <= nMAX + 1; N = N + 10000){

    es3_hit = hit_or_miss(func3, a3, b3, y3, N);
    fprintf(DEVHM_3, "%.16lf\n",fabs(TrueRes3-es3_hit));
    //printf("Esercizio 3 con Hit or Miss: %lf\nLa differenza col valore vero: %lf\n", es3_hit, fabs(TrueRes3-es3_hit)); 

    es4_hit = hit_or_miss(func4, a4, b4, y4, N);
    fprintf(DEVHM_4, "%.16lf\n",fabs(TrueRes4-es4_hit));
    //printf("Esercizio 4 con Hit or Miss: %lf\nLa differenza col valore vero: %lf\n", es4_hit, fabs(TrueRes4-es4_hit)); 

    es5_hit = hit_or_miss(func5, a5, b5, y5, N);
    fprintf(DEVHM_5, "%.16lf\n",fabs(TrueRes5-es5_hit));
    //printf("Esercizio 5 con Hit or Miss: %lf\nLa differenza col valore vero: %lf\n", es5_hit, fabs(TrueRes5-es5_hit)); 

    es3_s = montecarlo_sampling(func3, a3, b3, N);
    fprintf(DEVS_3, "%.16lf\n",fabs(TrueRes3-es3_s));

    es4_s = montecarlo_sampling(func4, a4, b4, N);
    fprintf(DEVS_4, "%.16lf\n", fabs(TrueRes4-es4_s));

    es5_s = montecarlo_sampling(func5, a5, b5, N);
    fprintf(DEVS_5, "%.16lf\n", fabs(TrueRes5-es5_s));

    fprintf(F_N, "%d\n", N);

    }


    
    return 0;
}