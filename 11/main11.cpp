//Generare esattamente numeri random con le seguenti distribuzion

#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>          
#include <cmath>
#include <vector>
#define FUNC(x) ((*func)(x))
#define nMAX 1e5

using namespace std;

double rand_range (double min, double max)
{
  return min + (max - min) * rand () / static_cast<double> (RAND_MAX);
}


double func1 (double x)
{
    double y;
    double coeff = (1-exp(-2));
    y = -log(1-x*coeff);
    return y;
}

double func2 (double x)
{
    double y;
    double coeff = -1.;
    y = 1 - log(-1*x + 1);
    return y;
}

double func3 (double x)
{
    double y;
    double coeff = -1.;
    y = sqrt(-log(coeff*x +1));
    return y;
}

int main ()
{   

    // numeri random con exp(-x) tra 0 e 2
    FILE * exp_02 = fopen("exp02.txt", "w");    
    double x1 = 0.;
    for(int i = 1; i <= nMAX + 1; i = i + 50){
        x1 = func1(rand_range(0.,1.)); //saranno distribuiti come un esponenziale tra 0 e 2
        fprintf(exp_02, "%.16lf\n", x1);
    }

    // numeri random con exp(-x) tra 1 e +inf
    FILE * exp_1inf = fopen("exp1inf.txt", "w");    
    double x2 = 0.;
    for(int i = 1; i <= nMAX + 1; i = i + 50){
        x2 = func2(rand_range(0.,1.)); //saranno distribuiti come un esponenziale tra 0 e 2
        fprintf(exp_1inf, "%.16lf\n", x2);
    }

    // numeri random con x*exp(-x^2) tra 0 e +inf
    FILE * exp_0inf = fopen("exp0inf.txt", "w");    
    double x3 = 0.;
    for(int i = 1; i <= nMAX + 1; i = i + 50){
        x3 = func3(rand_range(0.,1.)); //saranno distribuiti come un esponenziale tra 0 e 2
        fprintf(exp_0inf, "%.16lf\n", x3);
    }
    
    return 0;
}