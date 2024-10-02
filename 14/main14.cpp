//IMPORTANTE RIFARE I GRAFICI CON TSTEPS 2000 + IN SDCALA LOGARITMICA, EULERO LINEARE, RK2 QUADR., RK4 QUADRATICA.

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>

#define N 6 // <- dimensioni del sistema di Eq. Diff.
#define T_start  0.  // <- inizio di evoluzione del sistema di Eq. Diff.
#define T_end   2.   // <- fine di evoluzione del sistema di Eq. Diff.
#define T_StepsMax 2000 

double y[N];

/* Function prototypes */
double f0(double t, double *v);
double f1(double t, double *v);
/* integratore e stampa*/
void integratore(double (*pf[N])(double t, double *v), int T_steps, char method[]);
void stampa(int t, double *v, FILE * output);
/* passi elementari */
void Eulero(double t, double *v, double h, double (*pf[N])(double t, double *v));
void RK2(double t, double *v, double h, double (*pf[N])(double t, double *v));
void RK4(double t, double *v, double h, double (*pf[N])(double t, double *v));

FILE * EULEROFILE = fopen("Eulero.txt", "w"); 
FILE * RK2FILE = fopen("RK2.txt", "w"); 
FILE * RK4FILE = fopen("RK4.txt", "w"); 

int main() {

  double (*f[N])(double t, double *v)={f0,f1}; // <- definizione del sistema di Eq. Diff.
  int T_steps; //numero di passi in 2 unità temporali

  y[0]=0.; // THETA 0
  y[1]=1.; // VELOCITA' 0

  

//LE condizioni iniziali le reinizializzeremo tra un metodo e l'altro
  for(int i = 20; i <= T_StepsMax; i = i + 20){
   T_steps=i;  // <- numero di passi di integrazione
   y[0]=0.; // THETA 0
   y[1]=1.; // VELOCITA' 0
   integratore(f,T_steps,"Eulero");
   stampa(T_steps, y, EULEROFILE);
   y[0]=0.; // THETA 0
   y[1]=1.; // VELOCITA' 0
   integratore(f,T_steps,"RK2");
   stampa(T_steps, y, RK2FILE);
   y[0]=0.; // THETA 0
   y[1]=1.; // VELOCITA' 0
   integratore(f,T_steps,"RK4");
   stampa(T_steps, y, RK4FILE);
  }
  return 0;
 }


/* Definizione delle funzioni del sistema di eq. diff. */
double f0(double t, double *v)
{
  return v[1]; // ftheta = dtheta(t)/dt
}

double f1(double t, double *v)
{
  return -v[0]; // fz = -theta(t)
}

/* ***************************************************************** */
/* NB: la parte qui sotto non cambia con il sistema di eq. diff. !!! */
/* ***************************************************************** */

/* Definizione dell'integratore */

void integratore(double (*pf[N])(double t, double *v), int T_steps, char method[])
{
  void (*step)(double t, double *v, double h, double (*pf[N])(double t, double *v));
  double t=(double)T_start;
  double h = ((double)T_end - (double)T_start ) / (double)T_steps ; // definiamo l'incremento

  
  /* *********** scelta del Metodo *********** */  
  if(strcmp(method, "Eulero") == 0) {
    step = Eulero;
  } else if (strcmp(method, "RK2") == 0) {
    step = RK2;
  }  else if (strcmp(method, "RK4") == 0) {
    step = RK4;
  }  else printf("metodo %s non esiste!\n",method);

  /* *********** evoluzione *********** */  
  for(int k=0; k<T_steps; k++) 
    {
      step(t, y, h, pf);
      t = t + h;}
  return;
}


/* Funzione di stampa */
void stampa(int t, double *v, FILE * output)
{
  fprintf(output, "%d ",t);
  for(int n=0;n<N;n++) {
    if(n == 1) fprintf(output, "%.16lf", fabs(v[n]));
    else {fprintf(output, "%1.16lf ", v[n]);}
    
    };
  fprintf(output, "\n");
}
  


/* Definizione del passo elementare di integrazione */

void Eulero(double t, double *v, double h, double (*pf[N])(double t, double *v))
{
  double k1[N];

  for(int n=0;n<N;n++) k1[n] = h * pf[n](t,v);

  for(int n=0;n<N;n++) v[n] = v[n] + k1[n];
  
  return;
}
  

void RK2(double t, double *v, double h, double (*pf[N])(double t, double *v))
{
  double k1[N], k2[N], v_tmp[N];

  for(int n=0;n<N;n++) k1[n] = h * pf[n](t,v);

  for(int n=0;n<N;n++) v_tmp[n] = v[n] + k1[n] / 2.0;
  for(int n=0;n<N;n++) k2[n] = h * pf[n](t+h/2.,v_tmp);

  for(int n=0;n<N;n++) v[n] = v[n] + k2[n];
  
  return;
}
  

void RK4(double t, double *v, double h, double (*pf[N])(double t, double *v))
{
  double k1[N], k2[N], k3[N], k4[N], v_tmp[N];

  for(int n=0;n<N;n++) k1[n] = h * pf[n](t,v);

  for(int n=0;n<N;n++) v_tmp[n] = v[n] + k1[n] / 2.0;
  for(int n=0;n<N;n++) k2[n] = h * pf[n](t+h/2.,v_tmp);

  for(int n=0;n<N;n++) v_tmp[n] = v[n] + k2[n] / 2.0;
  for(int n=0;n<N;n++) k3[n] = h * pf[n](t+h/2.,v_tmp);
 
  for(int n=0;n<N;n++) v_tmp[n] = v[n] + k3[n];
  for(int n=0;n<N;n++) k4[n] = h * pf[n](t+h,v_tmp);

  for(int n=0;n<N;n++) v[n] = v[n] + (k1[n] + 2. * k2[n] + 2. * k3[n] + k4[n])/6.;
  
  return;
}