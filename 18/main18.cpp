// FARE I PLOT

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>

#define N 18 // <- dimensioni del sistema di Eq. Diff.
#define T_start  0.  // <- inizio di evoluzione del sistema di Eq. Diff.
#define T_end   40.   // <- fine di evoluzione del sistema di Eq. Diff.

double y[N]; //condizioni iniziali

#define m 0.3
#define m1 1.6
#define m2 0.4
#define m3 0.4

/* Function prototypes */
double f0(double t, double *v); //vx1
double f1(double t, double *v); //x1
double f2(double t, double *v); //vy1
double f3(double t, double *v); //y1
double f4(double t, double *v); //vz1
double f5(double t, double *v); //z1

double f6(double t, double *v); //vx2
double f7(double t, double *v); //x2
double f8(double t, double *v); //vy2
double f9(double t, double *v); //y2
double f10(double t, double *v); //vz2
double f11(double t, double *v); //z2

double f12(double t, double *v); //vx3
double f13(double t, double *v); //x3
double f14(double t, double *v); //vy3
double f15(double t, double *v); //y3
double f16(double t, double *v); //vz3
double f17(double t, double *v); //z3

/* integratore e stampa*/
void integratore(double (*pf[N])(double t, double *v), int T_steps, char method[]);
void stampa(double t, double *v, FILE *output);
/* passi elementari */
void RK4(double t, double *v, double h, double (*pf[N])(double t, double *v));

double distanza(double *v, int i, int j)
{
  // Inserisci i,j come la posizione della prima coordinata del vettore
  double diffx = v[i] - v[j];
  double diffy = v[i + 2] - v[j+2];
  double diffz = v[i + 4] - v[j+4];
  double distanza = sqrt(pow(diffx,2) + pow(diffy,2) + pow(diffz,2));
  return distanza;
}

double modulo(double *v, int i)
{
  double mod = sqrt(pow(v[i],2) + pow(v[i+2],2) + pow(v[i+4], 2));
  return mod;
}


//DEFINIAMO LE MASSE
double M1 = m1;
double M2 = m2;
double M3 = m3;

int main() {

  double (*f[N])(double t, double *v)={f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,
                                      f10,f11,f12,f13,f14,f15,f16,f17}; // <- definizione del sistema di Eq. Diff.
  int T_steps;

// y vettore delle condizioni iniziali, da inizializzare

  y[0] = 0; //vx1
  y[1] = 1; //x1
  y[2] = 0.4; //vy1
  y[3] = 0;  //y1
  y[4] = 0;  //vz1
  y[5] = 0;  //z1
  y[6] = 0; //vx2
  y[7] = -1; //x2
  y[8] = -0.8; //vy2 
  y[9] = 0;  //y2
  y[10] = 0.7;  //vz2
  y[11] = 0;  //z2
  y[12] = 0;  //vx3
  y[13] = 0;  //x3
  y[14] = -0.8;  //vy3
  y[15] = 0; //y3
  y[16] = -0.7; //vz3
  y[17] = 0; //z3
  
  T_steps=10000;  // <- numero di passi di integrazione
  integratore(f,T_steps,"RK4");
   
  return 0;
 }



/* Definizione delle funzioni del sistema di eq. diff. */

//pari = velocità
//dispari = posizione

double f0(double t, double *v) //vx1
{
  double diff12x = v[1] - v[7]; //entrano [i] dispari poiché posizioni
  double diff13x = v[1] - v[13]; 
      return -1*(M2*((diff12x)/(pow(fabs(distanza(v, 1, 7)),3))) +
                 M3*((diff13x)/(pow(fabs(distanza(v, 1, 13)),3))) 
                 ) ;
}

double f1(double t, double *v) 
{
  return v[0];
}

double f2(double t, double *v)
{
  double diff12y = v[3] - v[9];
  double diff13y = v[3] - v[15]; 
  return -1*(M2*((diff12y)/(pow(abs(distanza(v, 1, 7)),3))) +
                 M3*((diff13y)/(pow(abs(distanza(v, 1, 13)),3)))
                 ) ;
}

double f3(double t, double *v)
{
  return v[2];
}

double f4(double t, double *v)
{
  double diff12z = v[5] - v[11];
  double diff13z = v[5] - v[17]; 
  return -1*(M2*((diff12z)/(pow(fabs(distanza(v, 1, 7)),3))) +
                 M3*((diff13z)/(pow(fabs(distanza(v, 1, 7)),3)))
                 ) ;
}

double f5(double t, double *v)
{
  return v[4];
}

double f6(double t, double *v)
{

  double diff21x = v[7] - v[1];
  double diff23x = v[7] - v[13];
  return -1*(M1*((diff21x)/(pow(distanza(v, 7, 1),3))) +
                 M3*((diff23x)/(pow(distanza(v, 7, 13),3)))
                 ) ;
}

double f7(double t, double *v)
{
  return v[6];
}

double f8(double t, double *v)
{

  double diff21y = v[9] - v[3];
  double diff23y = v[9] - v[15];
  return -1*(M1*((diff21y)/(pow(distanza(v, 7, 1),3))) +
                 M3*((diff23y)/(pow(distanza(v, 7, 13),3)))
                 ) ;
}

double f9(double t, double *v)
{
  return v[8];
}

double f10(double t, double *v)
{

  double diff21z = v[11] - v[5];
  double diff23z = v[11] - v[17];
  return -1*(M1*((diff21z)/(pow(distanza(v, 7, 1),3))) +
                 M3*((diff23z)/(pow(distanza(v, 7, 13),3)))
                 ) ;
}

double f11(double t, double *v)
{
  return v[10];
}

double f12(double t, double *v)
{

  double diff31x = v[13] - v[1]; //31
  double diff32x = v[13] - v[7]; //32
  return -1*(M1*((diff31x)/(pow(distanza(v, 13, 1),3))) +
                 M2*((diff32x)/(pow(distanza(v, 13, 7),3)))
                 ) ;
}

double f13(double t, double *v)
{
  return v[12];
}

double f14(double t, double *v)
{

  double diff31y = v[15] - v[3];
  double diff32y = v[15] - v[9];
  return -1*(M1*((diff31y)/(pow(distanza(v, 13, 1),3))) +
                 M2*((diff32y)/(pow(distanza(v, 13, 7),3)))
                 ) ;
}

double f15(double t, double *v)
{
  return v[14];
}

double f16(double t, double *v)
{

  double diff31z = v[17] - v[5];
  double diff32z = v[17] - v[11];
  return -1*(m1*((diff31z)/(pow(distanza(v, 13, 1),3))) +
                 m2*((diff32z)/(pow(distanza(v, 13, 7),3)))
                 ) ;
}

double f17(double t, double *v)
{
  return v[16];
}




/* Definizione dell'integratore */

void integratore(double (*pf[N])(double t, double *v), int T_steps, char method[])
{
  void (*step)(double t, double *v, double h, double (*pf[N])(double t, double *v));
  double t=(double)T_start;
  double h = ((double)T_end - (double)T_start ) / (double)T_steps ;

  FILE *output = fopen(method, "w"); 
  FILE *potenziale = fopen("Potenziale.txt", "w"); 
  FILE *energiaFILE = fopen("Energia.txt", "w"); 
  
  /* *********** scelta del Metodo *********** */  
  if (strcmp(method, "RK4") == 0) {
    step = RK4;
  }  else printf("metodo %s non esiste!\n",method);

  /* *********** evoluzione *********** */  
  for(int k=0; k<T_steps; k++)
    {
      step(t, y, h, pf);
      t = t + h;
      stampa(t,y,output);

      double energy;
      double potential;
      double kinetic;
      //scriviamo l'energia potenziale, entreranno [i] dispari
      potential = -1.*( (M1*M2)/(distanza(y, 1, 7)) +  (M1*M3)/(distanza(y, 1, 13)) +  (M2*M3)/(distanza(y, 7, 13)) );
      //scriviamo l'energia cinetica, entreranno [i] pari
      kinetic = 0.5*(M1*pow(modulo(y, 0),2) + M2*pow(modulo(y,6),2) + M2*pow(modulo(y, 12),2)); //nel modulo entra la prima componente del vettore velocità, le altre vengono ricostruite nella funzione.
      energy = potential + kinetic;
      fprintf(potenziale, "%f %f\n", t, potential);
      fprintf(energiaFILE, "%f %f\n", t, energy);
    }
    fclose(output);
    fclose(energiaFILE);
  return;
}


/* Funzione di stampa */
void stampa(double t, double *v, FILE *output)
{
  fprintf(output, "%1.16lf ",t);
  for(int n=0;n<N;n++) fprintf(output, "%1.16lf ", v[n]);
  fprintf(output, "\n");
}
  


/* Definizione del passo elementare di integrazione */
  
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