#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <fstream>
#include <cmath>



using namespace std;



int main (int argc, char ** argv)
{
    // Salviamo la seconda soluzione dell'equazione prima in float, ovvero in singola precisione.

    float fPHI_2 = abs((-0.5) * (sqrtf( float(5) ) + 1.)); 
    double dPHI_2 = abs((-0.5) * (sqrt( double(5) ) + 1.)); 
    long double lPHI_2 = abs((-0.5) * (sqrtl( long(5) ) + 1.)); 

    float fPHI_1 = (0.5) * (sqrtf( float(5) ) - 1.); 
    double dPHI_1 = (0.5) * (sqrt( double(5) ) - 1.); 
    long double lPHI_1 = (0.5) * (sqrtl( long(5) ) - 1.); 

    float fchi_n_minus_2 = 1.;
    float fchi_n_minus_1 = fPHI_1;
    float fchi_n = 0.;

    double dchi_n_minus_2 = 1.;
    double dchi_n_minus_1 = dPHI_1;
    double dchi_n = 0.;

    long double lchi_n_minus_2 = 1.;
    long double lchi_n_minus_1 = lPHI_1;
    long double lchi_n = 0.;

    //  PER LA SECONDA PARTE

    long double eps = 1E-8 ; //è dell'ordine dei float, ma lavoro coi long double, pertanto mi aspetto di ottenere lo stesso riisultato per far capire che l'erorre è intrinseco
    long double due_chi_n_minus_2 = 1.;
    long double due_chi_n_minus_1 = lPHI_1 + eps*lPHI_2;
    long double due_chi_n = 0.;

    float FPHI2N = 0.;
    double DPHI2N = 0.;
    long double LPHI2N = 0.;


    // Stampiamo i valori dei PHI_2;

    printf("PHI_2 in FLOAT: %f \n", fPHI_2);
    printf("PHI_2 in DOUBLE: %lf \n", dPHI_2);
    printf("PHI_2 in LONG DOUBLE: %Lf \n", lPHI_2);

    // Stampiamo i valori dei PHI_1;

    printf("PHI_1 in FLOAT: %f \n", fPHI_1);
    printf("PHI_1 in DOUBLE: %lf \n", dPHI_1);
    printf("PHI_1 in LONG DOUBLE: %Lf \n", lPHI_1);

    
    // Definiamo il numero di passi che serviranno alla funzione di ricorsione

    int n = 0;
    int K = 100;
    int h = 1;

    // Salviamo i risultati su un file, tenendo conto del fatto che iteriamo su n

    FILE * DELTA_FLOAT = fopen("deltafloat.txt", "w");
    FILE * DELTA_DOUBLE = fopen("deltad.txt", "w");
    FILE * DELTA_LONG = fopen("deltalong.txt", "w");

    FILE * CHI_FLOAT = fopen("chif.txt", "w");
    FILE * CHI_DOUBLE = fopen("chid.txt", "w");
    FILE * CHI_LONG = fopen("chil.txt", "w");

    FILE * SECONDA_PARTE = fopen("secondaparte.txt", "w");

    //  Apriamo i file su cui stampare il valore di Phi2^N

    FILE * PHI2_FLOAT = fopen("FPHI2N.txt", "w");
    FILE * PHI2_DOUBLE = fopen("DPHI2N.txt", "w");
    FILE * PHI2_LONG = fopen("LPHI2N.txt", "w");

    cout << sizeof(long double) << endl;
    // Innestiamo un ciclo for per calcolare l'equazione a diversi n

    for( n = 2; n <= K; n = n + h ){
        // Definiamo la variabile chi_n nella quale depositeremo il risultato dell'equazione ricorsiva al passo n
        
        // Andiamo a calcolarci il valore di chi tenendo conto del fatto che la seconda potenza è affetta da errore di troncamento

        fchi_n = fchi_n_minus_2 - fchi_n_minus_1;
        dchi_n = dchi_n_minus_2 - dchi_n_minus_1;
        lchi_n = lchi_n_minus_2 - lchi_n_minus_1;
    
        // SECONDA PARTE
        
        due_chi_n = due_chi_n_minus_2 - due_chi_n_minus_1;
   
        // Ora abbiamo calcolato la "distanza" tra il valore calcolato con l'errore di troncamento e quello che assumiamo essere il valore vero

        float fdelta = abs(fchi_n - powf(fPHI_1, n));
        double ddelta = abs(dchi_n - pow(dPHI_1, n));
        long double ldelta =  abs(lchi_n - powl(lPHI_1, n));

       // NOTA IMPORTANTISSIMA: IN QUESTO CASO SO CHE LA DISCREPANZA VA COME EPS* PHI2, PASSANDO AL LOGARITMO HO CHE L'INTERCETTA è LOG(EPSILON), E CHE DEVO CALCOLARE CON UN FIT
       // log(delta) = log(eps) + n*log(abs(phi2)), dal fit posso trovare logeps e calcolarmi così la precisione.

        // Salviamo nei differenti file per plottare in MATLAB. I grafici su delta verranno plottati in scala logaritmica

        fprintf(DELTA_FLOAT, "%d %e\n", n, fdelta);
        fprintf(DELTA_DOUBLE, "%d %e\n", n, ddelta);
        fprintf(DELTA_LONG, "%d %e\n", n, (double) ldelta);

        // Vogliamo plottare anche i valori di chi^n.

        fprintf(CHI_FLOAT, "%d %e\n", n, fchi_n);
        fprintf(CHI_DOUBLE, "%d %e\n", n, dchi_n);
        fprintf(CHI_LONG, "%d %e\n", n, (double) lchi_n);
        fprintf(SECONDA_PARTE, "%d %e\n", n, (double) due_chi_n);

        // Poniamo le condizioni per ripetere il passo

        fchi_n_minus_2 = fchi_n_minus_1;
        fchi_n_minus_1 = fchi_n;

        dchi_n_minus_2 = dchi_n_minus_1;
        dchi_n_minus_1 = dchi_n;

        lchi_n_minus_2 = lchi_n_minus_1,
        lchi_n_minus_1 = lchi_n;

        due_chi_n_minus_2 = due_chi_n_minus_1;
        due_chi_n_minus_1 = due_chi_n;

    } //fine del ciclo for
    
    // NUOVO CICLO FOR CHE PARTE DA N = 0 E NON N = 2 PER PHI 2

    for( n = 1; n <= K; n = n + h ){
      // Andiamo a depositare i valori di PHI2^n

      FPHI2N = powf(fPHI_2, n);
      DPHI2N = pow(dPHI_2, n);
      LPHI2N = powl(lPHI_2, n);

      // Plottiamo anche PHI2^n

      fprintf(PHI2_FLOAT, "%d %e\n", n, FPHI2N);
      fprintf(PHI2_DOUBLE, "%d %e\n", n, DPHI2N);
      fprintf(PHI2_LONG, "%d %e\n", n, (double) LPHI2N);  //NOTA: IN QUESTO CASO NON VI SONO DIFFERENZE NOTEVOLI PERCHé L'ERRORE è DOMINANTE E LE CIFRE SIGNIFICATIVE NON VARIANO
    }
  return 0;
}



/* PUNTI CRITICI: iniziamo col plottare delta in scala logaritmica, dopodiché plottiamo chi per mostrare il momento in cui 
phi2 con errore diviene dominante sul valore di chi, e trovarlo in modo quantitativo (ad esempio, nel caso float mi aspetto che 
ciò avvenga a n = 16, ovvero quando chi vale 10^7, ovvero il numero massimo di mantissa di float), infatti epsilon è inversamente proporzionale
alla grandezza della mantissa */

