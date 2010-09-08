#include<stdio.h>
#include<stdlib.h> // has exit() function
#include<string.h>
#include<math.h>
#include<Rmath.h>


typedef struct {
  double **X, **cumX;
  double *Hvec, *av_Ni;
  double Hbar, muNi, sdNi;
} XMAT;



XMAT* Xdip(double **Xmat, int nrow, int ncol);
XMAT* Xhap(double **Xmat, int nrow, int ncol);


double NRroundit(double d, int dig);

float ran1(long *idum);
float ran2(long *idum);
void NRsort(int nr, double *arr);
