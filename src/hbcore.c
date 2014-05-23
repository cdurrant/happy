#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>

#include "hbcore.h"



XMAT* Xdip(double **Xmat, int nrow, int ncol) {

  int i,j;
  double new, dcheck, dncol, dnrow, maxH, Hbar;
  double mu_avNi, sd_avNi, var_avNi;

  double *Hvec=NULL, *av_Ni=NULL;
  double **newX=NULL, **cumX=NULL;

  XMAT *full=NULL;


  // memory allocation

  full = (XMAT*)calloc(1,sizeof(XMAT));

  newX = (double**)calloc(nrow,sizeof(double*));
  cumX = (double**)calloc(nrow,sizeof(double*));
  for (i=0; i<nrow; i++) {
    newX[i] = (double*)calloc(ncol,sizeof(double));
    cumX[i] = (double*)calloc(ncol,sizeof(double));
  }


  // convert founder strain expectations to probabilities and round

  for (j=0; j<nrow; j++) {
    for (i=0; i<ncol; i++) {
      new = Xmat[j][i]/2.0;
      new = NRroundit(new,4);
      newX[j][i] = new;
    }
  }

  // normalise probability vectors and calculate entropy, etc

  Hvec = (double*)calloc(nrow,sizeof(double));
  av_Ni = (double*)calloc(ncol,sizeof(double));

  dncol = (double)ncol;
  ncol = (int)ncol;

  dnrow = (double)nrow;
  nrow = (int)nrow;

  maxH = 0.0;
  for (j=0; j<ncol; j++) {
    maxH = maxH - (1.0/dncol)*(log10((1.0/dncol))/log10(2));
  }

  Hbar = 0.0;
  for (i=0; i<nrow; i++) {

    dcheck = 0.0;
    for (j=0; j<ncol; j++) {
      dcheck = dcheck + newX[i][j];
    }

    for (j=0; j<ncol; j++) {
      newX[i][j] = newX[i][j]/dcheck;
    }

    dcheck = 0.0;
    for (j=0; j<ncol; j++) {
      dcheck = dcheck + newX[i][j];
      cumX[i][j] = dcheck;
    }
    if ( (dcheck <=0.99999) || (dcheck >= 1.00001) ) {
      Rprintf("individual %i : dcheck = %e ERROR HMM probs do not sum to 1\n",i,dcheck);
    }

    Hvec[i] = 0.0;
    for (j=0; j<ncol; j++) {
      if ( newX[i][j] != 0.0 ) {
        Hvec[i] = Hvec[i] - (newX[i][j])*(log10(newX[i][j])/log10(2));
      }
      av_Ni[j] = av_Ni[j] + newX[i][j];
    }
    Hvec[i] = Hvec[i]/maxH;

    Hbar = Hbar + Hvec[i];
  }
  Hbar = Hbar/dnrow;

  mu_avNi = 0.0;
  for (i=0; i<ncol; i++) {
    mu_avNi = mu_avNi + av_Ni[i];
  }
  mu_avNi = mu_avNi/dncol;

  var_avNi = 0.0;
  for (i=0; i<ncol; i++) {
    var_avNi = var_avNi + (av_Ni[i] - mu_avNi)*(av_Ni[i] - mu_avNi);
  }
  var_avNi = var_avNi/dncol;
  sd_avNi = sqrt(var_avNi);


  // set up output structure


  (*full).X = newX;
  (*full).cumX = cumX;
  (*full).Hvec = Hvec;
  (*full).Hbar = Hbar;
  (*full).av_Ni = av_Ni;
  (*full).muNi = mu_avNi;
  (*full).sdNi = sd_avNi;

  return full;

}


XMAT* Xhap(double **Xmat, int nrow, int ncol) {

  int i,j;
  double new, dcheck, dncol, dnrow, maxH, Hbar;
  double mu_avNi, sd_avNi, var_avNi;

  double *Hvec=NULL, *av_Ni=NULL;
  double **newX=NULL, **cumX=NULL;

  XMAT *full=NULL;


  // memory allocation

  full = (XMAT*)calloc(1,sizeof(XMAT));

  newX = (double**)calloc(nrow,sizeof(double*));
  cumX = (double**)calloc(nrow,sizeof(double*));
  for (i=0; i<nrow; i++) {
    newX[i] = (double*)calloc(ncol,sizeof(double));
    cumX[i] = (double*)calloc(ncol,sizeof(double));
  }


  // convert founder strain expectations to probabilities and round

  for (j=0; j<nrow; j++) {
    for (i=0; i<ncol; i++) {
      new = Xmat[j][i];
      new = NRroundit(new,4);
      newX[j][i] = new;
    }
  }


  // normalise probability vectors and calculate entropy, etc

  Hvec = (double*)calloc(nrow,sizeof(double));
  av_Ni = (double*)calloc(ncol,sizeof(double));

  dncol = (double)ncol;
  ncol = (int)ncol;

  dnrow = (double)nrow;
  nrow = (int)nrow;

  maxH = 0.0;
  for (j=0; j<ncol; j++) {
    maxH = maxH - (1.0/dncol)*(log10((1.0/dncol))/log10(2));
  }

  Hbar = 0.0;
  for (i=0; i<nrow; i++) {

    dcheck = 0.0;
    for (j=0; j<ncol; j++) {
      dcheck = dcheck + newX[i][j];
    }

    for (j=0; j<ncol; j++) {
      newX[i][j] = newX[i][j]/dcheck;
    }

    dcheck = 0.0;
    for (j=0; j<ncol; j++) {
      dcheck = dcheck + newX[i][j];
      cumX[i][j] = dcheck;
    }
    if ( (dcheck <=0.99999) || (dcheck >= 1.00001) ) {
      Rprintf("individual %i : dcheck = %e ERROR HMM probs do not sum to 1\n",i,dcheck);
    }

    Hvec[i] = 0.0;
    for (j=0; j<ncol; j++) {
      if ( newX[i][j] != 0.0 ) {
        Hvec[i] = Hvec[i] - (newX[i][j])*(log10(newX[i][j])/log10(2));
      }
      av_Ni[j] = av_Ni[j] + newX[i][j];
    }
    Hvec[i] = Hvec[i]/maxH;

    Hbar = Hbar + Hvec[i];
  }
  Hbar = Hbar/dnrow;

  mu_avNi = 0.0;
  for (i=0; i<ncol; i++) {
    mu_avNi = mu_avNi + av_Ni[i];
  }
  mu_avNi = mu_avNi/dncol;

  var_avNi = 0.0;
  for (i=0; i<ncol; i++) {
    var_avNi = var_avNi + (av_Ni[i] - mu_avNi)*(av_Ni[i] - mu_avNi);
  }
  var_avNi = var_avNi/dncol;
  sd_avNi = sqrt(var_avNi);


  // set up output structure


  (*full).X = newX;
  (*full).cumX = cumX;
  (*full).Hvec = Hvec;
  (*full).Hbar = Hbar;
  (*full).av_Ni = av_Ni;
  (*full).muNi = mu_avNi;
  (*full).sdNi = sd_avNi;

  return full;
}


double NRroundit(double d, int dig) {

  double m = powl(10,dig);
  return round(d*m)/m;
}



#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float ran2(long *idum) {

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {                 // initialise
    if (-(*idum) < 1) *idum=1;      // be sure to prevent idum=0
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7; j>=0; j--) {     // load the shuffle table (after 8 warm-ups)
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }

  k=(*idum)/IQ1;                    // start here when not initialising
  *idum=IA1*(*idum-k*IQ1)-k*IR1;    // compute idum = (IA1*idum) % IM1 without overflows, by Schrage's \
method
  if (*idum < 0) *idum += IM1;
 k=idum2/IQ2;
 idum2=IA2*(idum2-k*IQ2)-k*IR2;    // compute idum2 = (IA2*idum) % IM2 likewise
 if (idum2 < 0) idum2 += IM2;
 j=iy/NDIV;                        // will be in the range 0..NTAB-1
 iy=iv[j]-idum2;                   // here idum is shuffled, idum and idum2 are combined to generate o\
utput
  iv[j] = *idum;
 if (iy < 1) iy += IMM1;
 if ((temp=AM*iy) > RNMX) return RNMX;
 else return temp;                 // because users don't expect endpoint values

}

