

#define _ISOC99_SOURCE
#define _GNU_SOURCE

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<search.h>
#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<Rmath.h>

#include"happy.h"
/*#include"cl.h"
#include"cmp.h"
#include"stats.h"
*/

#include "hbcore.h"
#include "hbrem.h"



SEXP hbrem( SEXP RX, SEXP HaploidInd, SEXP Ndip, SEXP Nind, SEXP Npost, SEXP Nbin, SEXP Ry ) {


  SEXP Rstats = R_NilValue;
  SEXP Reffectmeans = R_NilValue;
  SEXP ReffectSDs = R_NilValue;
  SEXP ReffectNis = R_NilValue;
  SEXP list = R_NilValue;

  int haploid_ind,nbin,ncol,nrow,nsim;

  double *y=NULL;
  double *pntRstats=NULL, *pntReffectmeans=NULL, *pntReffectSDs=NULL;
  double *pntReffectNis=NULL, *pntRX=NULL;


  PROTECT(RX = AS_NUMERIC(RX));
  pntRX = (double*)NUMERIC_POINTER(RX);

  PROTECT(HaploidInd = AS_INTEGER(HaploidInd));
  haploid_ind = (int)INTEGER_POINTER(HaploidInd)[0];

  PROTECT(Ndip = AS_INTEGER(Ndip));
  ncol = (int)INTEGER_POINTER(Ndip)[0];

  PROTECT(Nind = AS_INTEGER(Nind));
  nrow = (int)INTEGER_POINTER(Nind)[0];

  PROTECT(Npost = AS_INTEGER(Npost));
  nsim = (int)INTEGER_POINTER(Npost)[0];

  PROTECT(Nbin = AS_INTEGER(Nbin));
  nbin = (int)INTEGER_POINTER(Nbin)[0];

  /*  Rprintf( "nbin %d npost %d nrow %d ncol %d haploid_ind %d\n", nbin, nsim, nrow, ncol, haploid_ind); */
  PROTECT(Ry = AS_NUMERIC(Ry));
  y = (double*)NUMERIC_POINTER(Ry);

  PROTECT(Rstats = NEW_NUMERIC(45));  // allocates storage space
  pntRstats = (double*)NUMERIC_POINTER(Rstats);

  PROTECT(Reffectmeans = NEW_NUMERIC(ncol));   // allocates storage space
  pntReffectmeans = (double*)NUMERIC_POINTER(Reffectmeans);

  PROTECT(ReffectSDs = NEW_NUMERIC(ncol));   // allocates storage space
  pntReffectSDs = (double*)NUMERIC_POINTER(ReffectSDs);

  PROTECT(ReffectNis = NEW_NUMERIC(ncol));   // allocates storage space
  pntReffectNis = (double*)NUMERIC_POINTER(ReffectNis);


  int i,j,inc;
  long idum=0;
  float init;

  double **Xmat=NULL;

  XMAT *xmat=NULL;

  POST *SLsamp=NULL;
  SSTA *SLstats=NULL;


  Xmat = (double**)calloc(nrow,sizeof(double*));
  for (j=0; j<nrow; j++) {
    Xmat[j] = (double*)calloc(ncol,sizeof(double));
  }


  idum = -7;
  init = ran2(&idum);

  GetRNGstate();


  for (j=0; j<ncol; j++) {

    inc = j*nrow;

    for (i=0; i<nrow; i++) {
      Xmat[i][j] = pntRX[inc + i];
    }
  }


  if ( haploid_ind == 0 ) {

    xmat = Xdip(Xmat,nrow,ncol);

    SLsamp = single_locus_jointpostX(xmat,y,nsim,ncol,nrow,nbin,&idum);
    SLstats = single_locus_sumstatsX(xmat,SLsamp,y,nsim,ncol,nrow,nbin);

    pntRstats[0] = (*xmat).Hbar;
    pntRstats[1] = (*xmat).sdNi;
    pntRstats[2] = (*SLstats).BIC_qtl;
    pntRstats[3] = (*SLstats).BIC_null;
    pntRstats[4] = (*SLstats).BF;
    pntRstats[5] = (*SLstats).logBF;
    pntRstats[6] = (*SLstats).DIC_qtl;
    pntRstats[7] = (*SLstats).DIC_null;
    pntRstats[8] = (*SLstats).DIC_diff;
    pntRstats[9] = (*SLstats).pD_qtl;
    pntRstats[10] = (*SLstats).pD_null;
    pntRstats[11] = (*SLstats).modekT;
    pntRstats[12] = (*SLstats).ga;
    pntRstats[13] = (*SLstats).gb;
    pntRstats[14] = (*SLstats).modevar;
    pntRstats[15] = (*SLstats).medkT;
    pntRstats[16] = (*SLstats).medmu;
    pntRstats[17] = (*SLstats).medvar;
    pntRstats[18] = (*SLstats).mukT;
    pntRstats[19] = (*SLstats).mumu;
    pntRstats[20] = (*SLstats).muvar;
    pntRstats[21] = (*SLstats).HPDlkT[0];
    pntRstats[22] = (*SLstats).HPDukT[0];
    pntRstats[23] = (*SLstats).HPDlmu[0];
    pntRstats[24] = (*SLstats).HPDumu[0];
    pntRstats[25] = (*SLstats).HPDlvar[0];
    pntRstats[26] = (*SLstats).HPDuvar[0];
    pntRstats[27] = (*SLstats).HPDlkT[1];
    pntRstats[28] = (*SLstats).HPDukT[1];
    pntRstats[29] = (*SLstats).HPDlmu[1];
    pntRstats[30] = (*SLstats).HPDumu[1];
    pntRstats[31] = (*SLstats).HPDlvar[1];
    pntRstats[32] = (*SLstats).HPDuvar[1];
    pntRstats[33] = (*SLstats).HPDlkT[2];
    pntRstats[34] = (*SLstats).HPDukT[2];
    pntRstats[35] = (*SLstats).HPDlmu[2];
    pntRstats[36] = (*SLstats).HPDumu[2];
    pntRstats[37] = (*SLstats).HPDlvar[2];
    pntRstats[38] = (*SLstats).HPDuvar[2];
    pntRstats[39] = (*SLstats).HPDlkT[3];
    pntRstats[40] = (*SLstats).HPDukT[3];
    pntRstats[41] = (*SLstats).HPDlmu[3];
    pntRstats[42] = (*SLstats).HPDumu[3];
    pntRstats[43] = (*SLstats).HPDlvar[3];
    pntRstats[44] = (*SLstats).HPDuvar[3];

    for (i=0; i<ncol; i++) {
      pntReffectmeans[i] = (*SLstats).muT[i];
      pntReffectSDs[i] = (*SLstats).sdT[i];
      pntReffectNis[i] = (*xmat).av_Ni[i];
    }


    free((*SLsamp).kT);
    free((*SLsamp).var);
    free((*SLsamp).mu);
    for (j=0; j<nsim; j++) {
      free((*SLsamp).T[j]);
      free((*SLsamp).yTbarX[j]);
      free((*SLsamp).NiX[j]);
    }
    free((*SLsamp).T);
    free((*SLsamp).yTbarX);
    free((*SLsamp).NiX);
    free((*SLsamp).null_var);
    free((*SLsamp).null_mu);
    free((*SLsamp).lik_null);
    free((*SLsamp).lik_qtl);
    free(SLsamp);

    free((*SLstats).kThist);
    free((*SLstats).HPDlkT);
    free((*SLstats).HPDlmu);
    free((*SLstats).HPDlvar);
    free((*SLstats).HPDukT);
    free((*SLstats).HPDumu);
    free((*SLstats).HPDuvar);
    free((*SLstats).muT);
    free((*SLstats).sdT);
    free(SLstats);

    for (j=0; j<nrow; j++) {
      free((*xmat).X[j]);
      free((*xmat).cumX[j]);
    }
    free((*xmat).X);
    free((*xmat).cumX);
    free((*xmat).Hvec);
    free((*xmat).av_Ni);
    free(xmat);


  } else if ( haploid_ind == 1 ) {

    xmat = Xhap(Xmat,nrow,ncol);

    SLsamp = single_locus_jointpostX(xmat,y,nsim,ncol,nrow,nbin,&idum);
    SLstats = single_locus_sumstatsX(xmat,SLsamp,y,nsim,ncol,nrow,nbin);

    pntRstats[0] = (*xmat).Hbar;
    pntRstats[1] = (*xmat).sdNi;
    pntRstats[2] = (*SLstats).BIC_qtl;
    pntRstats[3] = (*SLstats).BIC_null;
    pntRstats[4] = (*SLstats).BF;
    pntRstats[5] = (*SLstats).logBF;
    pntRstats[6] = (*SLstats).DIC_qtl;
    pntRstats[7] = (*SLstats).DIC_null;
    pntRstats[8] = (*SLstats).DIC_diff;
    pntRstats[9] = (*SLstats).pD_qtl;
    pntRstats[10] = (*SLstats).pD_null;
    pntRstats[11] = (*SLstats).modekT;
    pntRstats[12] = (*SLstats).ga;
    pntRstats[13] = (*SLstats).gb;
    pntRstats[14] = (*SLstats).modevar;
    pntRstats[15] = (*SLstats).medkT;
    pntRstats[16] = (*SLstats).medmu;
    pntRstats[17] = (*SLstats).medvar;
    pntRstats[18] = (*SLstats).mukT;
    pntRstats[19] = (*SLstats).mumu;
    pntRstats[20] = (*SLstats).muvar;
    pntRstats[21] = (*SLstats).HPDlkT[0];
    pntRstats[22] = (*SLstats).HPDukT[0];
    pntRstats[23] = (*SLstats).HPDlmu[0];
    pntRstats[24] = (*SLstats).HPDumu[0];
    pntRstats[25] = (*SLstats).HPDlvar[0];
    pntRstats[26] = (*SLstats).HPDuvar[0];
    pntRstats[27] = (*SLstats).HPDlkT[1];
    pntRstats[28] = (*SLstats).HPDukT[1];
    pntRstats[29] = (*SLstats).HPDlmu[1];
    pntRstats[30] = (*SLstats).HPDumu[1];
    pntRstats[31] = (*SLstats).HPDlvar[1];
    pntRstats[32] = (*SLstats).HPDuvar[1];
    pntRstats[33] = (*SLstats).HPDlkT[2];
    pntRstats[34] = (*SLstats).HPDukT[2];
    pntRstats[35] = (*SLstats).HPDlmu[2];
    pntRstats[36] = (*SLstats).HPDumu[2];
    pntRstats[37] = (*SLstats).HPDlvar[2];
    pntRstats[38] = (*SLstats).HPDuvar[2];
    pntRstats[39] = (*SLstats).HPDlkT[3];
    pntRstats[40] = (*SLstats).HPDukT[3];
    pntRstats[41] = (*SLstats).HPDlmu[3];
    pntRstats[42] = (*SLstats).HPDumu[3];
    pntRstats[43] = (*SLstats).HPDlvar[3];
    pntRstats[44] = (*SLstats).HPDuvar[3];

    for (i=0; i<ncol; i++) {
      pntReffectmeans[i] = (*SLstats).muT[i];
      pntReffectSDs[i] = (*SLstats).sdT[i];
      pntReffectNis[i] = (*xmat).av_Ni[i];
    }

    free((*SLsamp).kT);
    free((*SLsamp).var);
    free((*SLsamp).mu);
    for (j=0; j<nsim; j++) {
      free((*SLsamp).T[j]);
      free((*SLsamp).yTbarX[j]);
      free((*SLsamp).NiX[j]);
    }
    free((*SLsamp).T);
    free((*SLsamp).yTbarX);
    free((*SLsamp).NiX);
    free((*SLsamp).null_var);
    free((*SLsamp).null_mu);
    free((*SLsamp).lik_null);
    free((*SLsamp).lik_qtl);
    free(SLsamp);

    free((*SLstats).kThist);
    free((*SLstats).HPDlkT);
    free((*SLstats).HPDlmu);
    free((*SLstats).HPDlvar);
    free((*SLstats).HPDukT);
    free((*SLstats).HPDumu);
    free((*SLstats).HPDuvar);
    free((*SLstats).muT);
    free((*SLstats).sdT);
    free(SLstats);

    for (j=0; j<nrow; j++) {
      free((*xmat).X[j]);
      free((*xmat).cumX[j]);
    }
    free((*xmat).X);
    free((*xmat).cumX);
    free((*xmat).Hvec);
    free((*xmat).av_Ni);
    free(xmat);

  }   // end of HAPLOID loop


  PutRNGstate();


  for (j=0; j<nrow; j++) {
    free(Xmat[j]);
  }
  free(Xmat);


  // create list structure to return

  PROTECT(list = allocVector(VECSXP,4));  // create list of 2 vectors

  SET_VECTOR_ELT(list, 0, Rstats);   //  attach stats vector to list
  SET_VECTOR_ELT(list, 1, Reffectmeans);  // attach effects vector to list
  SET_VECTOR_ELT(list, 2, ReffectSDs);  // attach effects vector to list
  SET_VECTOR_ELT(list, 3, ReffectNis);


  UNPROTECT(12);

  return(list);

}
