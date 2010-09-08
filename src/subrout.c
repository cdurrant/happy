#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>

#include"hbcore.h"
#include"hbrem.h"
#include"cmp.h"

POST* single_locus_jointpostX(XMAT *xmat,double *y,int nsim,int ncol,int nrow,int nbin,long *idum) {

  FILE *out=NULL;
  char check;
  char sampoutname[25]="SLXsamp_R.dat";

  int i,sim;
  double nu,kT;

  int **NiX=NULL;

  double *postkT=NULL, *postvar=NULL, *postmu=NULL, *nullvar=NULL, *nullmu=NULL;
  double *Lnull=NULL, *Lqtl=NULL;

  double **postT=NULL, **yTbarX=NULL;

  SIMX *selX=NULL;
  GRKT *kTdist=NULL;
  POST *single=NULL;


  single = (POST*)calloc(1,sizeof(POST));

  postkT = (double*)calloc(nsim,sizeof(double));
  postvar = (double*)calloc(nsim,sizeof(double));
  postmu = (double*)calloc(nsim,sizeof(double));

  postT = (double**)calloc(nsim,sizeof(double*));
  yTbarX = (double**)calloc(nsim,sizeof(double*));
  for (i=0; i<nsim; i++) {
    postT[i] = (double*)calloc(ncol,sizeof(double));
    yTbarX[i] = (double*)calloc(ncol,sizeof(double));
  }

  nullvar = (double*)calloc(nsim,sizeof(double));
  nullmu = (double*)calloc(nsim,sizeof(double));

  Lnull = (double*)calloc(nsim,sizeof(double));
  Lqtl = (double*)calloc(nsim,sizeof(double));

  NiX = (int**)calloc(nsim,sizeof(double));
  for (i=0; i<nsim; i++) {
    NiX[i] = (int*)calloc(ncol,sizeof(int));
  }


  for (sim=1; sim<=nsim; sim++) {

    selX = drawX(xmat,ncol,nrow,idum);
    kTdist = truegridkT(selX,y,ncol,nrow,nbin);

    for (i=0; i<ncol; i++) {
      NiX[sim-1][i] = (*selX).Ni[i];
      yTbarX[sim-1][i] = (*kTdist).yTbar[i];
    }

    postkT[sim-1] = drawkT(kTdist,idum);

    nu = ((*kTdist).N - 1.0);
    postvar[sim-1] = draw_knownvar(kTdist,(*selX).Ni,ncol,postkT[sim-1],nu,nbin);

    postmu[sim-1] = draw_knownmu(kTdist,(*selX).Ni,ncol,postkT[sim-1],postvar[sim-1],nbin);

    for (i=0; i<ncol; i++) {
      postT[sim-1][i] = draw_knownTi(kTdist,(*selX).Ni,postkT[sim-1],postvar[sim-1],postmu[sim-1],nbin,i);
    }


    nullvar[sim-1] = draw_nullvar(selX,y,nrow,nbin);
    nullmu[sim-1] = draw_nullmu(selX,y,nrow,nullvar[sim-1],nbin);


    Lqtl[sim-1] = qtl_lik(selX,y,postkT[sim-1],postvar[sim-1],postmu[sim-1],postT[sim-1],nrow,nbin);
    Lnull[sim-1] = null_lik(selX,y,nullvar[sim-1],nullmu[sim-1],nrow,nbin);


    free((*selX).X);
    free((*selX).Ni);
    free(selX);

    free((*kTdist).pdf);
    free((*kTdist).cdf);
    free((*kTdist).yTbar);
    free(kTdist);

  }  //  end of simulation loop



  //  out = fopen(sampoutname,"w");
  //  if (out == NULL) {
  //    printf("output file opening failed\n");
  //    exit(1);
  //  }

  //  fprintf(out," kT var mu T[1] T[2] T[3] T[4] T[5] T[6] T[7] T[8] T[9] T[10] T[11] T[12] T[13] T[14] T[15] T[16] T[17] T[18] T[19] yT[1] yT[2] yT[3] yT[4] yT[5] yT[6] yT[7] yT[8] yT[9] yT[10] yT[11] yT[12] yT[13] yT[14] yT[15] yT[16] yT[17] yT[18] yT[19] Ni[1] Ni[2] Ni[3] Ni[4] Ni[5] Ni[6] Ni[7] Ni[8] Ni[9] Ni[10] Ni[11] Ni[12] Ni[13] Ni[14] Ni[15] Ni[16] Ni[17] Ni[18] Ni[19]\n");

  //  for (sim=0; sim<nsim; sim++) {
  //    fprintf(out," %i  %f %f %f    %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",(sim+1),postkT[sim],postvar[sim],postmu[sim],postT[sim][0],postT[sim][1],postT[sim][2],postT[sim][3],postT[sim][4],postT[sim][5],postT[sim][6],postT[sim][7],postT[sim][8],postT[sim][9],postT[sim][10],postT[sim][11],postT[sim][12],postT[sim][13],postT[sim][14],postT[sim][15],postT[sim][16],postT[sim][17],postT[sim][18],yTbarX[sim][0],yTbarX[sim][1],yTbarX[sim][2],yTbarX[sim][3],yTbarX[sim][4],yTbarX[sim][5],yTbarX[sim][6],yTbarX[sim][7],yTbarX[sim][8],yTbarX[sim][9],yTbarX[sim][10],yTbarX[sim][11],yTbarX[sim][12],yTbarX[sim][13],yTbarX[sim][14],yTbarX[sim][15],yTbarX[sim][16],yTbarX[sim][17],yTbarX[sim][18],NiX[sim][0],NiX[sim][1],NiX[sim][2],NiX[sim][3],NiX[sim][4],NiX[sim][5],NiX[sim][6],NiX[sim][7],NiX[sim][8],NiX[sim][9],NiX[sim][10],NiX[sim][11],NiX[sim][12],NiX[sim][13],NiX[sim][14],NiX[sim][15],NiX[sim][16],NiX[sim][17],NiX[sim][18]);
  //  }

  //  check = fclose(out);
  //  if (check != 0) {
  //    printf("output file closing failed: function returned %i\n",check);
  //    exit(1);
  //  }



  (*single).kT = postkT;
  (*single).var = postvar;
  (*single).mu = postmu;
  (*single).T = postT;
  (*single).yTbarX = yTbarX;
  (*single).NiX = NiX;
  (*single).null_var = nullvar;
  (*single).null_mu = nullmu;
  (*single).lik_null = Lnull;
  (*single).lik_qtl = Lqtl;

  return single;

}


POST* single_locus_jointpost(SIMX *trueX,double *y,int nsim,int ncol,int nrow,int nbin,long *idum) {

  FILE *out=NULL;
  char check;
  char sampoutname[25]="SLsamp_R.dat";

  int i,sim;
  double nu, kTN, kTt;

  double *postkT=NULL, *postvar=NULL, *postmu=NULL, *Lqtl=NULL;
  double *nullvar=NULL, *nullmu=NULL, *Lnull=NULL;
  double **postT=NULL;

  GRKT *truekTdist=NULL;
  POST *single=NULL;


  single=(POST*)calloc(1,sizeof(POST));

  postkT = (double*)calloc(nsim,sizeof(double));
  postvar = (double*)calloc(nsim,sizeof(double));
  postmu = (double*)calloc(nsim,sizeof(double));

  postT=(double**)calloc(nsim,sizeof(double*));
  for (i=0; i<nsim; i++) {
    postT[i] = (double*)calloc(ncol,sizeof(double));
  }

  nullvar = (double*)calloc(nsim,sizeof(double));
  nullmu = (double*)calloc(nsim,sizeof(double));

  Lqtl = (double*)calloc(nsim,sizeof(double));
  Lnull = (double*)calloc(nsim,sizeof(double));


  // calculate grid pdf and cdf for kT
  truekTdist = truegridkT(trueX,y,ncol,nrow,nbin);

  kTt = (*truekTdist).t;
  kTN = (*truekTdist).N;

  for (sim=1; sim<=nsim; sim++) {

    postkT[sim-1] = drawkT(truekTdist,idum);

    //    nu=(nrow-1);
    nu = (kTN - 1.0);
    postvar[sim-1] = draw_knownvar(truekTdist,(*trueX).Ni,ncol,postkT[sim-1],nu,nbin);

    postmu[sim-1] = draw_knownmu(truekTdist,(*trueX).Ni,ncol,postkT[sim-1],postvar[sim-1],nbin);

    for (i=0; i<ncol; i++) {
      postT[sim-1][i] = draw_knownTi(truekTdist,(*trueX).Ni,postkT[sim-1],postvar[sim-1],postmu[sim-1],nbin,i);
    }

    nullvar[sim-1] = draw_nullvar(trueX,y,nrow,nbin);
    nullmu[sim-1] = draw_nullmu(trueX,y,nrow,nullvar[sim-1],nbin);

    Lqtl[sim-1] = qtl_lik(trueX,y,postkT[sim-1],postvar[sim-1],postmu[sim-1],postT[sim-1],nrow,nbin);
    Lnull[sim-1] = null_lik(trueX,y,nullvar[sim-1],nullmu[sim-1],nrow,nbin);

  }  //  end of simulation loop


  //  out = fopen(sampoutname,"w");
  //  if (out == NULL) {
  //    Rprintf("output file opening failed\n");
  //    exit(1);
  //  }

  //  fprintf(out," kT var mu T[1] T[2] T[3] T[4] T[5] T[6] T[7] T[8] T[9] T[10] T[11] T[12] T[13] T[14] T[15] T[16] T[17] T[18] T[19]\n");

  //  for (sim=0; sim<nsim; sim++) {
  //    fprintf(out," %i  %f %f %f  %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",(sim+1),postkT[sim],postvar[sim],postmu[sim],postT[sim][0],postT[sim][1],postT[sim][2],postT[sim][3],postT[sim][4],postT[sim][5],postT[sim][6],postT[sim][7],postT[sim][8],postT[sim][9],postT[sim][10],postT[sim][11],postT[sim][12],postT[sim][13],postT[sim][14],postT[sim][15],postT[sim][16],postT[sim][17],postT[sim][18]);
  //  }

  //  check = fclose(out);
  //  if (check != 0) {
  //    Rprintf("output file closing failed: function returned %i\n",check);
  //    exit(1);
  //  }



  (*single).kT = postkT;
  (*single).var = postvar;
  (*single).mu = postmu;
  (*single).T = postT;
  (*single).t = kTt;
  (*single).N = kTN;
  (*single).null_var = nullvar;
  (*single).null_mu = nullmu;
  (*single).lik_null = Lnull;
  (*single).lik_qtl = Lqtl;


  free((*truekTdist).pdf);
  free((*truekTdist).cdf);
  free((*truekTdist).yTbar);
  free(truekTdist);


  return single;

}


SSTA* single_locus_sumstats(POST *SLsamp, SIMX *trueX, double *y, int nsim, int ncol, int nrow, int nbin) {

  int i,j,index,count1,nprc,sum;

  double cred,dnrow,dnsim,prc;
  double mukT,medkT,mumu,medmu,muvar,medvar,modevar,modekT;
  double BF,BIC_qtl,BIC_null,logBF,DIC_null,DIC_qtl,DIC_diff;
  double avlik_qtl,avlik_null,Lqtl,Lnull,Lfocqtl,pD_null,pD_qtl;
  double ga,gb,varkT,mu_nullmu,mu_nullvar,mode_nullvar;
  double bkT,bmu,bvar,tkT,tmu,tvar,kTdiff,mudiff,vardiff;
  double kTb99,kTt99,kTb95,kTt95,kTb75,kTt75,kTb50,kTt50;
  double mub99,mut99,mub95,mut95,mub75,mut75,mub50,mut50;
  double varb99,vart99,varb95,vart95,varb75,vart75,varb50,vart50;

  int *kThst=NULL;
  double *stkT=NULL,*HPDlkT=NULL,*HPDukT=NULL;
  double *stmu=NULL,*HPDlmu=NULL,*HPDumu=NULL;
  double *stvar=NULL,*HPDlvar=NULL,*HPDuvar=NULL;
  double *muT=NULL,*sdT=NULL;

  SSTA *stats=NULL;


  stats = (SSTA*)calloc(1,sizeof(SSTA));

  HPDlkT = (double*)calloc(4,sizeof(double));
  HPDukT = (double*)calloc(4,sizeof(double));

  HPDlmu = (double*)calloc(4,sizeof(double));
  HPDumu = (double*)calloc(4,sizeof(double));

  HPDlvar = (double*)calloc(4,sizeof(double));
  HPDuvar = (double*)calloc(4,sizeof(double));

  kThst = (int*)calloc(1001,sizeof(int));

  stkT = (double*)calloc((nsim+1),sizeof(double));
  stmu = (double*)calloc((nsim+1),sizeof(double));
  stvar = (double*)calloc((nsim+1),sizeof(double));

  muT = (double*)calloc(ncol,sizeof(double));
  sdT = (double*)calloc(ncol,sizeof(double));


  nprc = 201;
  prc = 200.0;

  dnsim = (double)nsim;
  nsim = (int)nsim;

  
  // set up data vectors, histogram for kT, calculate means

  mukT = 0.0;
  mumu = 0.0;
  muvar = 0.0;
  mu_nullmu = 0.0;
  mu_nullvar = 0.0;
  avlik_null = 0.0;
  avlik_qtl = 0.0;

  for (i=0; i<nsim; i++) {

    stkT[i+1] = (*SLsamp).kT[i];
    mukT = mukT + (*SLsamp).kT[i];
    index = (int)((*SLsamp).kT[i]*prc);
    kThst[index] = kThst[index] + 1;

    stmu[i+1] = (*SLsamp).mu[i];
    mumu = mumu + (*SLsamp).mu[i];

    stvar[i+1] = (*SLsamp).var[i];
    muvar = muvar + (*SLsamp).var[i];

    for (j=0; j<ncol; j++) {
      if ( (*trueX).Ni[j] >= nbin ) {
	muT[j] = muT[j] + (*SLsamp).T[i][j];
      }
    }

    mu_nullmu = mu_nullmu + (*SLsamp).null_mu[i];
    mu_nullvar = mu_nullvar + (*SLsamp).null_var[i];

    avlik_null = avlik_null + (*SLsamp).lik_null[i];
    avlik_qtl = avlik_qtl + (*SLsamp).lik_qtl[i];
  }
  mukT = mukT/dnsim;
  mumu = mumu/dnsim;
  muvar = muvar/dnsim;
  for (j=0; j<ncol; j++) {
    if ( (*trueX).Ni[j] >= nbin ) {
      muT[j] = muT[j]/dnsim;
    }
  }
  mu_nullmu = mu_nullmu/dnsim;
  mu_nullvar = mu_nullvar/dnsim;
  avlik_null = avlik_null/dnsim;
  avlik_qtl = avlik_qtl/dnsim;

  modevar = (((*SLsamp).N - 3.0)/((*SLsamp).N + 1.0))*muvar;
  mode_nullvar = (((*SLsamp).N - 3.0)/((*SLsamp).N + 1.0))*mu_nullvar;

  varkT = 0.0;
  for (i=0; i<nsim; i++) {
    varkT = varkT + ((*SLsamp).kT[i] - mukT)*((*SLsamp).kT[i] - mukT);
    for (j=0; j<ncol; j++) {
      if ( (*trueX).Ni[j] >= nbin ) {
	sdT[j] = sdT[j] + ((*SLsamp).T[i][j] - muT[j])*((*SLsamp).T[i][j] - muT[j]);
      }
    }
  }
  varkT = varkT/dnsim;
  for (j=0; j<ncol; j++) {
    if ( (*trueX).Ni[j] >= nbin ) {
      sdT[j] = sdT[j]/dnsim;
      sdT[j] = sqrt(sdT[j]);
    }
  }

  ga = ((mukT*mukT*(1.0 - mukT)/varkT) - mukT);
  gb = ga*((1.0 - mukT)/mukT);
  modekT = (ga - 1.0)/(ga + gb - 2.0);
  if (modekT < 0.0) {
    modekT = 0.0;
  } else if (modekT > 1.0) {
    modekT = 1.0;
  }


  // plug-ins for pD

  Lnull = null_lik(trueX,y,mode_nullvar,mu_nullmu,nrow,nbin);
  Lqtl = qtl_lik(trueX,y,modekT,modevar,mumu,muT,nrow,nbin);

  pD_qtl = -2.0*(avlik_qtl - Lqtl);
  pD_null = -2.0*(avlik_null - Lnull);

  DIC_qtl = ( pD_qtl - 2.0*avlik_qtl );
  DIC_null = ( pD_null - 2.0*avlik_null );

  DIC_diff = ( DIC_null - DIC_qtl );

  // plug-ins for BIC

  Lfocqtl = qtl_Lfoc(trueX,y,modekT,modevar,mumu,nrow,ncol,nbin);

  BIC_qtl = ( (3.0*log((*SLsamp).N)) - (2.0*Lfocqtl) );
  BIC_null = ( (2.0*log((*SLsamp).N)) - (2.0*Lnull) );

  BF = exp(-(BIC_null - BIC_qtl)/2.0);
  logBF = -log10(BF);


  count1 = 0;
  for (i=0; i<nprc; i++) {
    if ( kThst[i] > count1 ) {
      count1 = kThst[i];
    }
  }


  // sort posterior sample
  // NRsort doesn't use 0 index - entries start at 1

  qsort( stkT+1, nsim, sizeof(double), dcmp );
  qsort( stmu+1, nsim, sizeof(double), dcmp );
  qsort( stvar+1, nsim, sizeof(double), dcmp );

//  NRsort(nsim,stkT);
//  NRsort(nsim,stmu);
//  NRsort(nsim,stvar);


  // calculate medians

  index = (int)(0.5*((float)nsim));
  cred = stkT[index] + stkT[index + 1];
  medkT = cred/2.0;
  cred = stmu[index] + stmu[index + 1];
  medmu = cred/2.0;
  cred = stvar[index] + stvar[index + 1];
  medvar = cred/2.0;


  // calculate HPD credible intervals

  kTdiff = 1.0;
  mudiff = ( stmu[nsim] -  stmu[1] );
  vardiff = ( stvar[nsim] -  stvar[1] );

  i = 1;
  while ( (i+(99*nsim/100)) < nsim ) {
    bkT = stkT[i];
    tkT = stkT[i+(99*nsim/100)];
    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTb99 = bkT;
      kTt99 = tkT;
    }
    bmu = stmu[i];
    tmu = stmu[i+(99*nsim/100)];
    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub99 = bmu;
      mut99 = tmu;
    }
    bvar = stvar[i];
    tvar = stvar[i+(99*nsim/100)];
    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb99 = bvar;
      vart99 = tvar;
    }
    i = i + 1;
  }


  kTdiff = 1.0;
  mudiff = ( stmu[nsim] -  stmu[1] );
  vardiff = ( stvar[nsim] -  stvar[1] );

  i = 1;
  while ( (i+(95*nsim/100)) < nsim ) {
    bkT = stkT[i];
    tkT = stkT[i+(95*nsim/100)];
    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTt95 = tkT;
      kTb95 = bkT;
    }
    bmu = stmu[i];
    tmu = stmu[i+(95*nsim/100)];
    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub95 = bmu;
      mut95 = tmu;
    }
    bvar = stvar[i];
    tvar = stvar[i+(95*nsim/100)];
    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb95 = bvar;
      vart95 = tvar;
    }
    i = i + 1;
  }


  kTdiff = 1.0;
  mudiff = ( stmu[nsim] -  stmu[1] );
  vardiff = ( stvar[nsim] -  stvar[1] );

  i = 1;
  while ( (i+(75*nsim/100)) < nsim ) {
    bkT = stkT[i];
    tkT = stkT[i+(75*nsim/100)];
    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTt75 = tkT;
      kTb75 = bkT;
    }
    bmu = stmu[i];
    tmu = stmu[i+(75*nsim/100)];
    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub75 = bmu;
      mut75 = tmu;
    }
    bvar = stvar[i];
    tvar = stvar[i+(75*nsim/100)];
    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb75 = bvar;
      vart75 = tvar;
    }
    i = i + 1;
  }


  kTdiff = 1.0;
  mudiff = ( stmu[nsim] -  stmu[1] );
  vardiff = ( stvar[nsim] -  stvar[1] );

  i = 1;
  while ( (i+(50*nsim/100)) < nsim ) {
    bkT = stkT[i];
    tkT = stkT[i+(50*nsim/100)];
    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTt50 = tkT;
      kTb50 = bkT;
    }
    bmu = stmu[i];
    tmu = stmu[i+(50*nsim/100)];
    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub50 = bmu;
      mut50 = tmu;
    }
    bvar = stvar[i];
    tvar = stvar[i+(50*nsim/100)];
    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb50 = bvar;
      vart50 = tvar;
    }
    i = i + 1;
  }



  HPDlkT[3] = kTb99;
  HPDukT[3] = kTt99;
  HPDlkT[2] = kTb95;
  HPDukT[2] = kTt95;
  HPDlkT[1] = kTb75;
  HPDukT[1] = kTt75;
  HPDlkT[0] = kTb50;
  HPDukT[0] = kTt50;

  HPDlmu[3] = mub99;
  HPDumu[3] = mut99;
  HPDlmu[2] = mub95;
  HPDumu[2] = mut95;
  HPDlmu[1] = mub75;
  HPDumu[1] = mut75;
  HPDlmu[0] = mub50;
  HPDumu[0] = mut50;

  HPDlvar[3] = varb99;
  HPDuvar[3] = vart99;
  HPDlvar[2] = varb95;
  HPDuvar[2] = vart95;
  HPDlvar[1] = varb75;
  HPDuvar[1] = vart75;
  HPDlvar[0] = varb50;
  HPDuvar[0] = vart50;


  free(stkT);
  free(stmu);
  free(stvar);

  (*stats).kThist = kThst;
  (*stats).HPDlkT = HPDlkT;
  (*stats).HPDukT = HPDukT;
  (*stats).modekT = modekT;
  (*stats).mukT = mukT;
  (*stats).medkT = medkT;
  (*stats).ga = ga;
  (*stats).gb = gb;

  (*stats).HPDlmu = HPDlmu;
  (*stats).HPDumu = HPDumu;
  (*stats).mumu = mumu;
  (*stats).medmu = medmu;

  (*stats).HPDlvar = HPDlvar;
  (*stats).HPDuvar = HPDuvar;
  (*stats).modevar = modevar;
  (*stats).muvar = muvar;
  (*stats).medvar = medvar;

  (*stats).muT = muT;
  (*stats).sdT = sdT;

  (*stats).pD_null = pD_null;
  (*stats).pD_qtl = pD_qtl;
  (*stats).DIC_null = DIC_null;
  (*stats).DIC_qtl = DIC_qtl;
  (*stats).DIC_diff = DIC_diff;
  (*stats).BIC_null = BIC_null;
  (*stats).BIC_qtl = BIC_qtl;
  (*stats).BF = BF;
  (*stats).logBF = logBF;


  return stats;

}


SSTA* single_locus_sumstatsX(XMAT *xmat, POST *SLsamp, double *y, int nsim, int ncol, int nrow, int nbin) {

  int i,j,index,count1,nprc,sum;

  double cred,dnbin,dnrow,dnsim,prc;

  double mukT,medkT,mumu,medmu,muvar,medvar,modevar,modekT;
  double BF,BIC_qtl,BIC_null,logBF,DIC_null,DIC_qtl,DIC_diff;
  double avlik_qtl,avlik_null,Lfocqtl,Lqtl,Lnull,pD_null,pD_qtl;
  double ga,gb,varkT,mu_nullmu,mu_nullvar,mode_nullvar,av_ysq,av_ybar;

  double bkT,bmu,bvar,tkT,tmu,tvar,kTdiff,mudiff,vardiff;
  double kTb99,kTt99,kTb95,kTt95,kTb75,kTt75,kTb50,kTt50;
  double mub99,mut99,mub95,mut95,mub75,mut75,mub50,mut50;
  double varb99,vart99,varb95,vart95,varb75,vart75,varb50,vart50;

  int *kThst=NULL;

  double *stkT=NULL,*HPDlkT=NULL,*HPDukT=NULL;
  double *stmu=NULL,*HPDlmu=NULL,*HPDumu=NULL;
  double *stvar=NULL,*HPDlvar=NULL,*HPDuvar=NULL;

  double *muT=NULL,*sdT=NULL,*nT=NULL,*av_yT=NULL;

  SSTA *stats=NULL;


  nprc = 201;
  prc = 200.0;


  stats = (SSTA*)calloc(1,sizeof(SSTA));

  HPDlkT = (double*)calloc(4,sizeof(double));
  HPDukT = (double*)calloc(4,sizeof(double));

  HPDlmu = (double*)calloc(4,sizeof(double));
  HPDumu = (double*)calloc(4,sizeof(double));

  HPDlvar = (double*)calloc(4,sizeof(double));
  HPDuvar = (double*)calloc(4,sizeof(double));

  kThst = (int*)calloc(nprc,sizeof(int));

  stkT = (double*)calloc((nsim+1),sizeof(double));
  stmu = (double*)calloc((nsim+1),sizeof(double));
  stvar = (double*)calloc((nsim+1),sizeof(double));

  muT = (double*)calloc(ncol,sizeof(double));
  sdT = (double*)calloc(ncol,sizeof(double));
  nT = (double*)calloc(ncol, sizeof(double));
  av_yT = (double*)calloc(ncol,sizeof(double));



  dnsim = (double)nsim;
  nsim = (int)nsim;

  dnrow = (double)nrow;
  nrow = (int)nrow;

  dnbin = (double)nbin;
  nbin = (int)nbin;


  mukT = 0.0;
  mumu = 0.0;
  muvar = 0.0;
  mu_nullmu = 0.0;
  mu_nullvar = 0.0;
  avlik_null = 0.0;
  avlik_qtl = 0.0;

  for (i=0; i<nsim; i++) {

    stkT[i+1] = (*SLsamp).kT[i];
    mukT = mukT + (*SLsamp).kT[i];
    index = (int)((*SLsamp).kT[i]*prc);
    kThst[index] = kThst[index] + 1;

    stmu[i+1] = (*SLsamp).mu[i];
    mumu = mumu + (*SLsamp).mu[i];

    stvar[i+1] = (*SLsamp).var[i];
    muvar = muvar + (*SLsamp).var[i];

    for (j=0; j<ncol; j++) {
      if ( (int)(*SLsamp).NiX[i][j] >= nbin ) {
	muT[j] = muT[j] + (*SLsamp).T[i][j];
	nT[j] = nT[j] + 1.0;
      }
    }

    mu_nullmu = mu_nullmu + (*SLsamp).null_mu[i];
    mu_nullvar = mu_nullvar + (*SLsamp).null_var[i];

    avlik_null = avlik_null + (*SLsamp).lik_null[i];
    avlik_qtl = avlik_qtl + (*SLsamp).lik_qtl[i];
  }

  mukT = mukT/dnsim;
  mumu = mumu/dnsim;
  muvar = muvar/dnsim;
  for(j=0; j<ncol; j++) {
    if ( nT[j] > 0.0 ) {
      muT[j] = muT[j]/nT[j];
    }
  }
  mu_nullmu = mu_nullmu/dnsim;
  mu_nullvar = mu_nullvar/dnsim;
  avlik_null = avlik_null/dnsim;
  avlik_qtl = avlik_qtl/dnsim;


  modevar = ((dnrow - 3.0)/(dnrow + 1.0))*muvar;
  mode_nullvar = ((dnrow - 3.0)/(dnrow + 1.0))*mu_nullvar;

  varkT = 0.0;
  for (i=0; i<nsim; i++) {
    varkT = varkT + ((*SLsamp).kT[i] - mukT)*((*SLsamp).kT[i] - mukT);
    for (j=0; j<ncol; j++) {
      if ( (int)(*SLsamp).NiX[i][j] >= nbin ) {
	sdT[j] = sdT[j] + ((*SLsamp).T[i][j] - muT[j])*((*SLsamp).T[i][j] - muT[j]);
      }
    }
  }
  varkT = varkT/dnsim;
  for (j=0; j<ncol; j++) {
    if ( nT[j] > 0.0 ) {
      sdT[j] = sdT[j]/nT[j];
      sdT[j] = sqrt(sdT[j]);
    }
  }

  ga = ((mukT*mukT*(1.0 - mukT)/varkT) - mukT);
  gb = ga*((1.0 - mukT)/mukT);
  modekT = (ga - 1.0)/(ga + gb - 2.0);
  if (modekT < 0.0) {
    modekT = 0.0;
  } else if (modekT > 1.0) {
    modekT = 1.0;
  }


  // plug-ins for pD

  av_ysq = 0.0;
  av_ybar = 0.0;
  for (i=0; i<nrow; i++) {
    av_ysq = av_ysq + (y[i]*y[i]);
    av_ybar = av_ybar + y[i];
  }
  av_ybar = av_ybar/dnrow;


  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++) {
      av_yT[j] = av_yT[j] + ((*xmat).X[i][j]*y[i]);
    }
  }
  for (j=0; j<ncol; j++) {
    if ( (*xmat).av_Ni[j] > 0.0 ) {
      av_yT[j] = av_yT[j]/(*xmat).av_Ni[j];
    }
  }

  Lnull = null_plug((*xmat).av_Ni,av_ysq,av_ybar,mode_nullvar,mu_nullmu,ncol,nrow);
  Lqtl = qtl_plug(av_yT,(*xmat).av_Ni,av_ysq,modekT,modevar,mumu,muT,ncol,nrow);

  pD_qtl = -2.0*(avlik_qtl - Lqtl);
  pD_null = -2.0*(avlik_null - Lnull);

  DIC_qtl = ( pD_qtl - 2.0*avlik_qtl );
  DIC_null = ( pD_null - 2.0*avlik_null );

  DIC_diff = ( DIC_null - DIC_qtl );

  // plug-ins for BIC

  Lfocqtl = qtl_LfocX(av_yT,(*xmat).av_Ni,av_ysq,av_ybar,modekT,modevar,mumu,ncol,nrow);

  BIC_qtl = ( (3.0*log(dnrow)) - (2.0*Lfocqtl) );
  BIC_null = ( (2.0*log(dnrow)) - (2.0*Lnull) );

  BF = exp(-(BIC_null - BIC_qtl)/2.0);
  logBF = -log10(BF);


  count1 = 0;
  for (i=0; i<nprc; i++) {
    if ( kThst[i] > count1 ) {
      count1 = kThst[i];
    }
  }


  // sort posterior sample
  // NRsort doesn't use 0 index - entries start at 1

  qsort( stkT+1, nsim, sizeof(double), dcmp );
  qsort( stmu+1, nsim, sizeof(double), dcmp );
  qsort( stvar+1, nsim, sizeof(double), dcmp );

//  NRsort(nsim,stkT);
//  NRsort(nsim,stmu);
//  NRsort(nsim,stvar);


  // calculate medians

  index = (int)(0.5*((float)nsim));
  nsim = (int)nsim;

  cred = stkT[index] + stkT[index + 1];
  medkT = cred/2.0;
  cred = stmu[index] + stmu[index + 1];
  medmu = cred/2.0;
  cred = stvar[index] + stvar[index + 1];
  medvar = cred/2.0;


  // calculate HPD credible intervals

  kTdiff = 1.0;
  mudiff = ( stmu[nsim] - stmu[1] );
  vardiff = ( stvar[nsim] - stvar[1] );

  i = 1;
  while ( (i+(99*nsim/100)) < nsim ) {

    bkT = stkT[i];
    tkT = stkT[i+(99*nsim/100)];

    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTb99 = bkT;
      kTt99 = tkT;
    }

    bmu = stmu[i];
    tmu = stmu[i+(99*nsim/100)];

    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub99 = bmu;
      mut99 = tmu;
    }

    bvar = stvar[i];
    tvar = stvar[i+(99*nsim/100)];

    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb99 = bvar;
      vart99 = tvar;
    }

    i = i + 1;
  }


  kTdiff = 1.0;
  mudiff = ( stmu[nsim] - stmu[1] );
  vardiff = ( stvar[nsim] - stvar[1] );

  i = 1;
  while ( (i+(95*nsim/100)) < nsim ) {

    bkT = stkT[i];
    tkT = stkT[i+(95*nsim/100)];

    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTt95 = tkT;
      kTb95 = bkT;
    }

    bmu = stmu[i];
    tmu = stmu[i+(95*nsim/100)];

    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub95 = bmu;
      mut95 = tmu;
    }

    bvar = stvar[i];
    tvar = stvar[i+(95*nsim/100)];

    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb95 = bvar;
      vart95 = tvar;
    }

    i = i + 1;
  }


  kTdiff = 1.0;
  mudiff = ( stmu[nsim] - stmu[1] );
  vardiff = ( stvar[nsim] - stvar[1] );

  i = 1;
  while ( (i+(75*nsim/100)) < nsim ) {

    bkT = stkT[i];
    tkT = stkT[i+(75*nsim/100)];

    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTt75 = tkT;
      kTb75 = bkT;
    }

    bmu = stmu[i];
    tmu = stmu[i+(75*nsim/100)];

    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub75 = bmu;
      mut75 = tmu;
    }

    bvar = stvar[i];
    tvar = stvar[i+(75*nsim/100)];

    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb75 = bvar;
      vart75 = tvar;
    }

    i = i + 1;
  }


  kTdiff = 1.0;
  mudiff = ( stmu[nsim] - stmu[1] );
  vardiff = ( stvar[nsim] - stvar[1] );

  i = 1;
  while ( (i+(50*nsim/100)) < nsim ) {

    bkT = stkT[i];
    tkT = stkT[i+(50*nsim/100)];

    if ( (tkT - bkT) < kTdiff ) {
      kTdiff = (tkT - bkT);
      kTt50 = tkT;
      kTb50 = bkT;
    }

    bmu = stmu[i];
    tmu = stmu[i+(50*nsim/100)];

    if ( (tmu - bmu) < mudiff ) {
      mudiff = (tmu - bmu);
      mub50 = bmu;
      mut50 = tmu;
    }

    bvar = stvar[i];
    tvar = stvar[i+(50*nsim/100)];

    if ( (tvar - bvar) < vardiff ) {
      vardiff = (tvar - bvar);
      varb50 = bvar;
      vart50 = tvar;
    }

    i = i + 1;
  }


  HPDlkT[3] = kTb99;
  HPDukT[3] = kTt99;
  HPDlkT[2] = kTb95;
  HPDukT[2] = kTt95;
  HPDlkT[1] = kTb75;
  HPDukT[1] = kTt75;
  HPDlkT[0] = kTb50;
  HPDukT[0] = kTt50;

  HPDlmu[3] = mub99;
  HPDumu[3] = mut99;
  HPDlmu[2] = mub95;
  HPDumu[2] = mut95;
  HPDlmu[1] = mub75;
  HPDumu[1] = mut75;
  HPDlmu[0] = mub50;
  HPDumu[0] = mut50;

  HPDlvar[3] = varb99;
  HPDuvar[3] = vart99;
  HPDlvar[2] = varb95;
  HPDuvar[2] = vart95;
  HPDlvar[1] = varb75;
  HPDuvar[1] = vart75;
  HPDlvar[0] = varb50;
  HPDuvar[0] = vart50;


  free(stkT);
  free(stmu);
  free(stvar);
  free(nT);
  free(av_yT);


  (*stats).kThist = kThst;
  (*stats).HPDlkT = HPDlkT;
  (*stats).HPDukT = HPDukT;
  (*stats).modekT = modekT;
  (*stats).mukT = mukT;
  (*stats).medkT = medkT;
  (*stats).ga = ga;
  (*stats).gb = gb;

  (*stats).HPDlmu = HPDlmu;
  (*stats).HPDumu = HPDumu;
  (*stats).mumu = mumu;
  (*stats).medmu = medmu;

  (*stats).HPDlvar = HPDlvar;
  (*stats).HPDuvar = HPDuvar;
  (*stats).modevar = modevar;
  (*stats).muvar = muvar;
  (*stats).medvar = medvar;

  (*stats).muT = muT;
  (*stats).sdT = sdT;

  (*stats).pD_null = pD_null;
  (*stats).pD_qtl = pD_qtl;
  (*stats).DIC_null = DIC_null;
  (*stats).DIC_qtl = DIC_qtl;
  (*stats).DIC_diff = DIC_diff;
  (*stats).BIC_null = BIC_null;
  (*stats).BIC_qtl = BIC_qtl;
  (*stats).BF = BF;
  (*stats).logBF = logBF;


  return stats;

}


double qtl_LfocX(double *av_yT, double *avNi, double av_ysq, double av_ybar, double kT, double var, double mu, int ncol, int nrow) {

  // integrated version of qtl_plug, compared to null_plug

  int i;
  double deno,dncol,dnrow,lik,prod1,sum1,sum2,sum3;


  dnrow = (double)nrow;
  nrow = (int)nrow;

  dncol = (double)ncol;
  ncol = (int)ncol;


  sum2 = 0.0;
  prod1 = 0.0;
  for (i=0; i<ncol; i++) {
    if ( avNi[i] > 0.0 ) {
      deno = ( 1.0 - kT + kT*avNi[i] );
      prod1 = prod1 + log(deno);
      sum2 = sum2 + ((avNi[i]*avNi[i]*(av_yT[i] - mu)*(av_yT[i] - mu))/deno);
    }
  }
  prod1 = prod1/2.0;

  sum1 = ( (dnrow*mu*(mu - 2.0*av_ybar)) + av_ysq );
  sum3 = (sum1 - (kT*sum2))/(2.0*var*(1.0 - kT));

  lik = ( (-dnrow/2.0)*log(2.0*M_PI) + ((dncol - dnrow)/2.0)*log(1.0 - kT) - (dnrow/2.0)*log(var) - prod1 - sum3 );

  return lik;

}


double qtl_plug(double *av_yT, double *avNi, double av_ysq, double kT, double var, double mu, double *T, int ncol, int nrow) {

  int i;
  double dnrow,lik,sum2;


  dnrow = (double)nrow;
  nrow = (int)nrow;


  sum2 = 0.0;
  for (i=0; i<ncol; i++) {
    if ( avNi[i] > 0.0 ) {
      sum2 = sum2 + ( (avNi[i])*(mu + T[i])*(mu + T[i] - 2.0*av_yT[i]) );
    }
  }
  sum2 = sum2 + av_ysq;

  lik = ( ((-dnrow/2.0)*log(2.0*M_PI)) - ((dnrow/2.0)*log(1.0 - kT)) - ((dnrow/2.0)*log(var)) - (sum2/(2.0*var*(1.0 - kT))) );

  return lik;

}


double qtl_lik(SIMX *trueX, double *y, double kT, double var, double mu, double *T, int nrow, int nbin) {

  int i,m;
  double Nsamp,lik,sum2;

  Nsamp = 0.0;
  sum2 = 0.0;
  for (i=0; i<nrow; i++) {
    m = ( (*trueX).X[i] - 1 );
    if ( (*trueX).Ni[m] >= nbin ) {
      sum2 = sum2 + (y[i] - mu - T[m])*(y[i] - mu - T[m]);
      Nsamp = Nsamp + 1.0;
    }
  }

  lik = ( ((-Nsamp/2.0)*log(2.0*M_PI)) - ((Nsamp/2.0)*log(1.0 - kT)) - ((Nsamp/2.0)*log(var)) - (sum2/(2.0*var*(1.0 - kT))) );


  return lik;

}


double qtl_Lfoc(SIMX *trueX, double *y, double kT, double var, double mu, int nrow, int ncol, int nbin) {

  int i,m;
  double deno,dn,Nsamp,lik,prod1,sum1,sum2,sum3,tsamp;

  double *yTbar=NULL;

  yTbar = (double*)calloc(ncol,sizeof(double));

  Nsamp = 0.0;
  sum1 = 0.0;
  for (i=0; i<nrow; i++) {
    m = ( (*trueX).X[i] - 1 );
    if ( (*trueX).Ni[m] >= nbin ) {
      sum1 = sum1 + ((y[i] - mu)*(y[i] - mu));
      yTbar[m] = yTbar[m] + y[i];
      Nsamp = Nsamp + 1.0;
    }
  }

  tsamp = 0.0;
  prod1 = 0.0;
  sum2 = 0.0;
  for (i=0; i<ncol; i++) {
    if ( (*trueX).Ni[i] >= nbin ) {
      dn = (double)(*trueX).Ni[i];
      deno = ( 1.0 - kT + kT*dn );
      tsamp = tsamp + 1.0;

      prod1 = prod1 + log(deno);
      yTbar[i] = yTbar[i]/dn;

      sum2 = sum2 + ((dn*dn*(yTbar[i] - mu)*(yTbar[i] - mu))/deno);
    }
  }
  prod1 = prod1/2.0;

  sum3 = (sum1 - (kT*sum2))/(2.0*var*(1.0 - kT));

  lik = ( (-Nsamp/2.0)*log(2.0*M_PI) + ((tsamp - Nsamp)/2.0)*log(1.0 - kT) - (Nsamp/2.0)*log(var) - prod1 - sum3 );

  free(yTbar);

  return lik;
}


double null_plug(double *avNi, double av_ysq, double av_ybar, double var, double mu, int ncol, int nrow) {

  int i;
  double dnrow,lik,sum1;

  dnrow = (double)nrow;
  nrow = (int)nrow;

  sum1 = ( (dnrow*mu*(mu - 2.0*av_ybar)) + av_ysq );

  lik = ( ((-dnrow/2.0)*log(2.0*M_PI)) - ((dnrow/2.0)*log(var)) - (sum1/(2.0*var)) );


  return lik;

}


double null_lik(SIMX *trueX, double *y, double var, double mu, int nrow, int nbin) {

  int i,j,m;
  double lik,sum1,Nsamp;

  Nsamp = 0.0;
  sum1 = 0.0;
  for (i=0; i<nrow; i++) {
    m = ( (*trueX).X[i] - 1 );
    if ( (*trueX).Ni[m] >= nbin ) {
      sum1 = sum1 + (y[i] - mu)*(y[i] - mu);
      Nsamp = Nsamp + 1.0;
    }
  }

  lik = ( ((-Nsamp/2.0)*log(2.0*M_PI)) - ((Nsamp/2.0)*log(var)) - (sum1/(2.0*var)) );


  return lik;

}


SIMX* drawX(XMAT *xmat, int ncol, int nrow, long *idum) {

  int i,j,m;
  double check,dncol,max,muNi,ran,sdNi,varNi;

  int *Xvec=NULL, *Ni=NULL;
  double *Xprob=NULL;

  SIMX *trueX=NULL;


  trueX = (SIMX*)calloc(1,sizeof(SIMX));

  Xvec = (int*)calloc(nrow,sizeof(int));
  Xprob = (double*)calloc(nrow,sizeof(double));

  Ni = (int*)calloc(ncol,sizeof(int));


  dncol = (double)ncol;
  ncol = (int)ncol;


  for (i=0; i<nrow; i++) {

    ran = (double)ran2(idum);
    m = 1;
    while ( ran > (*xmat).cumX[i][m-1] ) {
      m++;
    }

    if (m > ncol) {
      Rprintf("m = %i, ran = %f\n",m,ran);
      exit(1);
    }

    Xvec[i] = m;
    Xprob[i] = (*xmat).X[i][m-1];

    Ni[Xvec[i] - 1] = Ni[Xvec[i] - 1] + 1;
  }


  muNi = 0.0;
  for (i=0; i<ncol; i++) {
    muNi = muNi + (double)Ni[i];
  }
  muNi = muNi/dncol;

  varNi = 0.0;
  for (i=0; i<ncol; i++) {
    varNi = varNi + ((double)Ni[i] - muNi)*((double)Ni[i] - muNi);
  }
  varNi = varNi/dncol;

  sdNi = sqrt(varNi);


  free(Xprob);


  (*trueX).X = Xvec;
  (*trueX).Ni = Ni;
  (*trueX).sdNi = sdNi;

  return trueX;

}


GRKT* truegridkT(SIMX *trueX, double *y, int ncol, int nrow, int nbin) {

  int i,j,m,n,prc;

  double dprc,inc,inter,kT,rhosq,pow1,pow2,stmin,stmax;
  double deno,dn,Ngen,Nsamp,prod1,sum1,sum2,sum3,sum4,sum5,sum6;

  double *kT_pdf=NULL, *kT_cdf=NULL, *yTbar=NULL;

  GRKT *grid=NULL;


  prc = 200;
  dprc = 200.0;


  grid = (GRKT*)calloc(1,sizeof(GRKT));

  kT_pdf = (double*)calloc((prc+1),sizeof(double));
  kT_cdf = (double*)calloc((prc+1),sizeof(double));

  yTbar = (double*)calloc(ncol,sizeof(double));



  sum1 = 0.0;
  sum6 = 0.0;
  Nsamp = 0.0;
  for (i=0; i<nrow; i++) {
    m = ( (int)(*trueX).X[i] - 1 );
    if ( (int)(*trueX).Ni[m] >= nbin ) {
      yTbar[m] = yTbar[m] + y[i];
      sum1 = sum1 + (y[i]*y[i]);
      sum6 = sum6 + y[i];
      Nsamp = Nsamp + 1.0;
    }
  }
  sum6 = sum6/Nsamp;

  Ngen = 0.0;
  for (j=0; j<ncol; j++) {
    if ( (int)(*trueX).Ni[j] >= nbin ) {
      yTbar[j] = yTbar[j]/((double)(*trueX).Ni[j]);
      Ngen = Ngen + 1.0;
    }
  }


  kT = 0.0;
  inc = 1.0/dprc;
  stmax = -1000000.0;

  for (i=0; i<prc; i++) {

    sum2 = 0.0;
    prod1 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    for (j=0; j<ncol; j++) {
      if ( (int)(*trueX).Ni[j] >= nbin ) {

	dn = (double)(*trueX).Ni[j];
	deno = ( 1.0 - kT + (kT*dn) );

        sum2 = sum2 + (dn/deno);
	prod1 = prod1 + log(deno);
	sum3 = sum3 + (dn*dn*(yTbar[j])*(yTbar[j]))/deno;
	sum4 = sum4 + (dn*(yTbar[j]))/deno;
      }
    }


    inter = ( sum1 - (kT*sum3) - ((1.0-kT)*(sum4*sum4/sum2)) );


    sum2 = -log(sum2)/2.0;
    prod1 = -prod1/2.0;

    pow1 = ( Ngen - 1.0 )/2.0;
    pow2 = ( Nsamp - 1.0 )/2.0;

    kT_pdf[i] = ( (pow1*log(1-kT)) + sum2 + prod1 - (pow2*log(inter)) );

    if (kT_pdf[i] > stmax) {
      stmax = kT_pdf[i];
    }

    kT = kT + inc;

  }

  sum5 = 0.0;
  stmin = (stmax - 703.0);

  for (i=0; i<prc; i++) {
    kT_pdf[i] = (kT_pdf[i] - stmin);
    kT_pdf[i] = exp(kT_pdf[i]);
    sum5 = sum5 + kT_pdf[i];
  }

  kT_pdf[0] = kT_pdf[0]/sum5;
  kT_cdf[0] = kT_pdf[0];
  for (i=1; i<=prc; i++) {
    kT_pdf[i] = kT_pdf[i]/sum5;
    kT_cdf[i] = kT_cdf[i-1] + kT_pdf[i];
  }


  (*grid).pdf = kT_pdf;
  (*grid).cdf = kT_cdf;
  (*grid).yTbar = yTbar;
  (*grid).sum_ysq = sum1;
  (*grid).sum_y = sum6;
  (*grid).N = Nsamp;
  (*grid).t = Ngen;

  return grid;

}


double drawkT(GRKT *kTdist,long *idum) {

  int i,prc;
  double dprc,kT,ran;


  prc = 200;
  dprc = 200.0;


  ran=(double)ran2(idum);
  i=0;
  while ( ran > (*kTdist).cdf[i] ) {
    i++;
  }
  if ( i > (prc + 1) ) {
    Rprintf("error in draw of kT\n");
    exit(1);
  }

  kT = ((double)i)/dprc;


  return kT;

}


double draw_knownvar(GRKT *kTdist, int *Ni, int ncol, double kT, double nu, int nbin) {

  int j;

  double chi,scale,var;
  double deno,dn,sum2,sum3,sum4;


  if ( kT == 1.0 ) {

    var = 0.0;

  } else {

    chi = rchisq(nu);

    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    for (j=0; j<ncol; j++) {
      if ( (int)Ni[j] >= nbin ) {

	dn = (double)Ni[j];
	deno = ( 1.0 - kT + kT*dn );

        sum2 = sum2 + (dn/deno);
	sum3 = sum3 + (dn*dn*((*kTdist).yTbar[j])*((*kTdist).yTbar[j]))/deno;
	sum4 = sum4 + (dn*((*kTdist).yTbar[j]))/deno;
      }
    }


    // scale = rho squared multiplied by nu=degrees of freedom
    scale = ( ((*kTdist).sum_ysq)/(1.0 - kT) - (kT/(1.0 - kT))*sum3 - (sum4*sum4/sum2) );
    var = scale/chi;

  }

  return var;

}


double draw_nullvar(SIMX *selX, double *y, int nrow, int nbin) {

  int j,m;

  double chi,nu,scale,var;
  double dnrow,sum2,sum3,Nsamp;


  dnrow = (double)nrow;
  nrow = (int)nrow;

  Nsamp = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  for (j=0; j<nrow; j++) {
    m = ( (int)(*selX).X[j] - 1 );
    if ( (int)(*selX).Ni[m] >= nbin ) {
      sum2 = sum2 + (y[j]*y[j]);
      sum3 = sum3 + y[j];
      Nsamp = Nsamp + 1.0;
    }
  }
  sum3 = sum3/Nsamp;


  nu = (Nsamp - 1.0);
  chi = rchisq(nu);

  // scale = rho squared multiplied by nu=degrees of freedom
  scale = ( sum2 - (Nsamp*sum3*sum3) );
  var = scale/chi;


  return var;

}


double draw_knownmu(GRKT *kTdist, int *Ni, int ncol, double kT, double var, int nbin) {

  int j;

  double mu,nmean,nsigmasq,sigma;
  double deno,dn,sum2,sum4;


  sum2 = 0.0;
  sum4 = 0.0;
  for (j=0; j<ncol; j++) {
    if ( (int)Ni[j] >= nbin ) {

      dn = (double)Ni[j];
      deno = ( 1.0 - kT + kT*dn );

      sum2 = sum2 + (dn/deno);
      sum4 = sum4 + (dn*((*kTdist).yTbar[j]))/deno;
    }
  }

  nmean = sum4/sum2;
  nsigmasq = var/sum2;

  sigma = sqrt(nsigmasq);
  mu = rnorm(0.0, sigma);
  mu = mu + nmean;

  return mu;

}


double draw_nullmu(SIMX *selX, double *y, int nrow, double var, int nbin) {

  int j,m;

  double mu,nmean,nsigmasq,sigma;
  double dnrow,sum2,Nsamp;


  dnrow = (double)nrow;
  nrow = (int)nrow;


  sum2 = 0.0;
  Nsamp = 0.0;
  for (j=0; j<nrow; j++) {
    m = ( (int)(*selX).X[j] - 1 );
    if ( (int)(*selX).Ni[m] >= nbin ) {
      sum2 = sum2 + y[j];
      Nsamp = Nsamp + 1.0;
    }
  }
  sum2 = sum2/Nsamp;


  nmean = sum2;
  nsigmasq = var/Nsamp;

  sigma = sqrt(nsigmasq);
  mu = rnorm(0.0, sigma);
  mu = mu + nmean;

  return mu;

}


double draw_knownTi(GRKT *kTdist, int *Ni, double kT, double var, double mu, int nbin, int index) {

  double dn,deno;
  double nmean,nsigmasq,sigma,ti;


  ti = 0.0;
  if ( (int)Ni[index] >= nbin ) {

    dn = (double)Ni[index];
    deno = ( 1.0 - kT + kT*dn );

    nmean = (kT*dn*((*kTdist).yTbar[index] - mu))/deno;
    nsigmasq = (kT*(1.0 - kT)*var)/deno;

    sigma = sqrt(nsigmasq);
    ti = rnorm(0.0, sigma);
    ti = ti + nmean;

  }

  return ti;

}

