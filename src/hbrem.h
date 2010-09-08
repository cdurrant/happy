

typedef struct {
  int *X;           // genotype labels start from 1, not 0
  int *Ni;
  double sdNi;
} SIMX;

typedef struct {
  double *pdf, *cdf, *yTbar;
  double sum_ysq, sum_y, N, t;
} GRKT;

typedef struct {
  double *yTbar;
  double *kT, *var, *mu;
  double *null_var, *null_mu, *lik_null, *lik_qtl;
  double **T, **yTbarX;
  int **NiX;
  double t, N;
} POST;

typedef struct {
  int *kThist;
  double *HPDlkT, *HPDukT;
  double mukT, medkT, modekT, ga, gb;
  double *HPDlmu, *HPDumu;
  double mumu, medmu;
  double *HPDlvar, *HPDuvar;
  double muvar, medvar, modevar;
  double *muT, *sdT;
  double pD_null, pD_qtl, DIC_null, DIC_qtl, DIC_diff;
  double BIC_null, BIC_qtl, BF, logBF;
} SSTA;



SIMX* drawX(XMAT *xmat, int ncol, int nrow, long *idum);

GRKT* truegridkT(SIMX *selX, double *y, int ncol, int nrow, int nbin);
double drawkT(GRKT *kTdist,long *idum);

double draw_knownvar(GRKT *kTdist, int *Ni, int ncol,  double kT, double nu, int nbin);
double draw_knownmu(GRKT *kTdist, int *Ni, int ncol,  double kT, double var, int nbin);
double draw_knownTi(GRKT *kTdist, int *Ni, double kT, double var, double mu, int nbin, int index);

double draw_nullvar(SIMX *selX, double *y, int nrow, int nbin);
double draw_nullmu(SIMX *selX, double *y, int nrow, double var, int nbin);

double qtl_lik(SIMX *trueX, double *y, double kT, double var, double mu, double *T, int nrow, int nbin);
double null_lik(SIMX *trueX, double *y, double var, double mu, int nrow, int nbin);

double qtl_plug(double *av_yT, double *avNi, double av_ysq, double kT, double var, double mu, double *T, int ncol, int nrow);
double null_plug(double *avNi, double av_ysq, double av_ybar, double var, double mu, int ncol, int nrow);

double qtl_LfocX(double *av_yT, double *avNi, double av_ysq, double av_ybar, double kT, double var, double mu, int ncol, int nrow);
double qtl_Lfoc(SIMX *trueX, double *y, double kT, double var, double mu, int nrow, int ncol, int nbin);

POST* single_locus_jointpost(SIMX *trueX,double *y,int nsim,int ncol,int nrow,int nbin,long *idum);
POST* single_locus_jointpostX(XMAT *xmat,double *y,int nsim,int ncol,int nrow,int nbin,long *idum);
SSTA* single_locus_sumstats(POST *SLsamp, SIMX *trueX, double *y, int nsim, int ncol, int nrow, int nbin);
SSTA* single_locus_sumstatsX(XMAT *xmat, POST *SLsamp, double *y, int nsim, int ncol, int nrow, int nbin);
