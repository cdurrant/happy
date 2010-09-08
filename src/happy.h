#ifndef _HAPPY_H_ 
#define _HAPPY_H_

#include<stdio.h>

typedef enum { BOTH, LEFT, UNLINKED, RIGHT, LINKAGE_STATES } LINKAGESTATES;

typedef enum {UNKNOWN, MALE, FEMALE } GENDER;

#define NOT_FOUND  -1

#define ND_ALLELE "NA"


typedef struct {
  int markers;
  int *chrom1;
  int *chrom2;
} CHROM_PAIR;

typedef struct {
  int alleles;
  char *marker_name;
  char **allele_name;
  double *allele_freq; /* observed frequency of alleles */
  double **pr_AtoS; /* prob of strain s | allele a */
  double entropy;
  char chromosome[20];         /* the chromosome of the marker */
  double position;             /* estimate of the distance of the QTL from the left-hand end */
  double ProbSame;             /* prob of an observable recomb between this and the next marker */
  double **prior;              /* array of probabilities that the pair of QTL states are in */
  int *which_allele;
  int ND;                      /* index of the allele corresponding to "ND" (missing) */
} ALLELE_FREQ;

typedef struct {
  int strains;
  char **strain_name;
  int markers;
  int generations;
  double *Pr_ss;               /* array of transition probabilities for staying in same state */
  double *Pr_st;               /* array of transition probabilities for changing state */
  ALLELE_FREQ *af;
  double MinDist;               /* minimum distance apart for markers */
} ALLELES;

/* DP_MATRICES contains the forward and backward dynamic-programming matrices */
typedef struct {
  double ***Left;              /* forward DP matrix */
  double ***Right;             /* backward DP matrix */
  double *NonRecomb;            /* expected number of non-recombinant chromosomes */
} DP_MATRICES;


/* HAPLOID_DP_MATRICES contains the forward and backward dynamic-programming matrices for haploid genomes*/
typedef struct {
  double **Left;              /* forward DP matrix */
  double **Right;             /* backward DP matrix */
  double *NonRecomb;            /* expected number of non-recombinant chromosomes */
} HAPLOID_DP_MATRICES;

/* QTL_FIT contains all the data associated with fitting the QTL */

typedef struct {
  int locus;                   /* left-hand marker locus of current interest */
  double rss;                  /* the residual sum of squares */
  double fss;                  /* the fitting sum of squares */
  double F;                    /* the F-ratio */
  double pval;                 /* hte p-value of the F */
  double mean;                 /* the estimated mean trait response */
  double *trait;               /* array of estimated trait effects for each strain */
  double *trait_error;         /* array of estimated trait standard errors for each strain */
  double sigma;                /* estimated residual standard error */
  int *trait1;                 /* array of N predicted trait states for chrom1 */
  int *trait2;                 /* array of N predicted trait states for chrom2 */
  int left, right;             /* counts of the number of chromosomes allocated left and right */
  double **design_matrix;      /* alternative expectations of traits for regression */
  double *residual;            /* residduals from fit */
} QTL_FIT;

/* ANCESTRY is a matrix giving the expected proportion of each founder strain to be found in each subject; only used if these fractions are unequal */

typedef struct {
  int N;
  int S;
  char **strain_name;
  char **subject_name;
  double **prob;
  double ****pr_AtoS; /* subject-specific prob of strain s | allele a */
}  ANCESTRY;

/* QTL_DATA  is a portmanteau struct that contains pretty much all the data */

typedef struct {
  char *filename;                  /* Name of the data-set */
  int N;                       /* Number of individuals */
  int M;                       /* Number of markers */
  int S;                       /* Number of strains */
  char *missingCode;           /* missing allele code */
  int haploid;                  /* boolean indicating if data are haploid (== inbred)*/
  ALLELES *alleles;            /* pointer to ALLELES struct containing the founder strain info */
  ANCESTRY *an;               /* pointer to ANCESTRY struct containing the subject-specific ancestral strain probabilities (can be null)*/
  char **name;                 /* array of names of individuals */
  double *observed;            /* array of N observed trait values */
  CHROM_PAIR *genos;           /* array of N CHROM_PAIR structs containing raw marker genotypes */
  CHROM_PAIR *haplos;          /* array of N CHROM_PAIR structs containing deduced strain haplotypes */
  DP_MATRICES *dp_matrices;    /* array of N dynamic-programming matrices for computing the priors */
  HAPLOID_DP_MATRICES *haploid_dp_matrices;    /* array of N dynamic-programming matrices for computing the priors for haploid genomes*/
  QTL_FIT *fit;                /* array of QTL_FIT structures, for each marker locus */
  QTL_FIT *null_model;         /* QTL_FIT struct containing the null model fit */
  double drop;                 /* factor for reducing search space of prior configurations */
  int from_marker;
  int to_marker;
  int phase_known;             /* switch indicating whether the phase of the genotypes is known - ie they are haplotypes */
  int use_parents;            /* switch indicating whether pedigree data is available */
  int *mother;                 /* array of indices to mother , -1 if not present */
  int *father;                 /* array of indices to father, -1 if not present */
  char **family;                /* array of family names (null if no family) */
  int *sex;                    /* array of sex indices +1 (male) 0 (unknown) 2 (female) */
} QTL_DATA;



typedef struct  {
  double prior, posterior, cum;
} QTL_PRIOR;


typedef struct {
  double key, value;
}
KV;

typedef struct {
  char *key;
  int id;;
}
PARENT_KEY;


/* function prototypes */

CHROM_PAIR *new_chrom_pair( int markers );
ALLELES *input_allele_frequencies( FILE *fp, int generations, char *missingCode, double MinDist, int verbose );
QTL_PRIOR ***compute_qtl_priors(  QTL_DATA *qtl, QTL_PRIOR ***qp, int locus, double **prior  );
QTL_PRIOR **compute_haploid_qtl_priors(  QTL_DATA *qtl, QTL_PRIOR **qp, int locus  );
int qpcmp( const void *A, const void *B );
double fit_null_qtl_model( QTL_DATA *qtl_data );
void allocate_traits( QTL_DATA *q, QTL_PRIOR ***p, QTL_FIT *fit, int mode );
void fit_qtl( QTL_DATA *q, int locus, int verbose, int shuffles );
double fit_linear_additive_model( QTL_DATA *qtl, QTL_FIT *fit, int shuffles, int verbose );
QTL_FIT *allocate_qtl_fit( QTL_FIT *fit, int N, int strains );
void print_qtl_data ( QTL_DATA *q, QTL_FIT *fit, FILE *fp );
void qtl_fit_cp( QTL_FIT *fit1, QTL_FIT *fit2, int N, int S );
QTL_DATA *read_qtl_data( FILE *fp, char *name, ALLELES *a, int verbose, int use_parents, int ped_format, char *missing );
void write_qtl_data( FILE *fp, QTL_DATA *q );
int check_and_apply_ancestry(QTL_DATA *q );
double ***summed_dp_matrix( QTL_DATA *qtl, int individual,  double *p1, double *p2, int direction );
double **haploid_summed_dp_matrix( QTL_DATA *qtl, int individual, double *Pr_ss, double *Pr_st, int direction );
void create_summed_dp_matrices( QTL_DATA *q );
void create_haploid_summed_dp_matrices( QTL_DATA *q );
QTL_PRIOR ***allocate_qtl_priors( QTL_DATA *q );
QTL_PRIOR **allocate_haploid_qtl_priors( QTL_DATA *qtl );
int remove_partial_fit( QTL_DATA *q, char *marker, int verbose, int fail );
void permute_data( double *data , int N );
void permute_genotypes( QTL_DATA *q );
void pointwise_mapping( QTL_DATA *q, double step, int verbose );
void pointwise_interval_mapping_probabilities( QTL_DATA *q, int locus, double c, double **prior );
QTL_DATA *resample_qtl_data( QTL_DATA *q, QTL_DATA *r );
void bootstrap_analysis( QTL_DATA *q, int bootstrap, char *bootstart, char *bootstop, int verbose );
void sequential_fit( QTL_DATA *q );
int marker_index( const char *name, QTL_DATA *q, const int isIntervalModel );
int genotype_difference( QTL_DATA *q, int i, int j );
int pdump_prob_data( FILE *fp, int locus, QTL_DATA *q );
double ** additive_design_matrix( QTL_DATA *q, int locus );
SEXP happy( SEXP datafile, SEXP allelesfile, SEXP generations, SEXP phase, SEXP file_format, SEXP missing_code, SEXP do_dp, SEXP min_dist, SEXP haploid, SEXP anfilename );
SEXP getListElement(SEXP list, char *str);
QTL_DATA * validateParams ( SEXP handle, SEXP marker, int *locus, const int isIntervalModel);
SEXP happyprobs ( SEXP handle, SEXP marker );
SEXP happydesign( SEXP handle, SEXP marker, SEXP model );
double phaseProb( int a1, int a2, int m1, int m2, int p1, int p2, int NA );
SEXP happygenotype ( SEXP handle, SEXP marker );
SEXP happynonrecomb ( SEXP handle, SEXP marker );

ANCESTRY *read_subject_ancestries( FILE *fp, char *filename, int verbose );
void heterozygosity( QTL_DATA *q );
double subject_heterozygosity( QTL_DATA *q, int individual );
double marker_heterozygosity( QTL_DATA *q, int marker );


#endif
