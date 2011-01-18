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

#include"happy.h"
#include"readline.h"
#include"cmp.h"
#include"stats.h"

static QTL_DATA *qtldata[100];
static int nqtldata = 0;
int entrycmp( const void *a, const void *b );


SEXP happy( SEXP datafile, SEXP allelesfile, SEXP generations, SEXP phase, SEXP file_format, SEXP missing_code, SEXP do_dp, SEXP min_dist, SEXP haploid, SEXP ancestryfile ) {
  QTL_DATA *q = NULL;
  ALLELES *a = NULL;
  FILE *dfp=NULL, *afp=NULL, *anfp=NULL;
  const char *afilename, *anfilename;
  const char *dfilename;
  int gen = 0;
  double g;
  int verbose = 0;
  int use_parents = 0;
  int ped_format = 0;
  int phaseKnown = 0;
  SEXP strains,markers, chromosome, subjects, phenotypes, map, family;
  SEXP ans = R_NilValue;
  SEXP names, handle, class;
  char *MissingCode;
  int i;
  const char *PhaseStr;
  const char *File_FormatStr;
  int Do_dp, Haploid;
  double MinDist = 1.0e-5;

  if ( ! isString(datafile) || length(datafile) != 1 )
    error( "datafile is not a string");
  dfilename = CHAR(STRING_ELT(datafile,0));
  
  if ( ! (dfp = fopen( dfilename, "r" ) ) )
    error( "could not open data file" );

  if ( ! isString(allelesfile) || length(allelesfile) != 1 ) 
    error( "allelesfile is not a string");
  afilename = CHAR(STRING_ELT(allelesfile,0));

  if ( ! (afp = fopen( afilename, "r" ) ) )
    error( "could not open alleles file" );

  if ( isString(ancestryfile) && length(ancestryfile) == 1 ) {
    anfilename = CHAR(STRING_ELT(ancestryfile,0));
    if ( ! (anfp = fopen( anfilename, "r" ) ) )
      error( "could not open ancestry file" );
  }

  if ( ! isNumeric(generations) || length(generations) != 1 )
    error( "generations is not numeric");
  g = REAL(generations)[0];
  gen = (int)g; 

  if ( ! isString(phase) || length(phase) != 1 ) 
    error( "phase is not a string");
  PhaseStr = CHAR(STRING_ELT(phase,0));

  if ( ! isString(file_format) || length(file_format) != 1 ) 
    error( "file_format is not character(1)");
  File_FormatStr = CHAR(STRING_ELT(file_format,0));

  if ( ! isString(missing_code) || length(missing_code) != 1 )
    error( "missing_code is not character(1)");
  if ( strlen( CHAR(STRING_ELT(missing_code,0)) ) > 0 ) {
    MissingCode = (char*)CHAR(STRING_ELT(missing_code,0));
  } else {
    MissingCode = strdup(ND_ALLELE);
  }

  if ( ! isNumeric(do_dp) || length(do_dp) != 1 )
    error( "do_dp is not numeric(1)");
  Do_dp = INTEGER(do_dp)[0];

  if ( ! isNumeric(haploid) || length(haploid) != 1 )
    error( "haploid is not numeric(1)");
  Haploid = INTEGER(haploid)[0]; 

  if ( ! isNumeric(min_dist) || length(min_dist) != 1 )
    error( "min_dist is not numeric(1)");
  else if ( isNumeric(min_dist) )
    MinDist = (double)REAL(min_dist)[0];


  Rprintf( "mindist: %g\n", MinDist );
  Rprintf( "datafile %s allelesfile %s gen %d\n", dfilename, afilename, gen );
  Rprintf( "genotype phase: %s\n", PhaseStr);

  if ( ! strcmp( File_FormatStr, "ped") ) 
    ped_format = 1;
  else
    ped_format = 0;

  if ( ! strcmp( PhaseStr, "unknown" ) ) 
    use_parents = 0;
  else if ( ! strcmp( PhaseStr, "estimate" ) ) {
    use_parents = 1;
    ped_format = 1;
  }
  else if ( ! strcmp( PhaseStr, "known" ) ) {
    use_parents = 0;
    phaseKnown = 1;
  }

  if ( use_parents ) Rprintf( "using parental genotypes to help determine phase\n");

  a = input_allele_frequencies( afp, gen, MissingCode, MinDist, verbose );
  Rprintf( "a->markers %d\n", a->markers );
  q = read_qtl_data( dfp, (char*)dfilename, a,  verbose, use_parents, ped_format, MissingCode );
  q->an = read_subject_ancestries( anfp, (char*)anfilename, verbose );
  q->phase_known = phaseKnown;
  q->haploid = Haploid;
  if ( Haploid ) Rprintf( "assuming haploid(inbred) genotypes\n");
  if ( q->an ) 
    check_and_apply_ancestry( q );

  Rprintf( "dfile %s afile %s gen %d\n", dfilename, afilename, gen );

  
  if ( Do_dp ) {
    if ( q->haploid ) {
      /*      heterozygosity(q ); */
      create_haploid_summed_dp_matrices( q );
    }
    else
      create_summed_dp_matrices( q );
  }

  PROTECT(strains=allocVector(STRSXP,q->S));
  for(i=0;i<q->S;i++) {
    SET_STRING_ELT(strains,i, mkChar(a->strain_name[i]));
  }

  PROTECT(markers=allocVector(STRSXP,q->M));
  for(i=0;i<q->M;i++) {
    SET_STRING_ELT(markers,i, mkChar(a->af[i].marker_name));
  }

  PROTECT(chromosome=allocVector(STRSXP,q->M));
  for(i=0;i<q->M;i++) {
    SET_STRING_ELT(chromosome,i, mkChar(a->af[i].chromosome));
  }

  PROTECT(map=allocVector(REALSXP,q->M));
  for(i=0;i<q->M;i++) {
    REAL(map)[i] = a->af[i].position;
  }

  PROTECT(subjects=allocVector(STRSXP,q->N));
  for(i=0;i<q->N;i++) {
    SET_STRING_ELT(subjects,i, mkChar(q->name[i]));
  }

  PROTECT(family=allocVector(STRSXP,q->N));
  for(i=0;i<q->N;i++) {
  if ( q->family && q->family[i] )
    SET_STRING_ELT(family,i, mkChar(q->family[i]));
  else
    SET_STRING_ELT(family,i, mkChar(""));

  }

  PROTECT(phenotypes=allocVector(REALSXP,q->N));
  for(i=0;i<q->N;i++) {
    REAL(phenotypes)[i] = q->observed[i];
  }

  PROTECT( ans = allocVector( VECSXP, 8 ) );
  PROTECT(names=allocVector(STRSXP, 8 ) );

  PROTECT(handle = allocVector(INTSXP,1));
  qtldata[nqtldata] = q;
  INTEGER(handle)[0] = nqtldata++;

  
  /*  SET_VECTOR_ELT(ans,0,strains);
  SET_VECTOR_ELT(names,0,mkChar("strains"));
  SET_VECTOR_ELT(ans,1,markers);
  SET_VECTOR_ELT(names,1,mkChar("markers"));
  SET_VECTOR_ELT(ans,2,map);
  SET_VECTOR_ELT(names,2,mkChar("map"));
  SET_VECTOR_ELT(ans,3,subjects);
  SET_VECTOR_ELT(names,3,mkChar("subjects"));
  SET_VECTOR_ELT(ans,4,phenotypes);
  SET_VECTOR_ELT(names,4,mkChar("phenotypes"));
  SET_VECTOR_ELT(ans,5,handle);
  SET_VECTOR_ELT(names,5,mkChar("handle"));
  SET_VECTOR_ELT(ans,6,chromosome);
  SET_VECTOR_ELT(names,6,mkChar("chromosome"));
  SET_VECTOR_ELT(ans,7,family);
  SET_VECTOR_ELT(names,7,mkChar("family"));
  */
  SET_VECTOR_ELT(ans,0,strains);
  SET_STRING_ELT(names,0,mkChar("strains"));
  SET_VECTOR_ELT(ans,1,markers);
  SET_STRING_ELT(names,1,mkChar("markers"));
  SET_VECTOR_ELT(ans,2,map);
  SET_STRING_ELT(names,2,mkChar("map"));
  SET_VECTOR_ELT(ans,3,subjects);
  SET_STRING_ELT(names,3,mkChar("subjects"));
  SET_VECTOR_ELT(ans,4,phenotypes);
  SET_STRING_ELT(names,4,mkChar("phenotypes"));
  SET_VECTOR_ELT(ans,5,handle);
  SET_STRING_ELT(names,5,mkChar("handle"));
  SET_VECTOR_ELT(ans,6,chromosome);
  SET_STRING_ELT(names,6,mkChar("chromosome"));
  SET_VECTOR_ELT(ans,7,family);
  SET_STRING_ELT(names,7,mkChar("family"));

  setAttrib( ans, R_NamesSymbol, names );

  UNPROTECT(10);

  /* set the class */

  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("happy"));
  classgets(ans, class);
  UNPROTECT(1);


  return ans;
}


SEXP getListElement(SEXP list, char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}


/* copy a QTL_FIT structure */
void qtl_fit_cp( QTL_FIT *fit1, QTL_FIT *fit2, int N, int S ) {

  int i;

  fit1->locus = fit2->locus;
  fit1->rss = fit2->rss;
  fit1->fss = fit2->fss;
  fit1->F = fit2->F;
  fit1->pval = fit2->pval;
  fit1->mean = fit2->mean;
  fit1->sigma = fit2->sigma;
  fit1->left = fit2->left;
  fit1->right = fit2->right;

  for(i=0;i<N;i++) {
    fit1->trait1[i] = fit2->trait1[i];
    fit1->trait2[i] = fit2->trait2[i];
  }

  for(i=0;i<S;i++) {
    fit1->trait[i] = fit2->trait[i];
    fit1->trait_error[i] = fit2->trait_error[i];
  }
}


QTL_FIT *allocate_qtl_fit( QTL_FIT *fit, int N, int strains ) {

  if ( fit == NULL ) 
    fit = (QTL_FIT*)calloc(1,sizeof(QTL_FIT));

  fit->trait = (double*)calloc(strains,sizeof(double));
  fit->trait_error = (double*)calloc(strains,sizeof(double));
  fit->trait1 = (int*)calloc(N,sizeof(int));
  fit->trait2 = (int*)calloc(N,sizeof(int));

  return fit;
}


QTL_DATA *read_qtl_data( FILE *fp, char *name, ALLELES *a,  int verbose, int use_parents, int ped_format, char *missingCode ) {

  QTL_DATA *q = (QTL_DATA*)calloc(1,sizeof(QTL_DATA));
  int max_N = 10000;
  int m;
  int bufsize = 100+a->markers*20;
  char *buffer = (char*)calloc(bufsize,sizeof(char));
  char **mother = NULL;
  char **father= NULL;
  double NaN = nan("char-sequence");
  int *pcount;
  int nparents=0;
  q->alleles = a;
  strncpy(q->filename, name, MAX_LENGTH_FILENAME);
  q->N = 0;
  q->M = a->markers;
  q->S = a->strains;
  q->observed = (double*)calloc(max_N,sizeof(double));
  q->genos = (CHROM_PAIR*)calloc(max_N,sizeof(CHROM_PAIR));
  q->name = (char**)calloc(max_N,sizeof(char*));
  q->family = (char**)calloc(max_N,sizeof(char*));
  q->use_parents = use_parents;
  strncpy(q->missingCode,missingCode,MAX_LENGTH_MISSINGCODE);

  if ( use_parents || ped_format ) {
    q->sex = (int*)calloc(max_N, sizeof(int));
    mother = (char**)calloc(max_N, sizeof(char*));
    father = (char**)calloc(max_N, sizeof(char*));
    Rprintf( "Reading phenotype and genotype data from ped file %s\n", name );
  }
  else
    Rprintf( "Reading phenotype and genotype data from data file %s\n", name );
  while ( skip_comments( fp, buffer ) != EOF ) {
    char *str1, *str2;
    m = 0;
    int ok = 0;
    if ( q->N >= max_N ) {
      max_N *= 2;
      q->observed = (double*)realloc(q->observed,max_N*sizeof(double));
      q->genos = (CHROM_PAIR*)realloc(q->genos,max_N*sizeof(CHROM_PAIR));
      q->name = (char**)realloc(q->name,max_N*sizeof(char*));
      q->family = (char**)realloc(q->family,max_N*sizeof(char*));
      if ( use_parents ) {
	q->sex = (int*)realloc(q->sex, max_N*sizeof(int));
	mother = (char**)realloc(mother, max_N*sizeof(char*));
	father = (char**)realloc(father, max_N*sizeof(char*));
      }
    }
    q->genos[q->N].markers = q->M;
    q->genos[q->N].chrom1 = (int*)calloc(q->M,sizeof(int));
    q->genos[q->N].chrom2 = (int*)calloc(q->M,sizeof(int));

    if ( use_parents || ped_format ) {
      char *family = strtok( buffer, " 	" );
      char *id = strtok( NULL, " 	" );
      char *dad = strtok( NULL,  " 	" );
      char *mum = strtok( NULL,  " 	" );
      char *sex = strtok( NULL,  " 	" );
      char *pheno = strtok( NULL,  " 	" );
      

      if ( family && id && dad && mum && sex && pheno ) {
	char *endptr;
	q->family[q->N] = (char*)strdup(family);
	father[q->N] = (char*)strdup(dad);
	mother[q->N] = (char*)strdup(mum);
	q->sex[q->N] = atoi(sex);
	q->name[q->N] = (char*)strdup(id);
	/*	Rprintf( "name %s\n", q->name[q->N]); */
	q->observed[q->N] = strtod(pheno,&endptr);
	if (endptr == pheno)
	  q->observed[q->N] = NaN;
	ok = 1;
      }
    }
    else {
      char *id = strtok( buffer, " 	" );
      char *pheno = strtok( NULL,  " 	" );
      char *endptr;
      q->name[q->N] = (char*)strdup(id);
      q->observed[q->N] = strtod(pheno,&endptr);
      if (endptr == pheno)
	q->observed[q->N] = NaN;
      ok = 1;
    }

    if ( ok ) {
      if ( verbose >=2 ) 
	Rprintf("individual %s %.5f\n", q->name[q->N], q->observed[q->N] );
      while( (str1 = strtok( NULL, "	 " ) ) && (str2 = strtok( NULL, " 	" ) ) ) {
	if ( m >= q->M ) {
	  Rprintf( "ERROR: too many markers on line %d\n", q->N );
	  error("fatal HAPPY error");
	}
	if ( ! legal_string( str1, a->af[m].allele_name, a->af[m].alleles, &q->genos[q->N].chrom1[m] ) ) {
	  int k;
	  Rprintf( "ERROR: subject %s unknown allele1 %s for marker %d %s - legal values are ", q->name[q->N], str1, m, a->af[m].marker_name );
	  for(k=0;k<a->af[m].alleles;k++)
	    Rprintf( " %s", a->af[m].allele_name[k] );
	  Rprintf( "\n");
	  if ( ! legal_string( missingCode, a->af[m].allele_name, a->af[m].alleles, &q->genos[q->N].chrom1[m] ) ) {
	    Rprintf( "ERROR: subject %s unknown allele1 %s for marker %d %s - legal values are",  q->name[q->N], missingCode, m, a->af[m].marker_name );
	    error("fatal HAPPY error");
	  }
	}
	if ( strcmp( str1, a->af[m].allele_name[q->genos[q->N].chrom1[m]] ) ) {
	  Rprintf( "ERROR subject %s decoding allele %s %s\n",  q->name[q->N], str1, a->af[m].allele_name[q->genos[q->N].chrom1[m]] );
	}

	if ( ! legal_string( str2, a->af[m].allele_name, a->af[m].alleles, &q->genos[q->N].chrom2[m] ) ) {
	  int k;
	  Rprintf( "ERROR: subject %s unknown allele2 %s for marker %d %s - legal values are",  q->name[q->N], str2, m,  a->af[m].marker_name );
	  for(k=0;k<a->af[m].alleles;k++)
	    Rprintf( " %s", a->af[m].allele_name[k] );
	  Rprintf( "\n");
	  if ( ! legal_string( missingCode, a->af[m].allele_name, a->af[m].alleles, &q->genos[q->N].chrom2[m] ) ) {
	    Rprintf( "ERROR:subject %s  unknown allele2 %s for marker %d %s - legal values are",  q->name[q->N], missingCode, m, a->af[m].marker_name );
	    error("fatal HAPPY error");
	  }
	}
	if ( strcmp( str2, a->af[m].allele_name[q->genos[q->N].chrom2[m]] )) {
	  Rprintf( "ERROR subject %s decoding allele %s %s\n",  q->name[q->N], str2, a->af[m].allele_name[q->genos[q->N].chrom2[m]] );
	}
	/*Rprintf( "%d: %s %s\n", m, str1, str2 ); */
	a->af[m].allele_freq[q->genos[q->N].chrom1[m]]++;
	a->af[m].allele_freq[q->genos[q->N].chrom2[m]]++;
	m++;
      }
    }
    if ( m < q->M ) {
      Rprintf( "ERROR subject %s does not have enough alleles %d (%d) on line %d\n",  q->name[q->N], m, q->M, q->N );
      error("fatal HAPPY error");
    }
    q->N++;
  }

  if ( verbose>=2 ) {
    for(m=0;m<q->M;m++) {
      int al;
      ALLELE_FREQ *af = &a->af[m];
      Rprintf( "marker %s %.3f\n", af->marker_name, af->position );
      for(al=0;al<af->alleles;al++) 
	Rprintf( "%10s %5.0f\n", af->allele_name[al], af->allele_freq[al] );
    }
  }

  for(m=0;m<q->M;m++) {
    ALLELE_FREQ *af = &a->af[m];
    int al;
    for(al=0;al<af->alleles;al++) 
      af->allele_freq[al] /= (2*q->N);
  }

  for(m=0;m<q->M;m++) {
    ALLELE_FREQ *af = &a->af[m];
    int al, s;
    for(s=0;s<q->S;s++) {
      double p = 0.0;
      for(al=0;al<af->alleles;al++) 
	p += af->pr_AtoS[al][s]*af->allele_freq[al];
      if ( p > 1.0e-7 ) 
	af->entropy -= p*log(p);
    }
    /*    printf( "marker %d %s entropy %e\n", m+1, af->marker_name, af->entropy ); */
  }

  q->fit = (QTL_FIT*)calloc(q->M,sizeof(QTL_FIT));
  for(m=0;m<q->M;m++) 
    allocate_qtl_fit( &q->fit[m], q->N, q->S );
  

  Rprintf( "Number of individuals: %-5d\n", q->N );
  Rprintf( "Number of markers:     %-5d\n", q->M );
  Rprintf( "Number of strains:     %-5d\n", q->S );

  Rprintf( "Use Parents:           %s\n", q->use_parents ? "yes" : "no" );

#ifdef _USE_HASH_

  Rprintf( "status %d\n", use_parents || ped_format );
  if ( use_parents || ped_format ) {
    int i, nparents=0;
    PARENT_KEY e;
    PARENT_KEY *f;
    int both = 0;

    q->mother = (int*)calloc(q->N, sizeof(int));
    q->father = (int*)calloc(q->N, sizeof(int));
    pcount = (int*)calloc(q->N, sizeof(int));
    
    if ( hcreate(q->N) == 0 ) 
      error("Could not create hash table");

    for(i=0;i<q->N;i++) {
      e.key = q->name[i];
      e.data = (void*)(long)i;
      hsearch( e, ENTER );
    }

    for(i=0;i<q->N;i++) {
      e.key = mother[i];
      f = hsearch( e, FIND );
      if ( f != NULL )
	q->mother[i] = (int)((long)(f->data));
      else
	q->mother[i] = -1;

      e.key = father[i];
      f = hsearch( e, FIND );
      if ( f != NULL )
	q->father[i] = (int)((long)(f->data));
      else 
	q->father[i] = -1;
      if ( q->mother[i] > -1 && q->father[i] > -1 ) {
	both++;
	pcount[q->mother[i]]++;
	pcount[q->father[i]]++;
      }
    }

    for( i=0;i<q->N;i++)
      nparents += pcount[i]>0

    hdestroy();
    free(mother);
    free(father);
    free(pcount);

    Rprintf( "Number of subjects with two parents: %-5d\n", both );
    Rprintf( "Number of parents in nuclear families: %-5d\n", nparents );

  }
#else

  if ( use_parents || ped_format ) {
    int i;
    int both = 0;
    PARENT_KEY *sorted = (PARENT_KEY*)calloc(q->N, sizeof(PARENT_KEY));
    PARENT_KEY *m, *f;

    q->mother = (int*)calloc(q->N, sizeof(int));
    q->father = (int*)calloc(q->N, sizeof(int));
    pcount = (int*)calloc(q->N, sizeof(int));

    for(i=0;i<q->N;i++) {
      sorted[i].key = q->name[i];
      sorted[i].id = i;
    }

    qsort( sorted, q->N, sizeof(PARENT_KEY), entrycmp );


    for(i=0;i<q->N;i++) {
      PARENT_KEY mum, dad;
      mum.key = mother[i];
      dad.key = father[i];
      m = (PARENT_KEY*)bsearch( &mum, sorted, q->N, sizeof(PARENT_KEY), entrycmp );
      f = (PARENT_KEY*)bsearch( &dad, sorted, q->N, sizeof(PARENT_KEY), entrycmp );
      if ( m != NULL )
	q->mother[i] = m->id;
      else
	q->mother[i] = -1;

      if ( f != NULL )
	q->father[i] = f->id;
      else 
	q->father[i] = -1;
      if ( q->mother[i] > -1 && q->father[i] > -1 ) {
	both++;
	pcount[q->mother[i]]++;
	pcount[q->father[i]]++;
      }

    }

    for( i=0;i<q->N;i++)
      nparents += pcount[i]>0;

    free(mother);
    free(father);
    free(sorted);
    free(pcount);

    Rprintf( "Number of subjects with two parents: %-5d\n", both );
    Rprintf( "Number of parents in nuclear families: %-5d\n", nparents );

  }
#endif


  fit_null_qtl_model( q );

  free(buffer);
  return q;

}

double fit_null_qtl_model( QTL_DATA *q ) {

  int i, k;
  QTL_FIT *fit = q->null_model = (QTL_FIT*)calloc(1,sizeof(QTL_FIT));
  allocate_qtl_fit( fit, q->N, q->S );

  fit->mean = 0.0;
  for(i=0;i<q->N;i++) {
    fit->mean += q->observed[i];
    fit->sigma += q->observed[i]*q->observed[i];
  }
  fit->mean /= q->N;
  fit->sigma = ( fit->sigma - q->N*fit->mean*fit->mean )/(q->N-1);

  fit->rss = 0.0;
  for(i=0;i<q->N;i++) {
    double residual = q->observed[i] - fit->mean;
    fit->rss += residual*residual;
  }

  for(k=0;k<q->alleles->strains;k++) 
    fit->trait[k] = fit->trait_error[k] = 0.0;

  printf("null model mean %e var %e\n", fit->mean, fit->sigma );
  return fit->sigma;
}


ANCESTRY *read_subject_ancestries( FILE *fp, char *filename, int verbose ) {

  if ( fp ) {
    int subjects = 0;
    int strains = 0;
    char line[256];
    int line_no = 0;
    int i;
    int imax = 10000;

    Rprintf( "Reading subject ancestries from %s\n", filename );
    skip_comments(fp, line);
    line_no++;
    if ( sscanf( line, "subjects %d strains %d", &subjects, &strains ) == 2 ) {
      Rprintf( "subjects %d strains %d", subjects, strains );
      ANCESTRY *an = (ANCESTRY*)calloc(1, sizeof(ANCESTRY));
      an->N = subjects;
      an->S = strains;
      skip_comments( fp, line ); line_no++;
      if ( 0 == strncmp( line, "strain_names", strlen("strain_names") ) ) {
	int k;
	char *str=strtok(line,"	 ");
	an->strain_name = (char**)calloc(strains,sizeof(char*));
	for(k=0;k<strains;k++) {
	  if ( str=strtok(NULL," 	") ) {
	    an->strain_name[k] = (char*)strdup(str);
	  }
	  else {
	    Rprintf( "ERROR not enough strain names %d/%d\n", k, strains );
	    error("fatal HAPPY error");
	  }
	}
      }

      an->subject_name = (char**)calloc(subjects, sizeof(char*));
      an->prob = (double**)calloc(subjects,sizeof(double*));
      for(i=0;i<subjects;i++) {
	char *str;
	int k;
	double total = 1.0e-10;;
	line[0] = 0;
	skip_comments(fp, line);

	line_no++;
	str=strtok(line,"	 ");
	an->subject_name[i] = (char*)strdup(str);
	an->prob[i] = (double*)calloc(strains, sizeof(double));
	for(k=0;k<strains;k++) {
	  double x = 0;
	  if ( (str = (char*)strtok( NULL, "	 " )) && sscanf( str, "%lf", &x ) == 1 ) {
	    if ( x < 0.0 ) {
	      Rprintf( "ERROR negative ancestry probability %lf line %d, set to 0\n", x, line_no );
	      x = 0;
	    }
	    an->prob[i][k] = x;
	    total += x;
	  }
	  else {
	    Rprintf( "ERROR not a probability \"%s\" (token %d) in ancestry file line %d\n", str, k, line_no );
	    error( "fatal HAPPY error");
	  }
	}
	for(k=0;k<strains;k++)
	  an->prob[i][k] /= total;
      }
      return(an);
    }
  }
  return(NULL);
}

void heterozygosity( QTL_DATA *q ) {
  int id, m;
  ALLELES *A = q->alleles;

  for( id=0;id<q->N;id++) {
    double het = subject_heterozygosity(q, id);
    if ( het > 0.0 ) {
      Rprintf( "subject %20.20s heterozygosity %.4f\n", q->name[id], het );
    }
  }
 
  for( m=0;m<q->M;m++) {
    double het = marker_heterozygosity(q, m);
    if ( het > 0.0 ) {
      Rprintf( "marker %20.20s heterozygosity %.4f\n",  A->af[m].marker_name, het );
    }
  }
}

double subject_heterozygosity( QTL_DATA *q, int individual ) {

  double het = 0.0;
  int m;
  for(m=0;m<q->M;m++) {
    het += q->genos[individual].chrom1[m]!=q->genos[individual].chrom2[m];
  }
  het /= q->M;
  return(het);
}

double marker_heterozygosity( QTL_DATA *q, int marker ) {

  double het = 0.0;
  int i;
  for(i=0;i<q->N;i++) {
    het += q->genos[i].chrom1[marker]!=q->genos[i].chrom2[marker];
  }
  het /= q->N;
  return(het);
}

int check_and_apply_ancestry(QTL_DATA *q ) {
  ANCESTRY *an = q->an;

  if ( an ) {
    ALLELES *A = q->alleles;
    int i;
    if ( an->S != q->S ) {
      Rprintf( "ERROR number of strains in ancestry file %d unequal to number of strains in alleles file %d\n", an->S, q->S );
      error( "fatal HAPPY error");
    }
    else {
      int i;
      int bad = 0;
      for ( i=0; i<q->S; i++) {
	if ( strcmp( an->strain_name[i], A->strain_name[i] ) != 0) {
	  bad ++;
	  Rprintf( "ERROR strain at position %d name %s in ancestry differs from %s in alleles\n", i+1, an->strain_name[i], A->strain_name[i] );
	}
	if ( bad > 0 ) 
	  error( "fatal HAPPY error");
      }
      Rprintf( "Checked consistency of strain names between ancestry and alleles: OK\n" );
    }

    if ( an->N != q->N ) {
      Rprintf( "ERROR number of subjects in ancestry file %d unequal to number of subjects in alleles file %d\n", an->N, q->N );
      error( "fatal HAPPY error");
    }
    else {
      int i;
      int bad = 0;
      for ( i=0; i<q->N; i++) {
	if ( strcmp( an->subject_name[i], q->name[i] ) != 0) {
	  bad ++;
	  Rprintf( "ERROR subject at position %d name %s in ancestry differs from %s in data\n", i+1, an->subject_name[i], q->name[i] );
	}
	if ( bad > 0 ) 
	  error( "fatal HAPPY error");
      }
      Rprintf( "Checked consistency of subject names between ancestry and data: OK\n" );
    }

    an->pr_AtoS = (double****)calloc(an->N, sizeof(double***));
    for ( i=0; i<q->N; i++) {
      int m;
      an->pr_AtoS[i] = (double***)calloc(q->M, sizeof(double**));
      for(m=0;m<q->M;m++) {
	ALLELE_FREQ *af = &(A->af[m]);
	double **pr_AtoS = af->pr_AtoS;
	int s;
	int a;
	an->pr_AtoS[i][m] = (double**)calloc(af->alleles,sizeof(double*));
	for(a=0;a<af->alleles;a++) {
	  double total = 1.0e-10;
	  an->pr_AtoS[i][m][a] = (double*)calloc(q->S,sizeof(double));
	  for(s=0;s<q->S;s++) {
	    total += an->prob[i][s]*pr_AtoS[a][s];
	  }
	  for(s=0;s<q->S;s++) {
	    an->pr_AtoS[i][m][a][s] = an->prob[i][s]*pr_AtoS[a][s] / total;
	  }
	}
      }
    }

    return(1);
  }
  return(0);
}

ALLELES *input_allele_frequencies( FILE *fp, int generations, char *missingCode, double MinDist, int verbose ) {

  /* Example format:

  [ all probabilities are expressed as Pr(strain|allele) ]

  [rows are alleles, columns are strains]

  markers 5 strains 8 strain_name-1 strain-name2 ... strain-name8
  marker m1 4 81.3
  allele 1   0.5  0.0  0.5  0.0  0.0  0.0  0.0  0.0
  allele 2   0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
  allele 3   0.0  0.0  0.0  0.33 0.33 0.0  0.0  0.33 
  allele 4   0.0  0.0  0.0  0.0  0.0  0.5  0.5  0.0
  marker m2 4 85.5
  allele 1   0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
  allele 2   0.25 0.25 0.0  0.0  0.25 0.25 0.0  0.0
  allele 3   0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
  allele 4   0.0  0.0  0.0  0.0  0.0  0.0  0.5  0.5
  marker m3 5 86.0
  allele 1   0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0
  allele 2   0.0  0.0  0.5  0.0  0.5  0.0  0.0  0.0
  allele 3   1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
  allele 4   0.0  0.0  0.0  0.0  0.0  0.0  0.5  0.5
  allele 5   0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
  marker m4 2 87.5
  allele 1   0.2  0.2  0.2  0.0  0.0  0.0  0.2  0.2
  allele 2   0.0  0.0  0.0  0.33 0.33 0.33 0.0  0.0
  marker m5  4 87.8
  allele 1   0.0  0.0  0.0  0.5  0.0  0.0  0.0  0.5
  allele 2   0.0  0.0  0.0  0.0  0.5  0.5  0.0  0.0
  allele 3   0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
  allele 4   0.5  0.0  0.5  0.0  0.0  0.0  0.5  0.0

  Note the marker location in centimorgans can be omitted, in which case it defaults to 1.0 cM
  */


  /* long temporay strings to avoid memory errors by scanf */
  char line[10000];
  char marker_name[10000];
  char allele_name[10000];
  char chromosome_name[10000];
  int markers=-1;
  int strains=-1;
  int line_no = 0;
  int m, s, a;
  ALLELES *A=NULL;

  skip_comments( fp, line ); line_no++;
  if (strlen(line)>1000) {
   Rprintf( "Found length of line number %d to be too long (> 1000). Check line separator.\n", line_no);
   error( "Found length of line to be too long (> 1000). Check line separator.\n");
  }

  if ( sscanf( line, "markers %d strains %d", &markers, &strains ) == 2 ) {
    A = (ALLELES*)calloc(1,sizeof(ALLELES));
    A->markers = markers;
    A->strains = strains;
    A->generations = generations;
    A->af = (ALLELE_FREQ*)calloc(markers,sizeof(ALLELE_FREQ));
    skip_comments( fp, line ); line_no++;

    if ( 0 == strncmp( line, "strain_names", strlen("strain_names") ) ) {
      int k;
      char *str=strtok(line,"	 ");
      A->strain_name = (char**)calloc(strains,sizeof(char*));
      for(k=0;k<strains;k++) {
	if ( str=strtok(NULL," 	") ) {
	  A->strain_name[k] = (char*)strdup(str);
	}
	else {
	  Rprintf( "ERROR not enough strain names %d/%d\n", k, strains );
	  error("fatal HAPPY error");
	}
      }
    }

    for(m=0;m<markers;m++) {
      ALLELE_FREQ *af = &(A->af[m]);
      double total = 0.0;
      skip_comments( fp, line ); line_no++;

      if ( sscanf( line, "marker %s %d %s %lf", marker_name, &(af->alleles), chromosome_name, &(af->position) ) == 4 ) {
	strncpy(af->chromosome,chromosome_name,MAX_LENGTH_CHROMOSOME);
      }
      else if ( sscanf( line, "marker %s %d %lf", marker_name, &(af->alleles), &(af->position) ) == 3 ) {
	strncpy(af->chromosome,"unknown",MAX_LENGTH_CHROMOSOME);
      }
      else {
	Rprintf( "marker Parse ERROR, line %d %s\n", line_no, line );
	error("fatal HAPPY error");
      }

      strncpy(af->marker_name,marker_name,MAX_LENGTH_MARKER_NAME);

	if ( verbose>=2 ) 
	  Rprintf("marker %d %s %.3f %d\n", m, af->marker_name, af->position, af->alleles);
	af->pr_AtoS = (double**)calloc(af->alleles,sizeof(double*)); 
	for(s=0;s<af->alleles;s++)
	  af->pr_AtoS[s] = (double*)calloc(strains,sizeof(double)); 
	af->allele_name = (char**)calloc(af->alleles,sizeof(char*));
	af->which_allele = (int*)calloc(strains,sizeof(int));
	af->allele_freq = (double*)calloc(af->alleles,sizeof(double));
	/*	printf( "marker %d %s\n", m+1, marker_name ); */
	af->ND = -1; /* index of code for missing value */
	for(a=0;a<af->alleles;a++) {
	  skip_comments( fp, line ); line_no++;
	  if ( sscanf( line, "allele %s", allele_name ) != 1 ) {
	    Rprintf( "allele Parse ERROR, line %d %s\n", line_no, line );
	    error("fatal HAPPY error");
          }
	    char *str = (char*)strtok( line, "	 " );
	    af->allele_name[a] = (char*)strdup(allele_name);
	    str = (char*)strtok( NULL, " 	" );
	    if ( ! strncmp(allele_name, missingCode, MAX_LENGTH_MISSINGCODE) )
	      af->ND = a;
	    /*	    Rprintf( "missing code for %d %d\n", m, af->ND);*/
	    total = 0.0;
	    for(s=0;s<strains;s++) {
	      double x;
	      if ( (str = (char*)strtok( NULL, "	 " )) && sscanf( str, "%lf", &x ) == 1 ) {
		if ( x <= 0.0 ) 
		  x = 1.0e-7;
		else 
		  af->which_allele[s] = a;
		total += af->pr_AtoS[a][s] = x; /* prob allele a induces strain s  ie Pr(s|a) */
	      }
	      else {
		Rprintf( "probability Parse ERROR, line %d token :%s:\n", line_no, str );
		error("fatal HAPPY error");
	      }
	    }
	    for(s=0;s<strains;s++) {
	      af->pr_AtoS[a][s] /= total;
	    }
	} // EoFor alleles
    } // EoFor markers

    A->Pr_ss = (double*)calloc(markers,sizeof(double));
    A->Pr_st = (double*)calloc(markers,sizeof(double));

    for(m=1;m<markers;m++) {
      if ( !strcmp(A->af[m].chromosome, A->af[m-1].chromosome)) {
	double d = (A->af[m].position - A->af[m-1].position)/100.0;
	double lambda = generations*d; /* Poisson parameter */
	double q = exp(-lambda);
	double p = A->af[m].ProbSame = q + (1.0-q)/strains;
	A->Pr_ss[m-1] = p;
	A->Pr_st[m-1] = (1.0-p)/(strains-1.0);
	if ( verbose>=2 ) 
	  Rprintf("marker %d %s %.3f %.4f %.4f %.4f \n", m, A->af[m].marker_name, A->af[m].position, p, A->Pr_ss[m-1], A->Pr_st[m-1] );
      }
      else {
	A->af[m].ProbSame = 1.0/strains;
	A->Pr_ss[m-1] = 1.0/strains;
	A->Pr_st[m-1] = 1.0/strains;
      }
    }

    /* compute the prior probabilities that, if a QTL lies between
       marker m-1 and m, the linkage to the neighbouring markers for
       the two chromosomes are X,Y where X and Y must take exactly one
       of the values BOTH, LEFT, UNLINKED, RIGHT */
    A->MinDist = MinDist;
    for(m=0;m<markers-1;m++) {
      ALLELE_FREQ *af = &(A->af[m]);
      int s, t;
      double total;
      double P[LINKAGE_STATES];
      double d = (A->af[m+1].position - A->af[m].position)/100.0;
      double lambda, elambda, elambda2;


      if ( d < MinDist ) d = MinDist;
      lambda = generations*d;
      elambda = exp(-lambda);
      elambda2 = elambda*elambda;
      P[BOTH] = elambda;
      P[LEFT] = P[RIGHT] = (1-elambda)/lambda-elambda;
      P[UNLINKED] = 1.0-P[BOTH]-P[LEFT]-P[RIGHT];
      af->prior = (double**)calloc(LINKAGE_STATES,sizeof(double*));
      for(s=0;s<LINKAGE_STATES;s++)
	af->prior[s] = (double*)calloc(LINKAGE_STATES,sizeof(double));
      af->prior[BOTH][BOTH] = P[BOTH]*P[BOTH];
      af->prior[BOTH][RIGHT] = af->prior[RIGHT][BOTH] = af->prior[BOTH][LEFT] = af->prior[LEFT][BOTH] = P[BOTH]*P[LEFT];
      af->prior[BOTH][UNLINKED] = af->prior[UNLINKED][BOTH] = P[BOTH] - af->prior[BOTH][BOTH] - af->prior[BOTH][LEFT] - af->prior[BOTH][RIGHT];
      af->prior[LEFT][LEFT] = af->prior[RIGHT][RIGHT] = 0.5*(1-elambda2)/lambda - 2*elambda*(1-elambda)/lambda + elambda2;
      af->prior[LEFT][RIGHT] = af->prior[RIGHT][LEFT] = P[BOTH]*(1-2*P[LEFT]-P[BOTH]);
      af->prior[LEFT][UNLINKED] = af->prior[RIGHT][UNLINKED] = af->prior[UNLINKED][LEFT] = af->prior[UNLINKED][RIGHT] = 
	P[LEFT] - af->prior[BOTH][LEFT] - af->prior[LEFT][LEFT] - af->prior[RIGHT][LEFT];
      af->prior[UNLINKED][UNLINKED] = P[UNLINKED] - af->prior[UNLINKED][BOTH] - af->prior[UNLINKED][LEFT] - af->prior[UNLINKED][RIGHT];
      if ( verbose>=2 ) {
	Rprintf("priors for %d lambda %.5f d %.5f\n", m, lambda, d);
	for(s=BOTH;s<=RIGHT;s++) {
	  total = 0.0;
	  Rprintf("%2d  ", s);
	  for(t=BOTH;t<=RIGHT;t++) {
	    Rprintf(" %.5f", af->prior[s][t] );
	    total += af->prior[s][t];
	  }
	  Rprintf("   %.5f %.5f\n", total, P[s]);
	}
      }
    }
    
  }
  return A;
}



int KVcmp( const void *a, const void *b ) {

  const KV *A = (const KV*)a;
  const KV *B = (const KV*)b;
  double x = A->key-B->key;
  if ( x > 0.0 )
    return 1;
  else if ( x < 0.0 ) 
    return -1;
  else
    return 0.0;
}
 


/* compute the matrix of probabilities that, if the qtl is a distance
   c from the left end of the curent interval, the pair of linkage
   states for the two chromosomes takes the value s,t */

void pointwise_interval_mapping_probabilities( QTL_DATA *q, int locus, double c, double **prior ) {

  ALLELES *A = q->alleles;
  double d = (A->af[locus+1].position - A->af[locus].position)/100.0;
  double lambda = A->generations*d;
  double elambda = exp(-lambda);
  double elambdac = exp(-lambda*c);
  double elambdac1 = exp(-lambda*(1-c));
  int s, t;
  double P[LINKAGE_STATES];
  P[BOTH] = elambda;
  P[LEFT] = elambdac-elambda;
  P[RIGHT] = elambdac1 -elambda;
  P[UNLINKED] = 1.0-P[BOTH]-P[LEFT]-P[RIGHT];

  
  for(s=0;s<LINKAGE_STATES;s++)
    if ( P[s] < 0.0 ) P[s] = 0.0;

  for(s=0;s<LINKAGE_STATES;s++)
    for(t=0;t<LINKAGE_STATES;t++)
      prior[s][t] = P[s]*P[t];

}

QTL_PRIOR ***compute_qtl_priors(  QTL_DATA *qtl, QTL_PRIOR ***qp, int locus, double **P ) {

  /* compute, for each individual, the prior probabilities of the
     different strains at the flanking markers */

  int i;
  int m = locus;
  double total;
  double *LeftMargin = (double*)calloc(qtl->S,sizeof(double));
  double *RightMargin = (double*)calloc(qtl->S,sizeof(double));
  double s1 = 1.0/qtl->S;

  for(i=0;i<qtl->N;i++) {
    int s, t;
    DP_MATRICES *dp = &(qtl->dp_matrices[i]);
    double **Left = dp->Left[m];
    double **Right = dp->Right[m+1];
    dp->NonRecomb[m]=0.0; 
    for(s=0;s<qtl->S;s++) {
      LeftMargin[s] = RightMargin[s] = 0.0;
      for(t=0;t<qtl->S;t++) {
	LeftMargin[s] += Left[s][t];
	RightMargin[s] += Right[s][t];
      }
    }
    total = 0.0;
    for(s=0;s<qtl->S;s++) {
      for(t=0;t<qtl->S;t++) {
	total += qp[i][s][t].prior = 
	  Left[s][t]*Right[s][t]*P[BOTH][BOTH] + Left[s][t]*RightMargin[t]*P[LEFT][BOTH] + LeftMargin[t]*RightMargin[t]*P[UNLINKED][BOTH]*s1 + LeftMargin[t]*Right[s][t]*P[RIGHT][BOTH] +
	  Left[s][t]*RightMargin[s]*P[BOTH][LEFT] + Left[s][t]*P[LEFT][LEFT] + LeftMargin[t]*P[UNLINKED][LEFT]*s1 + LeftMargin[t]*RightMargin[s]*P[RIGHT][LEFT] +
	  LeftMargin[s]*RightMargin[s]*P[BOTH][UNLINKED]*s1 + LeftMargin[s]*P[LEFT][UNLINKED]*s1 + P[UNLINKED][UNLINKED]*s1*s1 + RightMargin[s]*P[RIGHT][UNLINKED]*s1 +
	  LeftMargin[s]*Right[s][t]*P[BOTH][RIGHT] + LeftMargin[s]*RightMargin[t]*P[LEFT][RIGHT] + RightMargin[t]*P[UNLINKED][RIGHT]*s1 + Right[s][t]*P[RIGHT][RIGHT]; 
	dp->NonRecomb[m] += 2*Left[s][t]*Right[s][t]*P[BOTH][BOTH] + Left[s][t]*RightMargin[t]*P[LEFT][BOTH] + LeftMargin[t]*RightMargin[t]*P[UNLINKED][BOTH]*s1 + LeftMargin[t]*Right[s][t]*P[RIGHT][BOTH] + Left[s][t]*RightMargin[s]*P[BOTH][LEFT] +  LeftMargin[s]*RightMargin[s]*P[BOTH][UNLINKED]*s1 + LeftMargin[s]*Right[s][t]*P[BOTH][RIGHT]; 
	/*	Rprintf("%d %d %lf\n", i, m, dp->NonRecomb[m]);*/
	
      }
    }
    /*    printf( "new total %d %.8f\n", i, total ); */
    for(s=0;s<qtl->S;s++) 
      for(t=0;t<qtl->S;t++) {
	qp[i][s][t].prior /= total;
	/*	if ( !(qp[i][s][t].prior < 1.0) ) 	qp[i][s][t].prior = s12; */
	/*	Rprintf( "%d %d %d %lf %d %d\n", i, s, t, qp[i][s][t].prior, qp[i][s][t].prior==NAN,  */
      }
    dp->NonRecomb[m] /= total;
  }

  free(LeftMargin);
  free(RightMargin);
  return qp;
}



QTL_PRIOR ***allocate_qtl_priors( QTL_DATA *qtl ) {
  
  QTL_PRIOR ***qp = (QTL_PRIOR***)calloc(qtl->N,sizeof(QTL_PRIOR**));
  int i, k;
  
  for(i=0;i<qtl->N;i++) {
    qp[i] = (QTL_PRIOR**)calloc(qtl->S,sizeof(QTL_PRIOR*));
    for(k=0;k<qtl->S;k++)
      qp[i][k] = (QTL_PRIOR*)calloc(qtl->S,sizeof(QTL_PRIOR));
  }
  
  return qp;
}



/* create the forward and backward Dynamic Programming matrices for the 
   sum of all possible haplotype reconstructions */
  
void create_summed_dp_matrices( QTL_DATA *q ) {

  int i;
  double *Pr_ss = q->alleles->Pr_ss;
  double *Pr_st = q->alleles->Pr_st;


  q->dp_matrices = (DP_MATRICES*)calloc(q->N,sizeof(DP_MATRICES));

  for(i=0;i<q->N;i++) {

    if ( i == 0 || genotype_difference( q, i, i-1 ) ) {
      q->dp_matrices[i].Left = summed_dp_matrix( q, i, Pr_ss, Pr_st, +1 ); 
      q->dp_matrices[i].Right = summed_dp_matrix( q, i, Pr_ss, Pr_st, -1 ); 
      q->dp_matrices[i].NonRecomb = (double*)calloc( q->M, sizeof(double)); 
    }
    else {
      /*      printf( "same genotype %d %d\n"); */
      q->dp_matrices[i].Left = q->dp_matrices[i-1].Left;
      q->dp_matrices[i].Right = q->dp_matrices[i-1].Right;
      q->dp_matrices[i].NonRecomb = q->dp_matrices[i-1].NonRecomb; 
    }
  }
}


double ***summed_dp_matrix( QTL_DATA *qtl, int individual, double *Pr_ss, double *Pr_st, int direction ) {

  double ***X;
  int start, stop, incr;
  CHROM_PAIR *genotypes = &(qtl->genos[individual]);
  ALLELES *A = qtl->alleles;
  int m, s, t;
  int offset;
  double total;
  int markers = genotypes->markers;
  int strains = A->strains;
  double **Pr_tr1;
  double **Pr_tr2;
  double root2 = sqrt(2.0);
  int mum = -1;
  int dad = -1;
  CHROM_PAIR *mgenotypes=NULL;
  CHROM_PAIR *fgenotypes=NULL;
  int use_parents = 0;


  /*  double dp=0.0;
  double ndp=1.0e-10;
  */
  if ( qtl->use_parents && (mum=qtl->mother[individual]) > -1 && (dad=qtl->father[individual]) > -1 ) {
    use_parents = 1;
    mgenotypes = &(qtl->genos[mum]);
    fgenotypes = &(qtl->genos[dad]);
  }


  /*  printf( "dp_matrix %d strains %d markers %d\n", direction, strains, markers ); */
  /*  Rprintf( "use_parents %d %d individual %d mum %d dad %d\n", use_parents,qtl->use_parents, individual, mum, dad ); */

  Pr_tr1 = (double**)calloc(strains,sizeof(double*));
  for(s=0;s<strains;s++) 
    Pr_tr1[s] = (double*)calloc(strains,sizeof(double));

  Pr_tr2 = (double**)calloc(strains,sizeof(double*));
  for(s=0;s<strains;s++) 
    Pr_tr2[s] = (double*)calloc(strains,sizeof(double));

  X = (double***)calloc(markers,sizeof(double**));
  for(m=0;m<markers;m++) {
    X[m] = (double**)calloc(strains,sizeof(double*));
    for(s=0;s<strains;s++) 
      X[m][s] = (double*)calloc(strains,sizeof(double));
  }

  if ( direction > 0 ) {
    start = 0;
    stop = markers-1;
    incr = +1;
    offset = 0;
  }
  else {
    start = markers-1;
    stop = 0;
    incr = -1;
    offset = -1;
  }

  if ( qtl->phase_known ) {
    double **a;
    if ( qtl->an ) 
      a = qtl->an->pr_AtoS[individual][start];
    else
      a = A->af[start].pr_AtoS;

    for(s=0;s<strains;s++)
      for(t=0;t<strains;t++) {
	int g1 = genotypes->chrom1[start];
	int g2 = genotypes->chrom2[start];
	X[start][s][t] = a[g1][s] * a[g2][t];
      }
  }
  else {
    double phaseP = 1.0;
    int g1 = genotypes->chrom1[start];
    int g2 = genotypes->chrom2[start];
    double **a;
    if ( use_parents ) {
      int m1 = mgenotypes->chrom1[start];
      int m2 = mgenotypes->chrom2[start];
      int p1 = fgenotypes->chrom1[start];
      int p2 = fgenotypes->chrom2[start];
      phaseP = phaseProb( g1, g2, m1, m2, p1, p2, A->af[m].ND );
    }
    if ( qtl->an ) 
      a = qtl->an->pr_AtoS[individual][start];
    else {
	ALLELE_FREQ *af = &(A->af[start]);
	a = af->pr_AtoS;
    }

    for(s=0;s<strains;s++)
      for(t=0;t<strains;t++) {

	X[start][s][t] = phaseP*a[g1][s] * a[g2][t] + (2.0-phaseP)*a[g2][s] * a[g1][t];
      }
  }

  for(m=start+incr;m!=stop;m+=incr) { 
    int g1 = genotypes->chrom1[m];
    int g2 = genotypes->chrom2[m];
    double pr_ss = Pr_ss[m+offset];
    double pr_st = Pr_st[m+offset];
    double **a;
    if ( qtl->an ) 
      a = qtl->an->pr_AtoS[individual][m];
    else
      a = A->af[m].pr_AtoS;

    int m1=0, m2=0, p1=0, p2=0;
    if ( use_parents ) {
      m1 = mgenotypes->chrom1[m];
      m2 = mgenotypes->chrom2[m];
      p1 = fgenotypes->chrom1[m];
      p2 = fgenotypes->chrom2[m];
    }
    for(s=0;s<strains;s++) {
      double norm1 = 1.0e-10;
      double norm2 = 1.0e-10;
      for(t=0;t<strains;t++) {
	if ( s == t ) {
	  norm1 += Pr_tr1[s][t] = pr_ss*a[g1][t];
	  norm2 += Pr_tr2[s][t] = pr_ss*a[g2][t];
	}
	else {
	  norm1 += Pr_tr1[s][t] = pr_st*a[g1][t];
	  norm2 += Pr_tr2[s][t] = pr_st*a[g2][t];
	}
      }
      if ( ! qtl->phase_known ) {
	norm1 *= root2;
	norm2 *= root2;
      }

      for(t=0;t<strains;t++) {
	Pr_tr1[s][t] /= norm1;
	Pr_tr2[s][t] /= norm2;
      }
    }
    
    total = 0.0;
    if ( qtl->phase_known ) {
      for(s=0;s<strains;s++)
	for(t=0;t<strains;t++) {
	  int S, T;
	  for(S=0;S<strains;S++) {
	    for(T=0;T<strains;T++) 
	      X[m][s][t] += X[m-incr][S][T] *  Pr_tr1[S][s] * Pr_tr2[T][t];
	  }
	}
    }
    else {
      
      double phaseP = 1.0;
      if ( use_parents ) {
	phaseP = 2*phaseProb( g1, g2, m1, m2, p1, p2, A->af[m].ND );
      }
      for(s=0;s<strains;s++)
	for(t=0;t<strains;t++) {
	  int S, T;
	  for(S=0;S<strains;S++) {
	    for(T=0;T<strains;T++) 
	      X[m][s][t] += X[m-incr][S][T] * ( Pr_tr1[S][s] * Pr_tr2[T][t] * phaseP + Pr_tr2[S][s] * Pr_tr1[T][t] * (2.0-phaseP) );
	  }
	}
    }
  }

  for(s=0;s<strains;s++) {
    free(Pr_tr1[s]);
    free(Pr_tr2[s]);
  }
  free(Pr_tr1);
  free(Pr_tr2);

  return X;
}	      





int genotype_difference( QTL_DATA *q, int i, int j ) {

  int d=0;
  if ( i >=0 && i < q->N && j >= 0 && j < q-> N ) {
    int m;
    for(m=0;m<q->M;m++) 
      d += (q->genos[i].chrom1[m] != q->genos[j].chrom1[m]) + 
	(q->genos[i].chrom2[m] != q->genos[j].chrom2[m]);
  }  
  else {
    d = -1;
  }
  /*  printf( "comparing %d %d diff %d\n", i, j, d); */
  return d;
}






int marker_index( const char *name, QTL_DATA *q, const int isIntervalModel )
{
    ALLELE_FREQ *af = q->alleles->af;
    int i;
    int numMarkers = isIntervalModel
            ? q->M-1
            : q->M;

    for( i=0; i<numMarkers ; i++)
    {
        if ( ! strcmp( name, af[i].marker_name ) )
        return i;
    }

    return NOT_FOUND;
}

QTL_DATA* validateParams(
        SEXP handle,
        SEXP marker,
        int *locus,
        const int isIntervalModel)
{

    QTL_DATA *q = NULL;
    int id=0;

    *locus = -1;

    /*
     * validate handle
     */
    if ( isInteger(handle) )
    {
        id = INTEGER(handle)[0];
    }
    else if ( isNumeric(handle))
    {
        id = (int)REAL(handle)[0];
    }
    else
    {
        error("attempt to extract locus using non-number handle");
    }

    if ( id >= 0 && id < nqtldata )
    {
        q = qtldata[id];
        if ( q == NULL ) error( "no QTL data");
    }
    else
    {
        error("attempt to extract locus using invalid handle %d" ,id);
    }

    /*
     * validate marker
     */

    if ( isString(marker) )
    {
      const char *string = CHAR(STRING_ELT(marker,0));
      int i = marker_index( string, q, isIntervalModel );
      if (NOT_FOUND == i)
        {
	  error("could not find locus named %s", string);
        }
      *locus = i;
    }
    else if ( isInteger(marker) || isNumeric(marker) )
    {
        int m = isInteger(marker)
                ? INTEGER(marker)[0]
                : (int)REAL(marker)[0];

        int upper = (isIntervalModel)
                ? q->M-1
                : q->M;

        m--; /* to start indexing at 0 */
        if ( m >=0 && m < upper )
        {
            *locus = m;
        }
        else
        {
            error("no such locus %d", m);
        }
    }
    else
    {
        error("locus must be specified as a number or a string");
    }

    return q;
}






SEXP happyprobs ( SEXP handle, SEXP marker ) {

  int locus = -1;
  QTL_DATA *q = validateParams( handle, marker, &locus, 0 );
  SEXP Matrix = R_NilValue;

  if ( locus >= 0 && q->dp_matrices != NULL) {
    ALLELE_FREQ *af = &(q->alleles->af[locus]);
    QTL_PRIOR ***p;
    int S2 = q->S*(q->S+1)/2;

    int i;

  
    /* get the prior probabilities of the strain state combinations at the flanking marker */

    p = allocate_qtl_priors( q );
    compute_qtl_priors( q, p, locus, af->prior );

    PROTECT( Matrix = allocMatrix( REALSXP, q->N, S2) );

    for(i=0; i<q->N; i++) {
      int j, k,m=0;
      for(j=0;j<q->S;j++) {
	for(k=0;k<j;k++,m++) {
	  REAL(Matrix)[i+q->N*m]  = 2*p[i][j][k].prior;
	}
	REAL(Matrix)[i+q->N*m]  = p[i][j][j].prior;
	m++;
      }
    }
    UNPROTECT(1);

    for(i=0;i<q->N;i++) {
      int s1;
      for(s1=0;s1<q->S;s1++)
	free(p[i][s1]);
      free(p[i]);
    }
    free(p);

  }
  return Matrix;
}

SEXP happyprobs2 ( SEXP handle, SEXP marker, SEXP symmetrize ) { 
  /* returns a matrix giving for each row the probability that the corresponding individual is descened form a pair of strains. happyprobs2 differs from happyprobs only in the row of probabilities return is complete - ie due to symmetry it contains the off-diagonal elements twice. this is useful if the probabilities have been computed using parental information to estimate the phase of the haplotypes. */
 

  int locus = -1;
  QTL_DATA *q = validateParams( handle, marker, &locus, 1 );
  SEXP Matrix = R_NilValue;
  int Symmetrize;
  if ( ! isNumeric(symmetrize) || length(symmetrize) != 1 )
    error( "symmetrize is not numeric(1)");
  Symmetrize = (int)REAL(symmetrize)[0]; 

  if ( locus >= 0 && q->dp_matrices != NULL) {
    ALLELE_FREQ *af = &(q->alleles->af[locus]);
    QTL_PRIOR ***p;
    int i;

  
    /* get the prior probabilities of the strain state combinations at the flanking marker */

    p = allocate_qtl_priors( q );
    compute_qtl_priors( q, p, locus, af->prior );

    if ( Symmetrize ) {
      int S2 = q->S*(q->S+1)/2;
      PROTECT( Matrix = allocMatrix( REALSXP, q->N, S2) );

      for(i=0; i<q->N; i++) {
	int j, k,m=0;
	for(j=0;j<q->S;j++) {
	  for(k=0;k<j;k++,m++) {
	    REAL(Matrix)[i+q->N*m]  = p[i][j][k].prior + p[i][k][j].prior;
	  }
	  REAL(Matrix)[i+q->N*m]  = p[i][j][j].prior;
	  m++;
	}
      }
      UNPROTECT(1);
    }
    else {
      int S2 = q->S*q->S;
      PROTECT( Matrix = allocMatrix( REALSXP, q->N, S2) );

      for(i=0; i<q->N; i++) {
	int j, k,m=0;
	for(j=0;j<q->S;j++) {
	  for(k=0;k<q->S;k++,m++) {
	    REAL(Matrix)[i+q->N*m]  = p[i][j][k].prior;
	  }
	}
      }
      UNPROTECT(1);
    }
    for(i=0;i<q->N;i++) {
      int s1;
      for(s1=0;s1<q->S;s1++)
	free(p[i][s1]);
      free(p[i]);
    }
    free(p);
    
  }
  return Matrix;
}

SEXP happynonrecomb ( SEXP handle, SEXP marker ) {

  int locus = -1;
  QTL_DATA *q = validateParams( handle, marker, &locus, 0 );
  SEXP nonrecomb = R_NilValue;

  if ( locus >= 0 ) {
    ALLELE_FREQ *af = &(q->alleles->af[locus]);
    QTL_PRIOR ***p;
    int i;

  
    /* get the prior probabilities of the strain state combinations at the flanking marker */
  
    p = allocate_qtl_priors( q );
    compute_qtl_priors( q, p, locus, af->prior );

    PROTECT( nonrecomb = allocVector( REALSXP, q->N ) );
    for(i=0;i<q->N;i++) {
      REAL(nonrecomb)[i] = q->dp_matrices[i].NonRecomb[locus];
    }
    UNPROTECT(1);
  

    for(i=0;i<q->N;i++) {
      int s1;
      for(s1=0;s1<q->S;s1++)
	free(p[i][s1]);
      free(p[i]);
    }
    free(p);
  }

  return(nonrecomb);
}

SEXP happygenotype ( SEXP handle, SEXP marker ) {

  int locus = -1;
  QTL_DATA *q = validateParams( handle, marker, &locus, 0 );
  SEXP Genotype = R_NilValue;

  if ( locus >= 0 ) {
    ALLELE_FREQ *af = &(q->alleles->af[locus]);
    int i;
    PROTECT( Genotype = allocMatrix( STRSXP, q->N, 2 ) );

    for(i=0;i<q->N;i++) {
      CHROM_PAIR cp =q->genos[i];
      char *g1 = af->allele_name[cp.chrom1[locus]];
      char *g2 = af->allele_name[cp.chrom2[locus]];
      if ( !strcmp( g1, "NA" ) || ! strcmp( g2, "NA") ) {
       SET_STRING_ELT(Genotype,i, R_NaString);
       SET_STRING_ELT(Genotype,i+q->N, R_NaString);
      }
      else {
	SET_STRING_ELT(Genotype,i, mkChar(g1));
	SET_STRING_ELT(Genotype, i+q->N, mkChar(g2));
      }
    }
    UNPROTECT(1);
  }
  return Genotype;
}


SEXP happydesign( SEXP handle, SEXP marker, SEXP model ) {

  SEXP Design = R_NilValue;
  char *mod=NULL;
  int locus = -1;
  QTL_DATA *q = validateParams( handle, marker, &locus, 1 );


  if ( isString(model) ) {
    mod = (char*)CHAR(STRING_ELT(model,0));
  }

  if ( locus >= 0 && q->dp_matrices != NULL ) {
    ALLELE_FREQ *af = &(q->alleles->af[locus]);
    QTL_PRIOR ***p;
    int i;

  
    /* get the prior probabilities of the strain state combinations at the flanking marker */
  
    p = allocate_qtl_priors( q );
    compute_qtl_priors( q, p, locus, af->prior );

    /* allocate and instantiate the design matrix */
  
   
    if ( mod == NULL || ! strcmp( mod, "additive") ) {
      PROTECT( Design = allocMatrix( REALSXP, q->N, q->S ) );
      for(i=0; i<q->N; i++) {
	int j;
	for(j=0;j<q->S;j++) {
	  REAL(Design)[i+q->N*j] = 0.0;
	}
      }
      
      for(i=0; i<q->N; i++) {
	int j, k;
	for(j=0;j<q->S;j++) {
	  for(k=0;k<q->S;k++) {
	    REAL(Design)[i+q->N*j]  += p[i][j][k].prior;
	    REAL(Design)[i+q->N*k]  += p[i][j][k].prior;
	  }
	}
      }

      UNPROTECT(1);
    }
    else if ( mod != NULL && ! strcmp ( mod, "full" ) ) {
      int dim = q->S*(q->S+1)/2;
      PROTECT( Design = allocMatrix( REALSXP, q->N, dim ) );
      
      for(i=0; i<q->N; i++) {
	int j, k, n=0;
	for(j=0;j<q->S;j++) 
	  REAL(Design)[i+q->N*n++]  = p[i][j][j].prior;
	for(j=0;j<q->S;j++) 
	  for(k=0;k<j;k++,n++) {
	    REAL(Design)[i+q->N*n]  = p[i][j][k].prior + p[i][k][j].prior;
	  }
      }
      UNPROTECT(1);      
    }
    else if ( mod != NULL && ! strcmp ( mod, "full.asymmetric" ) ) {
      int dim = q->S*q->S;
      PROTECT( Design = allocMatrix( REALSXP, q->N, dim ) );
      
      for(i=0; i<q->N; i++) {
	int j, k, n=0;
	for(j=0;j<q->S;j++) 
	  for(k=0;k<j;k++,n++) {
	    REAL(Design)[i+q->N*n]  = p[i][j][k].prior;
	  }
      }
      UNPROTECT(1);      
    }
    else {
      warning( "unknown model %s", mod );
    }
  
    for(i=0;i<q->N;i++) {
      int s1;
      for(s1=0;s1<q->S;s1++)
	free(p[i][s1]);
      free(p[i]);
    }
    free(p);
  }
  else {
    /*    warning("Error - locus index %d out of range\n", locus ); */
  }
  
  return Design;
}

double phaseProb( int a1, int a2, int m1, int m2, int p1, int p2, int NA ) {
  double Q12, Q21;
  double T;

  if ( a1 == NA || a2 == NA || m1 == NA || m2 == NA || p1 == NA || p2 == NA )
    return 0.5;
  
  Q12 = (a1==m1)*(a2==p1) + (a1==m2)*(a2==p1) + (a1==m1)*(a2==p2) + (a1==m2)*(a2==p2);
  Q21 = (a2==m1)*(a1==p1) + (a2==m2)*(a1==p1) + (a2==m1)*(a1==p2) + (a2==m2)*(a1==p2);
  
  T = Q12 + Q21;

  if ( T>0 ) {
    Q12 /= T;
  }
  else
    Q12 = 0.5;

  return Q12;
}  
    

/* HAPLOID GENOMES */
SEXP haploid_happydesign( SEXP handle, SEXP marker ) {

  SEXP Design = R_NilValue;
  char *mod=NULL;
  int locus = -1;
  QTL_DATA *q = validateParams( handle, marker, &locus, 1 );

  if ( locus >= 0 && q->haploid_dp_matrices != NULL ) {
    QTL_PRIOR **p;
    int i;

    /* get the prior probabilities of the strain state combinations at the flanking marker */
  
    p = allocate_haploid_qtl_priors( q );
    compute_haploid_qtl_priors( q, p, locus );

    /* allocate and instantiate the design matrix */
  
   
    PROTECT( Design = allocMatrix( REALSXP, q->N, q->S ) );
    for(i=0; i<q->N; i++) {
      int j;
      for(j=0;j<q->S;j++) {
	REAL(Design)[i+q->N*j] = 0.0;
      }
    }
      
    for(i=0; i<q->N; i++) {
      int j;
      for(j=0;j<q->S;j++) {
	REAL(Design)[i+q->N*j]  = p[i][j].prior;
      }
    }

    UNPROTECT(1);
    for(i=0;i<q->N;i++) {
      free(p[i]);
    }
    free(p);
  }

  return Design;
}


void create_haploid_summed_dp_matrices( QTL_DATA *q ) {

  int i;
  double *Pr_ss = q->alleles->Pr_ss;
  double *Pr_st = q->alleles->Pr_st;


  q->haploid_dp_matrices = (HAPLOID_DP_MATRICES*)calloc(q->N,sizeof(HAPLOID_DP_MATRICES));

  for(i=0;i<q->N;i++) {

    if ( i == 0 || genotype_difference( q, i, i-1 ) ) {
      q->haploid_dp_matrices[i].Left = haploid_summed_dp_matrix( q, i, Pr_ss, Pr_st, +1 ); 
      q->haploid_dp_matrices[i].Right = haploid_summed_dp_matrix( q, i, Pr_ss, Pr_st, -1 ); 
      q->haploid_dp_matrices[i].NonRecomb = (double*)calloc( q->M, sizeof(double)); 
    }
    else {
      /*      printf( "same genotype %d %d\n"); */
      q->haploid_dp_matrices[i].Left = q->haploid_dp_matrices[i-1].Left;
      q->haploid_dp_matrices[i].Right = q->haploid_dp_matrices[i-1].Right;
      q->haploid_dp_matrices[i].NonRecomb = q->haploid_dp_matrices[i-1].NonRecomb; 
    }
  }
}


double **haploid_summed_dp_matrix( QTL_DATA *qtl, int individual, double *Pr_ss, double *Pr_st, int direction ) {

  double **X, **a;
  int start, stop, incr;
  CHROM_PAIR *genotypes = &(qtl->genos[individual]);
  ALLELES *A = qtl->alleles;
  int m, s, t;
  int offset;
  double total;
  int markers = genotypes->markers;
  int strains = A->strains;
  double **Pr_tr1;

  Pr_tr1 = (double**)calloc(strains,sizeof(double*));
  for(s=0;s<strains;s++) 
    Pr_tr1[s] = (double*)calloc(strains,sizeof(double));

  X = (double**)calloc(markers,sizeof(double*));
  for(m=0;m<markers;m++) {
    X[m] = (double*)calloc(strains,sizeof(double));
  }

  if ( direction > 0 ) {
    start = 0;
    stop = markers-1;
    incr = +1;
    offset = 0;
  }
  else {
    start = markers-1;
    stop = 0;
    incr = -1;
    offset = -1;
  }

  if ( qtl->an ) 
    a = qtl->an->pr_AtoS[individual][start];
  else
    a = A->af[start].pr_AtoS;
  for(s=0;s<strains;s++) {
    int g1 = genotypes->chrom1[start];
    X[start][s] = a[g1][s];
  }

  for(m=start+incr;m!=stop;m+=incr) { 

    int g1 = genotypes->chrom1[m];
    double pr_ss = Pr_ss[m+offset];
    double pr_st = Pr_st[m+offset];
    double **a;
    if ( qtl->an ) 
      a = qtl->an->pr_AtoS[individual][m];
    else
      a = A->af[m].pr_AtoS;

    for(s=0;s<strains;s++) {
      double norm1 = 1.0e-10;
      for(t=0;t<strains;t++) {
	if ( s == t ) {
	  norm1 += Pr_tr1[s][t] = pr_ss*a[g1][t];
	}
	else {
	  norm1 += Pr_tr1[s][t] = pr_st*a[g1][t];
	}
      }

      for(t=0;t<strains;t++) {
	Pr_tr1[s][t] /= norm1;
      }
    }
    
    for(s=0;s<strains;s++) {
      int S;
      total = 1.0e-10;
      for(S=0;S<strains;S++) {
	X[m][s] += X[m-incr][S] *  Pr_tr1[S][s];
      }
      total += X[m][s];
      /*      for(S=0;S<strains;S++) {
	X[m][s] /= total;
	} */
    }
  }

  for(s=0;s<strains;s++) {
    free(Pr_tr1[s]);
  }
  free(Pr_tr1);

  return X;
}	      


QTL_PRIOR **compute_haploid_qtl_priors(  QTL_DATA *qtl, QTL_PRIOR **qp, int locus  ) {

  /* compute, for each individual, the prior probabilities of the
     different strains at the flanking markers */

  int i;
  int m = locus;
  double total;
  double P[LINKAGE_STATES];
  ALLELES *A = qtl->alleles;
  double d = (A->af[locus+1].position - A->af[locus].position)/100.0;
  double lambda, elambda;
  double MinDist = qtl->alleles->MinDist;

  if ( d < MinDist ) d = MinDist;
  lambda = A->generations*d;
  elambda = exp(-lambda);
  P[BOTH] = elambda;
  P[LEFT] = P[RIGHT] = (1-elambda)/lambda-elambda;
  P[UNLINKED] = 1.0-P[BOTH]-P[LEFT]-P[RIGHT];

  for(i=0;i<qtl->N;i++) {
    int s;
    HAPLOID_DP_MATRICES *dp = &(qtl->haploid_dp_matrices[i]);
    double *Left = dp->Left[m];
    double *Right = dp->Right[m+1];
    dp->NonRecomb[m]=0.0; 

    total = 0.0;
    for(s=0;s<qtl->S;s++) {
      total += qp[i][s].prior = 
	Left[s]*Right[s]*P[BOTH] + Left[s]* P[LEFT] + Right[s] * P[RIGHT]  + P[UNLINKED]; 
    }
    for(s=0;s<qtl->S;s++) 
      qp[i][s].prior /= total;

    dp->NonRecomb[m] /= total;
  }

  return qp;
}

QTL_PRIOR **allocate_haploid_qtl_priors( QTL_DATA *qtl ) {
  
  QTL_PRIOR **qp = (QTL_PRIOR**)calloc(qtl->N,sizeof(QTL_PRIOR*));
  int i;
  
  for(i=0;i<qtl->N;i++) {
    qp[i] = (QTL_PRIOR*)calloc(qtl->S,sizeof(QTL_PRIOR));
  }
  
  return qp;
}

int entrycmp( const void *a, const void *b ) {
  const PARENT_KEY *A = (PARENT_KEY*)a;
  const PARENT_KEY *B = (PARENT_KEY*)b;

  return strcmp( (const char*)A->key, (const char*)B->key );
}

    
