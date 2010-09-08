/* CMP contains standrd comparison functions for use with qsort */

#define _GNU_SOURCE
#include<string.h>
#include<stdio.h>
#include<ctype.h>
#include"cmp.h"

int icmp( const void *A, const void *B) {

  const int *a = (const int*)A;
  const int *b = (const int*)B;

  return *a-*b;
}

int Icmp( const void *A, const void *B) {

  const int **a = (const int**)A;
  const int **b = (const int**)B;

  return **a - **b;
}

int fcmp( const void *A, const void *B) {

  const float *a = (const float*)A;
  const float *b = (const float*)B;

  float x = *a - *b;
  if ( x > 0.0 )
    return 1;
  else if ( x < 0.0 )
    return -1;
  else
    return 0;
}

int Fcmp( const void *A, const void *B) {

  const float **a = (const float**)A;
  const float **b = (const float**)B;

  float x = **a - **b;
  if ( x > 0.0 )
    return 1;
  else if ( x < 0.0 )
    return -1;
  else
    return 0;
}

int Rstrcmp( const void *A, const void *B) {

  const char *a = (const char*)A;
  const char *b = (const char*)B;

  /* string comparison with strings reversed, case insensitive */


  int la = strlen(a)-1;
  int lb = strlen(b)-1;
  int n;
  while ( la && lb )
    {
      if (n = ( ((int)a[la--]) - ((int)b[lb--]) ) ) {
	return n;
      }
    }
  return la-lb;
}

int Strcmp( const void *A, const void *B) {

  const char **a = (const char**)A;
  const char **b = (const char**)B;

  return strcmp(*a,*b);
}


/* a version of strcmp which works with null strings */
	  
int nStrcmp( const void *A, const void *B) {

  const char *a = (const char*)A;
  const char *b = (const char*)B;

  if ( a && b )
    return strcmp(a,b);
  else if ( a )
    return 1;
  else if ( b ) 
    return -1;
  else 
    return 0;
}

int SStrcmp( const void *A, const void *B) {

  const char ***a = (const char***)A;
  const char ***b = (const char***)B;

  return Strcmp((const void*)*a,(const void*)*b);
}

int uscmp( const void *A, const void *B) {

  const unsigned short *a = (const unsigned short*)A;
  const unsigned short *b = (const unsigned short*)B;

  return (int)(*a-*b);
}


int dcmp( const void *A, const void *B) {

  const double *a = (const double*)A;
  const double *b = (const double*)B;

  double x = *a - *b;
  if ( x > 0.0 )
    return 1;
  else if ( x < 0.0 )
    return -1;
  else
    return 0;
}
