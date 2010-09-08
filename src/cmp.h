/* headers for simple cmp functions, for use with qsort and bsearch */

#ifndef _CMP_H_
#define _CMP_H_

int icmp( const void *, const void *);  /* int   */
int Icmp(  const void *, const void *);   /* int** */
int dcmp( const void *, const void *);  /* double* */
int fcmp( const void *, const void *);  /* float* */
int Fcmp( const void *, const void *);  /* float** */
int Rstrcmp( const void *, const void *); /* Reversed strcmp */
int Strcmp( const void *, const void *); /* for char** */
int SStrcmp( const void *, const void *); /* for char *** !!!! */
int uscmp( const void *, const void *); /* unsigned short */

#endif
