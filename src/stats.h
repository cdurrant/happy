/* misc statistical procedures */

#ifndef _STATS_H_
#define _STATS_H_ 

double rank_lin_regression( double *x, double *y, int from, int to, double *intercept, double *slope, double *sigma, double *t_slope );
double lin_regression( double *x, double *y, int from, int to, double *intercept, double *slope, double *sigma, double *t_slope, double *stderr_slope, double *stderr_intercept );
double *replace_by_ranks( double *array, int start, int stop );
double durbin_watson_test( double *x, double *y, int from, int to, double slope, double intercept );


double perm_test( int N, int M, int **table, int *seed, int permutations );
double chi_stat( int N, int M, int **table, int *margin1, int *margin2, int total );
int **reduce_table( int *N, int *M, int **table );


double erfcc( double x );
double normal_tail( double z );

double betai(double a, double b, double x);
double betacf (double a, double b, double x);
double gammln (double xx);

#endif
