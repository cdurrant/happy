#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"stats.h"
#include"cmp.h"

/* misc statistical procedures */



double rank_lin_regression( double *x, double *y, int from, int to, double *intercept, double *slope, double *sigma, double *t_slope ) {

  double *rankx = replace_by_ranks( x, from, to );
  double *ranky = replace_by_ranks( y, from, to );
  double c;
  double e_slope, e_intercept;

  c =  lin_regression( rankx, ranky, 0, to-from+1, intercept, slope, sigma, t_slope, &e_slope, &e_intercept);
  free(rankx);
  free(ranky);
  return c;
}

double lin_regression( double *x, double *y, int from, int to, double *intercept, double *slope, double *sigma, double *t_slope, double *stderr_slope, double *stderr_intercept ) {

  double s_x, s_y, ss_x, ss_y, ss_xy;
  int k;
  double N=to-from+1, R;
  
  s_x = s_y = ss_x = ss_y = ss_xy = 0.0;
  for(k=from;k<=to;k++) {
    s_x += x[k];
    ss_x += x[k]*x[k];
    s_y += y[k];
    ss_y += y[k]*y[k];
    ss_xy += y[k]*x[k];
  }  

  s_x /= N;
  s_y /= N;
  ss_x = ss_x-s_x*s_x*N;
  ss_y = ss_y-s_y*s_y*N;
  ss_xy = ss_xy-s_x*s_y*N;

  *slope = ss_xy/ss_x;
  *intercept = s_y - *slope * s_x;

  *sigma = sqrt( (ss_y - *slope *ss_xy)/(N-2) );
  *t_slope = *slope*sqrt(ss_x)/(*sigma);

  *stderr_slope = (*sigma)/sqrt(ss_x);
  *stderr_intercept = *sigma*sqrt((1.0/N+s_x*s_x/ss_x));
  R = ss_xy/sqrt(ss_x*ss_y);
  return R; /* correlation coefficient */
}

double *replace_by_ranks( double *array, int start, int stop ) {

  int len = stop-start+1;
  double *rank = (double*)calloc( len, sizeof(double) );
  double **ptr = (double**)calloc( len, sizeof(double*) );
  int n;

  for(n=0;n<len;n++) {
    rank[n] = array[n+start];
    ptr[n] = &rank[n];
  }

  qsort( ptr, len, sizeof(double*), Fcmp );

  for(n=0;n<len;n++) {
    *ptr[n] = n;
  }

  free(ptr);
  return rank;
}

    

double durbin_watson_test( double *x, double *y, int from, int to, double slope, double intercept ) {

  int len = to-from+1;
  double *diff;
  double sigma=0, dsigma=0, z;
  double dw;
  int n;

  diff = (double*)calloc(len,sizeof(double));

  fprintf(stderr,"here 1 intercept =%.4f slope = %.4f\n", intercept, slope);
  for(n=from;n<=to;n++) {
    diff[n-from] = y[n]-intercept-slope*x[n];
/*    printf("%d %7.3f %7.3f %7.3f\n",n,diff[n-from], x[n], y[n]); */
  }
  fprintf(stderr,"here 2\n");


  sigma = diff[0]*diff[0];
  for(n=1;n<len;n++) {
    sigma += diff[n]*diff[n];
    z = diff[n]-diff[n-1];
    dsigma += z*z;
  }
  fprintf(stderr,"here 3\n");


  dw = dsigma/sigma;
  printf("durbin-watson: sigma = %.4f dsigma = %.4f dw = %.4f\n", sigma, dsigma, dw );

  {
    double p, q, pp, pq, qp, qq, z, e, len1, chi;

    p = q = 0.0;
    for(n=0;n<len;n++) {
      if ( diff[n] > 0 ) 
	p++;
      else 
	q++;
    }
    p /= len;
    q /= len;
    
    pp = pq = qp = qq = 0.0;
    for(n=1;n<len;n++) {
      if ( diff[n] > 0 ) {
	if ( diff[n-1]>0 ) 
	  pp++;
	else
	  pq++;
      }
      else {
	if ( diff[n-1]>0 ) 
	  qp++;
	else
	  qq++;
      }
    }
    
    len1 = len-1;
    
    e = q*q*len1;  z = ( qq - e );  chi = z*z/e;
    e = p*q*len1;  z = ( pq - e );  chi += z*z/e;
    e = q*p*len1;  z = ( qp - e );  chi += z*z/e;
    e = p*p*len1;  z = ( pp - e );  chi += z*z/e;
    printf("chisq = %.3f\n", chi );
  }  


  free(diff);


  return dw;
}



double normal_tail( double z ) {

  return erfcc( z/M_SQRT2 )/2.0;

}

double erfcc( double x ) {
  double t, z, ans;

  z = fabs(x);
  t = 1.0/(1.0+0.5*z);
  ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.096778418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  return x > 0.0 ? ans : 2.0-ans;
}
