\name{gfit}
\alias{gaussian.iterate}
\alias{gaussian.null}
\alias{gaussian.loop}
\alias{gaussian.init}
\alias{gaussian.fn}
\alias{gaussian.gr}
\alias{gfit}

\title{Fit a Gaussian Mixture Model to an object returned by happy()}
\description{ gfit() fits a QTL model to a happy() object. The model is a mixture of Gaussians, each with a different mean, and corresponds loosely to the "full" model in hfit(). The difference is that hfit() fits the observed phenotype values to the expected phenotypes under a full model, whereas gfit() uses maximum likelihood to fit the observed phenotype values to a mixture of Gaussians, each with a different mean but common variance. The other functions in this suite are not usually called directly.

The statistical model fitted is as follows.    Consider first the case of
  fitting a QTL at a locus, \eqn{L}. Let \eqn{\bf y}{\bold{y}}  be the
  vector of trait values. Let \eqn{\bf X_L}{\bold{X(L)}} be the
  design matrix for fitting a QTL at the locus \eqn{L}. Let
  \eqn{\bf t_L}{\bold{ t(L)}} be the vector of parameters to be estimated at the
  locus. For an additive QTL, the paramters are the strain effect sizes;
  for a full interaction model there is a paramter for every possible
  strain combination. Then the
  one-QTL model is

  \deqn{\bf E(y) = X_L t_L}{\bold{E(y) = X(L).t}}

    There are \eqn{S(S-1)/2} parameters to
  be estimated in a full model allowing for interactions between the
  alleles within the locus, so the \eqn{i,j}'th element of the design matrix \eqn{\bold{X}} is related to the strain probabilities thus:

  \deqn{\bold{X_{Lij} = F_{iLst}}}{\bold{X(Lij) = F(iLst)}}, where
  \deqn{j(s,t) = \min(s + S(t-1), t + S(s-1))}

  In the function hfit(), the observed phenotypes are regressed directly on the expected trait values. This is not an optimal procedure becuase the data are really a mixture:

  \deqn{ y_i \tilde \sum_{st} F_{iLst} f( ( y_i - \beta_{Lst} )/2\sigma_L^2) }

where \eqn{f(x)} is a standard Gaussian density. The \eqn{\bold{\beta_L}} is a vector of mean trait values for the strain combinations. The parameters \eqn{\bold{\beta_L} , \sigma_L} are estimated by maximum likelihood, and the test for the presence of a QTL at locus \eqn{L} is equivalent to the test that all the \eqn{\beta_{st}=\mu}, when the model collapses to a single Gaussian distribution. 

The model-fitting is implemented in the function gfit() by an iterative process, rather like a simplified version of EM. Is is slower than hfit(), and generally gives similar results as far as overall QTL detection is concered,m but gives more accurate parameter estimates. The log-likelihood for the data is

\deqn{ L = \sum_{i} \log ( \sum_j p_{ij} \frac{\exp(-\frac{(y_i-\beta_j)^2}{2\sigma^2})}{\sqrt{2\pi \sigma^2}})  }


\deqn{ = \sum_i \log ( \sum_j p_{ij} \exp(-\frac{(y_i-\beta_j)^2}{2\sigma^2})) - \frac{N \log(2\pi\sigma^2)}{2}  }

Differentiating wrt to the parameters gives

\deqn{ \frac{\partial L}{\partial \sigma^2} = \sum_i \frac{\sum_j p_{ij} (y_i-\beta_j)^2 \exp(-\frac{(y_i-\beta_j)^2}{2\sigma^2}) }{ 2\sigma^4 \sum_j p_{ij} \exp(-\frac{(y_i-\beta_j)^2 }{2\sigma^2})}  - \frac{N}{2\sigma^2}  }

\deqn{ \frac{\partial L}{\partial \beta_j } = - \sum_i \frac{ p_{ij}   \frac{(y_i-\beta_j) }{\sigma^2} \exp( -\frac{(y_i-\beta_j)^2}{2\sigma^2})}{ \sum_j e_{ij}} }
\deqn{ = \frac{1}{\sigma^2} \left(  - \sum_i \frac{y_i e_{ij} }{\sum_j e_{ij}} + \beta_j \frac{\sum_i e_{ij} }{\sum_j e_{ij}} \right) }

write 
\deqn{ w_{ij} = \frac{p_{ij} \exp(-\frac{(y_i-\beta_j)^2}{2\sigma^2}) }{ \sum_j \frac{p_{ij} \exp(-\frac{(y_i-\beta_j)^2}{2\sigma^2})}}}

then the mle satisfies

\deqn{ \hat{\sigma^2} = \frac{1}{N} \sum_i \frac{\sum_j p_{ij}(y_i-\beta_j)^2 \exp(-\frac{(y_i-\beta_j)^2}{2\sigma^2})}  {\sum_j p_{ij}\exp(-\frac{(y_i-\beta_j)^2}{2\sigma^2})}}





\deqn{ \hat{\sigma^2} = \frac{1}{N} \sum_i \sum_j \hat{w}_{ij}(y_i-\hat{\beta}_j)^2 }




\deqn{ \hat{\beta_j} = \frac{\sum_i  \hat{e}_{ij} y_i}{\sum_i \hat{e}_{ij}} }

and the log-likelihood is

\deqn{ \hat{L} =  \sum_i(\log \sum_j \hat{e}_{ij}) - \frac{N \log(2\pi\hat{\sigma}^2)}{2} }


}

\usage{
gfit( h,eps=1.0e-4, shuffle=FALSE, method="optim" )
gaussian.loop( d, maxit=100, eps=1.0e-3, df=NULL ) 
gaussian.null( n, sigma2 ) 
gaussian.init( d ) 
gaussian.iterate( d, params ) 
gaussian.fn( p, d=NULL )
gaussian.gr( p, d=NULL )

}
\arguments{
	\item{h}{ an object returned by a previous call to happy()}
	\item{shuffle}{ boolean indicating whether the shuffle the phenotypes to perform a permutation test}
	\item{method}{The optimisation algorithm. Default is to use R's "optim" function, which uses derivative information. All other values of this argument will use an EM type iteration.}
	\item{d}{ a list comprising two elements d, probs}
	\item{maxit}{ the maximum number of iterations in the ML fitting}
	\item{eps}{ the terminatation accuracy in the model fitting : the log likelihood must change by less than eps in successive iterations}
	\item{df}{ the degress of freedom to use. If NULL then this is computed as the rank of the data}
	\item{n}{ the number of observations with non-missing phenotypes}
	\item{sigma2} {the variance of the errors in the data}
	\item{params}{ a list with two components, beta = the group means and sigma = the standard deviation}
	\item{p}{ vector of paramters. For internal use only}
}

\value{ gfit() returns a matrix with columns "marker", "LogL", "Null", "Chi", "df", "Pval", "LogPval". Each row of the column describes the fit of the model for thecorresponding  marker interval.


	gaussian.loop() fits the model to a single marker and returns a list with the same elements as in hfit()

	gaussian.iterate() performs a single iteration of the fitting process and returns a list with the updated LogL, beta, sigma, dbeta and dsigma

	gaussian.init() intialises the parameters under the Null model, ie to the case where the means are all identical and the variance is the overal variance.

	gaussian.null() returns the log-likelihood under the Null model

	gaussian.fn() and gaussian.gr() are the function and gradient required by the optim function.
}

\examples{
## An example session:
# initialise happy
\dontrun{h <- happy('HS.data','HS.alleles')}
# fit all the markers
\dontrun{f <- gfit(h)}
}

\seealso{ happy{}, hprob{} }
\keyword{models}  
\author{Richard Mott}


