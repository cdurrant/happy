\name{hfit}
\alias{hfit}
\alias{hfit.sequential}
\alias{pfit}
\alias{normalise}
\title{Fit a model to an object returned by happy() }
\description{ hfit() fits a QTL model to a happy() object, for a set of markers
  specified. The model can additive or full (ie allowing for dominance
  effects). The test is a partial F-test. In the case of the full model
  two tests are performed: the full against the null, and the full
  against the additive.

  hfit.sequential() performs an automated search for multiple QTL,
  fitting marker intervals in a sequential manner, and testing for a QTL
  conditional upon the presence of previously identified QTL. This is
  essentially forward selection of variables in multiple regression, and
  very similar to composite interval mapping.

  pfit() is a conveneince function to fit several univariate phenotypes
  to the same genotype data.

  normalise() is a convenience function to convert vector of phenotype
  data into a set of standard Gaussian deviates: the values are first
  ranked and then the ranks replaced by the corresponding percentiles in
  a standard Normal distribution. This may be used to help map traits that
  are strongly non-normal (or use the permute argument in hfit()).
  
}
\usage{
hfit( h, markers=NULL, model='additive', mergematrix=NULL,
covariatematrix=NULL, verbose=FALSE, phenotype=NULL, family='gaussian', permute=0 )

hfit.sequential( h, threshold=2,  markers=NULL, model='additive',
mergematrix=NULL, covariatematrix=NULL, verbose=FALSE, family='gaussian')

pfit( h, phen, markers=NULL, model='additive', mergematrix=NULL,
covariatematrix=NULL, verbose=FALSE, family='gaussian' )

normalise(values)
}
\arguments{
  \item{h}{ an object returned by a previous call to happy()}
  \item{markers}{ a vector of marker intervals to test. The markers can
    either be specified by name or by index. Default is NULL, in which
    case all the markers are fitted (same as setting markers=h$markers) }
  \item{model}{ specify the type of model to be fit. Either 'additive',
    where the contrinutions of each allele at the locus are assumed to
    act additively, or 'full', in which a term for every possible
    combination of alleles is included. The default 'additive' mimics
    the behaviour of the original C HAPPY software.}
  \item{mergematrix}{ specify a mergematrix object ( returned by mergematrices() ) which describes which founder
    strains are to be merged. This is used to test whether merging
    strains reduces statistical significance (see mergematrices())}
  \item{covariatematrix}{Optional additional matrix of covariates
    to include in the mode. These may be additional markers (terurned by
    hdesign) or covariates such as sex, age etc. }
  \item{verbose}{ control whether to print the results of the fits
    to the screen, or work silently (the default)}
  \item{threshold}{ the logP threshold used in hfit.sequential to
    decide whether to include a marker interval in the QTL model. The
    default is 2, ie a marker interval must have a partial F statistic
    with P-value <0.01 (=logP 2) to be included.}
  \item{family}{ the distribution of errors in the data. The default is
    'gaussian'. This variable controls the type of model fitting. In the
    gauusian case a standard linear model is fitted using lm(). Otherwise
    the data are fitted as a generalised linear model using glm(), when
    the value of family must be one of the distributions hangled by
    glm(), such as 'binomial', 'gamma'. See family() for the full range
    of models.}
    \item{permute}{The number of permutations to perform. Default is
      0, i.e. no permutation testing is done, and statistical
      significance is assessed by ANOVA (or Analysis of Deviance if
      family != 'gaussian') If permute>0 then statistical significance is
      assessed based on permuting the phenotypes between individuals,
      repeating the model fit, and finding the top-scoring marker
      interval. The emprical distribution of the max logP values is then
      used to assess statistical significance. This technique is useful
      for non-normally distributed phenotypes and for estimating
      region-wide significance levels Note that permutation testing is very slow.}
  \item{phen}{a data.frame containing additional phenotypes. Each column
    should be numeric (only used by pfit())}.
  \item{phenotype}{ An optional  vector containing the phenotype
    values. Used to override the default phenotype in h$phenotype
    (hfit() only ). }
  \item{values}{a numeric vector of phenotype values to transform into
    normal deviates (normalise() only) }
}


\value{hfit() returns a list. The following components of the list are
  of interest:
    \item{table}{a table with the log-P values of the F statistics. The
      table contains rows, one per marker interval. The columns are the
      negative base-10 logarithms of the F-test P-values that there is
      no QTL in the marker interval. In the case of model='full', the
      partial F-test that the full model is no better than the additive
      is also given.

      In the special case of model='additive' and verbose=TRUE the
      effects of all estimable strains are compared with a T-test,
      taking into account the correlations between these
      estimates. However, it should be noted that estimates of
      individual strain effects may be hard to interpret when some
      combinations of strains are indistinguishable, and it is possible
      for the overall F-statistic to be very significant whilst none of
      the strains appear to be significant, based on their
      T-statistics. The F-statistic is a better indicator of the true
      overall fit of the model.  
    }
    \item{permdata}{ a list containing the results of the permutation
      analysis, or NULL if permute=0. The list contains the following
      elements:
	\itemize{
      \item{N}{ The number of permutations }
      \item{permutation.dist}{ A vector containing sorted ANOVA logP values from the N
	permutations. These values can be used to estimate the shape of
	the null distribution, and plotted e.g. using hist().}
      \item{permutation.pval}{A data table containing the permutation
	p-values for each marker interval. The columns in the datatable
	give the position in cM, the marker name (left-hand marker in
	the interval), the original ANOVA logp, the permutation pval for
	this logp, and the log permutation P-value. Bothe global (ie
	region-wide) and pointwise pvalues are given. The Global pvalue for
	a marker interval is the fraction of times that the logP for the
	interval (either additive or full, depending on the model
	specified) is exceeded by the maximum logP in all intervals for
	permuted data. The pointwise pvalue is the fraction of
	permutation logP at the marker interval that exceed the logP for
	that interval.
      }
	}
      }

	 
    The object returned by hfit() is suitable for plotting with
    happyplot()


    pfit() returns a list of hfit() objects, the n'th being the fit for
    the n'th column (phenotype) in phen. 
  }
  \examples{
## An example session:
# initialise happy
\dontrun{h <- happy('Hs.data','HS.alleles')}
# fit all the markers with an additive model
\dontrun{f <- hfit(h)}
# plot the results
\dontrun{happyplot(f)}
# fit a non-additive model
\dontrun{ff <- hfit(h, model='full')}
# view the results
\dontrun{write.table(ff,quote=F)}
# plot the results
\dontrun{happyplot(ff)}
# use noramlised trait values
\dontrun{ff <- hfit(h,phenotype=normalise(h$phenotypes))}
# permutation test with 1000 permutations 
\dontrun{ff <- hfit(h, model='full', permute=1000)}

}
\seealso{ happy{} }
\keyword{models}  
\author{Richard Mott}

