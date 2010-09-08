\name{mergeprepare}
\alias{mergeprepare}
\alias{mergefit}
\alias{fastmergefit}
\alias{condmergefit}
\title{ Perform tests to determine whether individual
  polymorphisms could have given rise to a QTL}
\description{

    mergeprepare() reads in datafiles descrbing the locations
  and strain distribution patterns of polymorphisms (SNPs or otherwise) which
  have not necessarily been genotyped. The following tasks are performed:
  \enumerate{
    \item the polymorphism data are read in from testmarkerfile.For each
    polymorphism the corresponding sketon marker interval is determined,
    based on their coordinates. Only
    those polymorphisms lying inside a skeleton marker interval are retained.
    \item the coordinates (typically in bp rather than cM) of the genotyped markers are read in from
    markerposfile. Note that these coordinates are distinct from those
    in the cM map in h\$map used in happy(). Only those markers listed in
    markerposfile that are also in h\$markers are retained - the rest
    are discarded. The retained markers are referred to as 'skeleton'
    markers as they define a framework of genotype data that can us used
    to test the significance of other polymorphisms.
  }

    mergefit() tests each of the polymorphisms to see if it
  could be a QTL. It performs the following operations on each polymorphism:
  \enumerate{
    \item The founder strains are merged together based on the strin
    distribution pattern for that polymorphism.
    \item The merged data are used to fit a QTL in the corresponding
    skeleton marker interval
    \item The unmerged data are used to fir a QTL in the corresponding
    skeleton marker interval.
    \item The fits of the merged and unmerged data are compared with a
    partial F-test. If the unmerged data are significant but the merged
    data are not then there is evidence to reject the polymorphism as
    being associated with the trait.
  }

  fastmergefit() is a convenience function which perfroms a complete
  analysis without making a prior call to happy().

  condmergefit() performs a conditional analysis in which each variant
  is fitted conditional upon every other variant being included in
  turn. This is VERY SLOW.
}
\usage{
mergeprepare( h, markerposfile, testmarkerfile, verbose=FALSE )
mergefit( h, mergedata,  model='additive', covariatematrix=NULL,
verbose=FALSE )
fastmergefit( datafile, allelesfile, markerposfile,
testmarkerfile, generations=200, model='additive', verbose=FALSE )
condmergefit( h, mergedata, model='additive', covariatematrix=NULL,
verbose=FALSE )
}
\arguments{
\item{h}{ an object returned by a previous call to happy()}
\item{markerposfile}{ the name of a text file containing the names and
  locations of the genotyped markers. Contains two names columns
  'marker' and 'POSITION'}
\item{testmarkerfile}{ the name of a text file containign the names,
  positions and strain/allele distribution patterns for each
  polymorphism to be tested. Contains two columns 'marker' and
  'POSITION' plus an additional named column for each of the strains
  listed in h\$strains - \emph{the column names and strain names must
    match exactly}.
}
\item{verbose}{switch to control the level of ouput sent to the
  screen}
\item{mergedata}{ an object created by a previous call to
  mergeprepare() }
\item{model}{ determine the type of model to be fitted -
  either 'additive' or 'full'.

  For the additive model it is assumed that
  the contribution to the phenotype from each chromosome is additive, ie
  if the founder strains at the locus being tested are \eqn{s,t} then
  the expected phenotype will be of the form
  \eqn{T_s+T_t}{T(s)+T(t)}.

  For the full model the expected phenotype
  will be of the form \eqn{T_{st}}{T(s,t)}.

  Analysis of variance is used to test for differences between the
  estimated effects \eqn{T_s, T_{st}}{T(s), T(s,t)}.
  
  The additive model is a submodel of the full, so for model='full' in addition a partial
  F-test is performed to test if the full model explains more variance
  than the additive.
}
\item{covariatematrix}{ an optional design matrix which can be
  used to include additional terms in the model, such as other markers
  (using the matrix returned by hdesign()) and/or other covariates such
  as sex, age etc }
\item{datafile}{the name of a genotype datafile to be passed to happy()}
\item{allelesfile}{the name of the corresponding alleles datafile to be
  passed to happy() }
\item{generations}{ the number of generations to be passed to happy() }
}
\value{
  mergeprepare() returns a list with the following named elements:
  \itemize{
    \item{markerpos}{the positions of the markers}
    \item{interval}{an array. interval[m] contains the index of the
      genotyped marker interval in which the polymorhism p is located,
      or NULL if it is outside all genotyped intervals. 
    }
    \item{markers}{}
    \item{testmarkerdata}{details about the polymorphisms to be tested}
  }

  mergefit() and fastmergefit() return an object, called say 'fit', suitable for plotting using
  mergeplot(). It contains a named element 'table' containing the log-P
  values as in hfit(), which can be printed using
  \code{write.table(fit$table)}.

  condmergefit() returns a table with columns  "position", "interval",
  "sdp", "logPself", "logPmax", "logPmaxPosition" .
  
}
\examples{
## An example session:
# initialise happy
\dontrun{h <- happy('Hs.data','HS.alleles')}
# prepare the merge files
\dontrun{prep <- mergeprepare('markers.positions','testmarkers.txt')}
# run the merge fit
\dontrun{fit <- mergefit( h, prep )}
# alternative, and equivalent, use of fastmergefit():
\dontrun{fit <- fastmergefit( 'Hs.data','HS.alleles',
'markers.positions','testmarkers.txt' )}
# plot the results
\dontrun{mergeplot( fit, prep )}
}
\seealso{ happy(), mergeplot() }
\keyword{models}
\author{Richard Mott}
