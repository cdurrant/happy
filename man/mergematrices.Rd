\name{mergematrices}
\alias{mergedpositionmatrix}
\alias{mergematrices}
\title{Construct matrices used to merge together founder strains}
\description{mergematrices() creates a list containing two matrices
  suitable for pre-multiplying with an additive or full  happy marker
  design matrix, in order to produce matrices with certain columns
  combined. These reduced matrices are used to test whether the
  specified merge reduces the significance of the fit.
  This function is not usually called directly but is used by mergfit()
  and hfit()
  megedpositionmatrix() will return either the merged design matrix or
  the mergematrices object corresponding to an object returned by mergeprepare()

 
}
\usage{
 mergematrices( strains, mergelist=NULL, verbose=FALSE )
 mergedpositionmatrix( h, position, prepmerge, model='additive',
verbose=FALSE, design=TRUE )
}
\arguments{
  \item{strains}{character array of strain names}
  \item{mergelist}{ a list of lists of strains describing how the
    strains are grouped together. For instance
    \code{mergelist <- list( A=list('AJ', 'BALB', 'AKR'),
      T=list('RIII','I', 'DBA', 'C57', 'C3H') )}
    divides the strains into two groups corresponding to the alleles A,
    T (the allele names are not important).
    It is essential that the all strain names match all the values in
    strains.}
  \item{verbose}{switch to determine whether to tell what
    is happening.}
  \item{h}{ an object returned by a previous call to happy() }
  \item{position}{ the coordinate of the polymorphism to be tested, ie
an entry in prepmerge$testmarkerdata$POSITION}
  \item{prepmerge}{ an object returned by mergeprepare()}
  \item{model}{ the type of model to be fitted - 'additive'
    or 'full'}
  \item{design}{ switch to make mergepositionmatrix return the
    mergematrix object rather than the merged design matrix}
}
\seealso{happy(), mergefit(), hfit(), mergelist()}
\value{
  mergematrices() and mergepositionmatrix() return an object comprising a list with two elements:
  \itemize{
    \item{amat}{ the matrix to apply to an additive-model design
     matrix}
  \item{imat}{ the matrix to apply to a full-model (interaction) design
    matrix}
}
}
\keyword{models}
\author{Richard Mott}
