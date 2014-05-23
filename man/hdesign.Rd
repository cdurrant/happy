\name{hdesign}
\alias{hdesign}
\alias{hprob}
\alias{hgenotype}
\title{Extract design matrix or genotypes for a specific marker interval from a happy
  object}
\description{hdesign() will call C to extract the design matrix to fit a
  QTL to a marker interval. hprob() will call C to extract a raw probability matrix. hgenotype() will return the raw genotype data for a marker}
\usage{
hdesign( h, marker, model='additive', mergematrix=NULL )
hprob( h, marker=NULL )
hgenotype ( h, marker, collapse=FALSE, sep="" )
}
\arguments{
  \item{h}{ an object returned by a previous call to happy()}
  \item{marker}{ either a character string giving the name of the marker
    or the index of the marker in the array h\$markers}
  \item{model}{either 'additive' (default) or 'full'. The additive
    design matrix returns an array with S columns, where S is the number
    of founder strains in the HS. The full design matrix returns a
    matrix with S(S-1)/2 columns, one for each combination of strains}
  \item{mergematrix}{an object returned by mergematrices, used to define
    sets of strains that are to be merged together. This is accomplished
    by adding the corresponding columns in the original design matrix.}
  \item{collapse}{ a boolean variable indicating whether to collapse the
    alleles into a single genotype.}
  \item{sep}{ the text to be used to separate the alleles if collapsed.}
}
\value{ hdesign() returns a design matrix \eqn{d_{ij}}{d[i,j]}, in which the \eqn{i}th row corresponds to
  the subject \eqn{i}, and the \eqn{j}th column to the corresponding strain or
  combination of strains or merged strains.
 hprob() returns a matrix \eqn{p_{ix}{p[i,x]}}, in which the \eqn{i} the row corresponds to the subject \eqn{i}, and the \eqn{x=s*S+t} th column contains the probability that the ancestral strains are \eqn{s,t} where \eqn{S} is the total number of strains.
 hgenotype() returns a \eqn{Nx2} matrix \eqn{g_{ix}{g[i,x]}} in which
 the \eqn{i} th row corresponds to the subject \eqn{i}, and column 1
 contains the first allele and column 2 the second allele at the marker
 specified, or (if \code{collapse=TRUE}) a vector of genotypes with the
 alleles pasted together.
}
\seealso{ happy(), hfit() }
\examples{
\dontrun{h <- happy( 'HS.data', 'HS.alleles', generations=200 )}
\dontrun{d <- hdesign( h, 1 ) ## the first marker interval}
\dontrun{d <- hdesign( h, 'D1MIT264') ## the marker interval with left-hand marker D1MIT264}
\dontrun{d <- hdesign( h, 'D1MIT264', model='full' ) ## ditto with full design matrix}
}

\keyword{models}
\author{Richard Mott}
