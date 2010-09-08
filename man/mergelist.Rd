\name{mergelist}
\alias{mergelist}
\title{Create an object descrbing how to merge strains together}
\description{
 mergelist() is a convenience function which creates a list object
 suitable for use with mergematrices()
}
\usage{
 mergelist( strains, alleles )
}
\arguments{
  \item{strains}{ a character vector of strain names}
  \item{alleles}{a character matrix with one row of strain/allele combinations. There
    must be a named column in the matrix corresponding to every strain
    name in strains. The value of the element is the allele for that strain}
}
\value{
a list of lists of strains describing how the
strains are grouped together. For instance
\code{mergelist <- list( A=list('AJ', 'BALB', 'AKR'),
T=list('RIII','I', 'DBA', 'C57', 'C3H') )}
divides the strains into two groups corresponding to the alleles A,
T (the allele names are not important).
It is essential that the all strain names match all the values in
strains.

The object should be used as an input parameter to mergematrices()
}
\seealso{ mergematrices()}
\author{Richard Mott}
\keyword{models}
