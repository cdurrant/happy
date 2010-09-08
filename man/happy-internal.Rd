\name{happy-internal}
\alias{happy-internal}
\title{Internal Happy Functions}
\description{Internal functions for happy. These are not normally called
  by the user}
\alias{comparelist}
\alias{matrixSquared}
\alias{mfit}
\alias{twofit}
\alias{condfit}
\alias{strain.effects}
\alias{glmfit}
\alias{sdp}
\alias{hprob2}
\alias{hnonrecomb}
\alias{h.sum.prob2}

\usage{

matrixSquared( matrix1, matrix2 )
twofit( happy, marker1, marker2, merge1=NULL, merge2=NULL, model =
'additive', verbose=TRUE, family='gaussian' )
mfit( happy, markers, model='additive', mergematrix=NULL,
covariatematrix=NULL, verbose=TRUE , family='gaussian', variants=NULL )
condfit( happy, markers, condmarker, merge=NULL, condmerge=NULL,
model='additive',condmodel='additive', epistasis=FALSE, verbose=TRUE, family='gaussian' )
strain.effects( h, fit, family='gaussian' )
glmfit( formula=NA, family='gaussian' )
sdp( strains, alleles )
}
\author{Richard Mott}
\keyword{internal}
