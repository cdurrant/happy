\name{epistasis}
\alias{epistasis}
\alias{epistasispair}
\title{Analysis of Epistasis between Markers}
\description{epistasis() will test for a statistical interaction between
  two sets of
  markers within the happy framework. The markers should be sufficiently far
  apart that they are unlinked (in practice 10cM for a 30 generation HS
  is sufficient). A partial F-test is performed to test if a model
  allowing for interactions fits better than a model in which each
  marker's contribution is additive between loci. Note that the effect
  of each marker within a locus can be either additive or full. Merging
  of strain is permitted.

  epistasispair() is the same as epistasis() except that only one pair
  of markers is tested.
}
\usage{
epistasis( h, markers1, markers2, merge1=NULL, merge2=NULL,
model='additive', verbose=FALSE, family='gaussian' )
epistasispair( h, marker1, marker2, merge1=NULL, merge2=NULL,
model='additive', verbose=FALSE, d1=NULL, d2=NULL, main1=0, main2=0, family='gaussian' )
}
\arguments{
  \item{h}{an object returned by a previous call to happy()}
  \item{markers1}{an array of marker names or indices}
  \item{markers2}{an array of marker names or indices}.
  In epistasis() Every marker
    listed in marker1 is tested for interaction with every marker in
    marker2.If marker2 in NULL then every marker in marker1 is tested
    against every other marker in the array.
  \item{marker1}{ a single marker name or index}
  \item{marker2}{ a single marker name or index}
  In epistasispair() marker1 is tested for interaction with marker2.
  \item{merge1}{an optional merge object (returned by mergematrices())
    determining how the strains should be merged together for the
    markers listed in marker1}
  \item{merge2}{an optional merge object (returned by mergematrices())
    determining how the strains should be merged together for the
    markers listed in marker2}
  \item{model}{the type of model fitted at each locus. Either 'additive'
    or 'full'}
  \item{verbose}{switch controlling output to screen}
  \item{d1}{optional design matrix for the main effect of the first
    marker (saves computation time)}
  \item{main1}{optional log-P-value for the main effect of the first
    marker. NOTE: If d1 is not NULL then main1 \emph{must} be set }
  \item{d2}{optional design matrix for the main effect of the second 
    marker (saves computation time). }
  \item{main2}{optional log-P-value for the main effect of the second
    marker. NOTE: If d2 is not NULL then main2 \emph{must} be set }
  \item{family}{ The distribution of errors in the data. The default is
    'gaussian'. This variable controls the type of model fitting. In the
    gauusian case a standard linear model is fitted using lm(). Otherwise
    the data are fitted as a generalised linear model using glm(), when
    the value of family must be one of the distributions hangled by
    glm(), such as 'binomial', 'gamma'. See family() for the full range
    of models.}
}
\value{
  epistasis() returns a matrix with columns named 'marker1', 'marker2', 'main1', 'main2',
  'main1+main2', 'main1*main2', 'main1.main2'. marker1 and marker2 are
  the names of the markers being compared in a given row, the remaining
  values are the  ANOVA log-P-values of the main
  effects (main1 and main2), the combined additive effect (main1+main2),
  the addtive plus interaction (main1*main2) and the partial F of the
  interaction (main1.main2) after allowing for main1+main2.
  epistasispair() returns a list with the same fields.}
\keyword{models}
\author{Richard Mott}

  
