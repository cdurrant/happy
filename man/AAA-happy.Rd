\name{Happy}
\alias{introduction}
\alias{happy}
\alias{happy.matrices}
\alias{happy.save}
\title{Quantitative Trait Locus analysis in Heterogeneous
  Stocks}
\description{happy is an \R interface into the HAPPY C package
 for fine-mapping Quantitative Trait Loci (QTL) in mosaci crosses such
 as Heterogenous
 Stocks (HS). HAPPY uses a multipoint analysis which offers significant
 improvements in statistical power to detect QTLs over that achieved by
 single-marker association.  An HS is an advanced intercross between
 (usually
 eight) founder inbred strains of mice. HS are suitable for
 fine-mapping QTL. The happy package is an extension of the
 original C program happy; it uses the C code to compute the
 probability of descent from each of the founders, at each
 locus position, but the happy packager allows a much richer
 range of models to be fit to the data. 

 happy() is used to initialise input files and perform dynamic
 programming in C. Model fitting is then performed by subsequent calls
 to hfit() etc. Input file foramt is described at
 \url{http://www.well.ox.ac.uk/happy}
 }

\details{
\bold{Biological Background}
 
 Most phenotypes of medical importance can be measured quantitatively,
 and in many cases the genetic contribution is substantial, accounting
 for 40\% or more of the phenotypic variance. Considerable efforts have
 been made to isolate the genes responsible for quantitative genetic
 variation in human populations, but with little success, mostly because
 genetic loci contributing to quantitative traits (quantitative trait
 loci, QTL) have only a small effect on the phenotype. Association
 studies have been proposed as the most appropriate method for finding
 the genes that influence complex traits. However, family-based studies
 may not provide the resolution needed for positional cloning, unless
 they are very large, while environmental or genetic differences between
 cases and controls may confound population-based association studies.

These difficulties have led to the study of animal models of human
traits. Studies using experimental crosses between inbred animal strains
have been successful in mapping QTLs with effects on a number of
different phenotypes, including behaviour, but attempts to fine-map QTLs
in animals have often foundered on the discovery that a single QTL of
large effect was in fact due to multiple loci of small effect positioned
within the same chromosomal region. A further potential difficulty with
detecting QTLs between inbred crosses is the significant reduction in
genetic heterogeneity compared to the total genetic variation present in
animal populations: a QTL segregating in the wild need not be present in
the experimental cross.

In an attempt to circumvent the difficulties encountered with inbred
crosses, we have been using a genetically heterogeneous stock (HS) of
mice for which the ancestry is known. The heterogeneous stock was
established from an 8 way cross of C57BL, BALB/c, RIII, AKR, DBA/2, I, A
and C3H/2 inbred strains. Since its foundation 30 years ago, the stock
has been maintained by breeding from 40 pairs and, at the time of this
experiment, was in its 60th generation. Thus each chromosome from an HS
animal is a fine-grained genetic mosaic of the founder strains, with an
average distance between recombinants of 1/60 or 1.7 cM.

Theoretically, the HS offers at least a 30 fold increase in resolution
for QTL mapping compared to an F2 intercross. The high level of
recombination means that fine-mapping is possible using a relatively
small number of animals; for QTLs of small to moderate effect, mapping
to under 0.5 cM is possible with fewer than 2,000 animals. The large
number of founders increases the genetic heterogeneity, and in theory
one can map all QTLs that account for progenitor strain genetic
differences. Potentially, the use of the HS offers a substantial
improvement over current methods for QTL mapping.


\bold{Problem Statement and Requirements}

\enumerate{
  
    \item HAPPY is designed to map QTL in Heterogeneous Stocks (HS), ie populations founded from known inbred lines, which have interbred over many generations. No pedigree information is required.
    \item Obviously, phenotypic values for the trait must be known for all individuals. It is preferable that these are normally distributed because HAPPY uses Analysis of Variance F statistics to test for linkage (however, a permutation test can be used instead).
    \item For each genotyped marker, it is necessary to know the ancestral alleles in the inbred founders (which by definition must be homozygous), and the genotypes from the individuals in the final generation.
    \item The chromosomal position in centiMorgans of each marker must be known.
    \item Missing data are accomodated provided these are due to random failures in the genotyping and not selective genotyping based on the trait values (however, it is permissible to selectively genotype all the markers provided the same individuals are genotyped at each locus). 
  }
  
\bold{What HAPPY does}

HAPPY's analyis is essentially two stage; ancestral haplotype reconstruction using dynamic programming, followed by QTL testing by linear regression:
\itemize{
    \item Assume that at a QTL, a pair of chromosomes originating from
    the progenitor strains, labelled \eqn{s,t} contribute an unknown  amount
    \eqn{T_{st}}{T(st)} to the phenotype. In the special case where the contribution from each chromosome is
    additive at the locus then \eqn{T_{st} = T_s + T_t}{T(s,t) = T(s)+T(t)},say.
    \item a test for a QTL is equivalent to testing for differences between the \eqn{T}'s.
    \item A dynamic-programming algorithm is used to compute the
    probability \eqn{F_{iLst}}{F(n,s,t)} that a given individual \eqn{i} has the ancestral
    alleles \eqn{s, t} at locus labelled \eqn{L}, conditional upon all the genotype
    data for the individual. Then the expected phenotype is
    \deqn{ y = \sum_{st} T_{st} F_{iLst}}{ y = Sum (st) T(s,t)F(i, L,s,t)},
    and the \eqn{T}'s are estimated by a linear regression of the observed phenotypes on these expected values across all individuals, followed by an analysis of variance to test whether the progenitor estimates differ significantly.
    \item The method's power depends on the ability to distinguish
    ancestral haplotypes across the interval; clearly the power will be
    lower if all markers in a region have the same type of
    non-informative allele distribution, but the markers can share
    information where there is a mixture.
    \item   All inference is based
    on regression of the phenotypes on the probabilities
    of descent from the founder loci, \eqn{F_{nst}}{F(n,s,t)}.

  
  }


  Although the models are presented here in the linear model framework (ie least-squares estimation,
  with ANOVA F-tests), it is of course straighforward to extend them to \R's
  generalised linear model framework. Multivariate analysis is also possible.

    It is
  straighforward to fit models involving the effects of multiple loci
  and of covariates. It is easiest to see this by rewriting the problem
  in standard linear modelling notation. Consider first the case of
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
  alleles within the locus, and \eqn{S-1} parameters in an additive
  model. For the full model, the \eqn{i,j}'th element of the design matrix \eqn{\bold{X}} is related to the strain probabilities thus:

  \deqn{\bold{X_{Lij} = F_{iLst}}}{\bold{X(Lij) = F(iLst)}}, where
  \deqn{j(s,t) = \min(s + S(t-1), t + S(s-1)}
  and for the additive model
  \deqn{\bold{X_{Lij} = \sum_s F_{iLsj}}}{\bold{X(Lij) = \sum(s) F(iLst)}} 

  \bold{ More complex models}
  
  To add covariates to the model (for instance sex or  age ) we add
  additional columns \eqn{\bold{C}} to the design matrix:

  \deqn{ \bf E(y) = \left[ X_L \| C \right] ( t_L \| c ) }{\bold{E(y) = [X(L)|C].[t|c]}}

  where \eqn{\bf C} is a design matrix representing the covariates of
  interest, and \eqn{\bf c}{\bold{c}} are the parameters to be
  estimated.\eqn{\bf (t_L\|c)}{\bold{(t(L)|c)}} represents the vector formed by
  adjoining the vectors \eqn{t_L}{t(L)} and \eqn{c}. Note
  that at present \eqn{\bold{C}} must be a numeric matrix: factors must
  be explicitly converted into columns of dummy variables.

  Similarly to fit an additional locus \eqn{K} we adjoin the design
  matrix \eqn{\bold{X_K}}{\bold{X(K)}}, for example:
  
\deqn{ \bf E(y) = \left[ X_L \| X_K \| C \right] ( t_L \| t_K \|c ) }{\bold{E(y) = [X(L)|X(K)|C].[t(L)|t(K)|c]}}

  (this is essentially \bold{composite interval mapping}). The
  happy package allows the inclusion of arbitrary covariate matrices,
  which can include other loci;
  new loci are then tested to see if they significantly improve the fit
  conditional upon the presence of the covariates. In this way we can
  analyse any number of linear combinations of loci and covariates.

  \bold{Epistasis}, or the interaction between loci, is supported as well. At
  present the package can test for interactions between unlinked loci,
  but not linked loci. The test compares the fit between the sum of the
  additive contributions from each locus and the interaction. This is
  accomplished as follows: Let \eqn{X_L, X_K}{X(L),X(K)} be the design
  matrices for the loci \eqn{L,K}. Let \eqn{m_L}{m(L)} be the number of columns in
  \eqn{X_L}{X(L)} Then form a matrix \eqn{X_{LK}}{X(LK)} whose
  \eqn{m_Lm_K}{m(L)m(K)} columns are
  formed by multiplying the elements in each pair of columns in the
  original matrices.

  
  \bold{Merging Strains}

  An important feature of the happy package is the suite of functions to merge
  strains together. The models described above (particularly the full
  interaction models) have the disadvantage
  that the fits sometimes involve a larger number of parameters, with
  many degrees of freedom. This is particularly true for full
  non-additive models and for epistasis. For
  example in an 8-strain HS, 28 df are required to fir a full model for
  a single locus. Large numbers of degrees of freedom have two problems:
  firstly the models may become overspecified, and secondly even if
  there are plenty of degrees of freedom for the residual error, the
  power to detect an effect is diluted.

  A partial solution is to note that since most polymorphisms are
  diallelic (eg SNPs), it makes sense to group the strains according to
  their alleles at some polymorphic locus. This corresponds to operating
  with design matrices in which certain columns are combined by adding
  their corresponding elements together. A diallelic merge reduces the
  number of degrees of freedom dramatically: only 3 df (instead of 28df) are required to
  fit a full model at a locus (and only 1 df instead of 7df for the
  addtive model), and an epistatic interaction between two
  merged loci will involve only 3df (additive) or 8df (full).  

}


}

\references{ Mott R, Talbot CJ, Turri MG, Collins AC, Flint J. 	
	A method for fine mapping quantitative trait loci in outbred animal stocks.
	Proc Natl Acad Sci U S A. 2000 Nov 7;97(23):12649-54. }
      
      \usage{
happy( datafile, allelesfile, generations=200, standardise=FALSE, 
	phase="unknown", file.format="happy", missing.code="NA",
do.dp=TRUE, min.dist=1.0e-5, mapfile=NULL, ancestryfile=NULL, haploid=FALSE )
happy.matrices( h )
happy.save( h, file )

}

\arguments{
  \item{datafile}{ name of the text file containing the genotype and
    phenotype data in HAPPY format}
  \item{allelesfile} { name of the text file containing the
    allele/strain data in HAPPY format}
  \item{generations}{ the number of breeding generations in the HS}
  \item{standardise}{ if TRUE then the phenotype is transformed to have
    mean 0 and variance 1. This does not affect the identification of QTL
    but can make the interpretation of effects easier.}

  \item{phase}{ If phase=="unknown" then the phase of the genotypes is
  unknown and no attempt is made to infer it. If
  phase="estimate" then it is estimated using parental genotype data
  when available. If phase="known" then it is assumed the phase of the
  input genotypes is correct i.e. the first and second alleles in each
  genotype for an individual are on the respectively the first and
  second chromosomes.  Where phase is known this setting
  should increase power, but it will cause erroneous output if it is
  set when the data are unphased. If phase="estimate" then
  file.format="ped" is assumed automatically, because the input data file must be
  in ped-file format in order to specify parental information.  }

   \item{file.format}{ The format of the genotype file. Either "happy"
   (the default) or "ped". "happy" files do not contain any pedigree
   information. They are structured so that one record corresponds to
   an individual. The first two fields are the subject id (unique) and the
   phenotype value. The remaining fields are the N genotypes for the
   subject (where N is the number of markers specified in the alleles
   file), arranged in 2N fields all separated by spaces. 

   The "ped"
   file format is similar except that in place of the two columns "id"
   and "phenotype" in the original "happy" file format there should be
   six columns "family", "id", "mother", "father", "gender",
   "phenotype". Note that the "family" field can be constant so long
   as the id's are all unique. The resultant name of each subject is
   constructed as "family.id", which must be unique.  The genotype
   data then follow as in "happy" format. Note that if
   phase="estimate" then file.format="ped" automatically; it is only
   necessary to use this option when the input file format is "ped"
   but it is desired not to make any use of pedigree information.}

 \item{missing.code}{ The code for a missing allele in the input
   file. Defaults to "NA". Note that old HAPPY files use "ND"}
 \item{do.dp}{ A switch that turns off the dynamic programming part of
   happy. By default dynamic programmig is performed. The only reason to
   turn this off is when only the genotypes are required.}
 \item{min.dist}{ The minimum genetic distance (in centiMorgans) allowed
   between adjacent merkers. Markers positioned closer than min.dist in
   the input file are treated as being min.dist cM apart. This prevents
   problems with markers at the same position, and which HAPPY cannot
   process.}
 \item{mapfile}{ Optional name of text file containing the physical
   base-pair coordinates of the markers (the alleles file only contains
   the genetic map in centiMorgans). The file format has three columns
   named "marker", "chromosome" and "bp". This file is not required
   unless genome cache objects are to be made (see save.genome()).}
 \item{haploid}{ A boolean variable indicating if the genomes should be
   interpreted as haploid, ie. homozygous at every locus. This option is
   used for the analysis of both truly haploid genomes and for
   recombinant inbred lines where all genotypes should be
   homozygotes. Note that the format of the genotype file (the .data
   file) is unchanged, but only the first allele of each genotype is
   used in the analysis.The default value for this option is FALSE,
   i.e. the genomes are assumued to be diploid and heterozygous. }
 \item{ancestryfile}{ An optional file name that is used to provide
   subject-specific ancestry information. More Soon...}.
  \item{h}{ An object of class "happy"}
  \item{file}{ Name of file in which to save data }
}
\value{
  happy() returns an object of type happy, which should be passed onto
  model-fitting functions such as hfit(). A happy object 'h' is a list with
  a number of useful members:
  \item{strains}{ a character vector containing the names of the founder strains}
  \item{markers}{ a character vector containing the names of the markers, in map
    order}
  \item{map}{ a numeric vector containing the map coordinates in centiMorgans
    of the markers}
  \item{subjects}{ a character vector containing the subject names }
  \item{phenotype}{ a numeric vector containing the subject
    phenotypes}
  \item{handle}{a numeric index used internally by the C-code. Do not
    change.}
  \item{matrices}{ a list of matrices used in model fitting (only
    created after a call to happy.matrices())
  \item{use.pedigrees}{ boolean variable indicating whether pedigree information was used to help determine the phase of the genotypes}
  \item{phase.known}{ boolean variable indicating whether or not the phase of the genotypes is assumed to be known}
}

  happy.save() will save a happy object to a file so that it can be
  re-used in a later session with the load() command.
  
  happy.matrices() is not normally called directly - its function is to
  copy all the dynamic-programming matrices created by a call to happy()
  from the underlying C memory space into R objects. The object returned
  is still a happy object, but with an additional component 'matrices'. It can be used in
  exactly the same way as a normal happy object except that the
  underlying C memory is no longer used. When a happy object is saved
  using happy.save() is is first converted by a call to
  happy.matrices(), and when a happy object is reloaded using load() it
  uses matrices stored in R memory. Thus these functions are a useful way
  to save computing time - the dynamic programming step need only be
  performed once, the data can be persisted to disk and then re-used (e.g. to analyse multiple
  phenotypes) at a later date.

  }

\seealso{ hfit(), mergefit(), happyplot() }
\examples{
\dontrun{h <- happy('HS.data', 'HS.alleles', generations=200)}
\dontrun{ happy.save(h,'h.Rdata')}
\dontrun{ load('h.Rdata')}
}
\keyword{models}  
\author{Richard Mott}

