\name{cache}
\alias{save.genome}
\alias{save.happy}
\alias{load.genome}
\alias{load.markers}
\alias{the.chromosomes}

\title{Save HAPPY design matrices and genotypes to disk for rapid
  reloading}
\description{ \code{save.genome()} will persist the happy design matrices or
  genotypes from a series of happy objects to disk as a collection of R delayed data packages (as implemented in the
  package \code{g.data}). \code{load.genome()} "reloads" the data, although the
  matrices are not actually loaded  into memory until
  used. \code{load.markers()} loads in a specific set of design matrices or
  genotypes, as defined by their marker names. These functions are very
  usefiul when access to a random selection of loci across the genome is
  required, and when it would be impossible for reasons of space to load many entire HAPPY objects
  into memory. \code{save.happy()}  saves a single
  happy object as a delayed data package. \code{the.chromosomes()} is a
  conveniemce funtion that generates a character vector of chromosome names.}
}
\usage{
save.genome( gdir, sdir, prefix, chrs=NULL, file.format="ped",
mapfile=NULL,ancestryfile=NULL, generations=50, phase="unknown", haploid=FALSE )
genome <- load.genome( sdir, use.X=TRUE, chr=the.chromosomes(use.X=use.X) )
marker.list <- load.markers( genome, markers )
save.happy( h, pkg, dir, model="additive" )
}
\arguments{
  \item{gdir}{ Path to the directory containing the genotype
    (.alleles and either .data  or .ped ) input files required to instantiate happy
    objects. This directory wil1 typically contain a pair of files for
    each chromosome of the genome of interest}. 
  \item{sdir}{Path to the directory where the data will be saved by \code{save.genome}, and
    read back by \code{load.genome()}.}
  \item{prefix}{ Text fragment used to define the file names sought by
    \code{save.genome()}.  An attempt is made to find files in
    \code{gdir} named like \code{chrN.prefix.*} where N is the
    chromosome number (1...20, X, Y), as defined in \code{chrs}.}
  \item{chrs}{ List of chromosome numbers to be processed.}
  \item{use.X}{ logical to determine whether to use X-chromsome data, in
    load.genome().}
  \item{file.format}{Defines the input genotype file format, either
    "ped" (Ped file format) or "happy" ( HAPPY .data file format).}
  \item{mapfile}{ Name of a text file containing the physicla (base
    pair) map for the genome. It contains three columns named "marker",
    "chromosome" and "bp". Every marker in the .alleles files should be
    listed in the file.}
  \item{generations}{ The number of generations since the HS was founded
    (see happy()).}
  \item{genome}{ An object returned by \code{load.genome()}.}
  \item{markers}{ A vector of marker names. These names will be searched
    for in the \code{genome} object, and if found, their corresponding
    data retrieved.}
 \item{haploid}{ A boolean variable indicating if the genomes should be
   interpreted as haploid, ie. homozygous at every locus. This option is
   used for the analysis of both truly haploid genomes and for
   recombinant inbred lines where all genotypes should be
   homozygotes. Note that the format of the genotype file (the .data
   file) is unchanged, but only the first allele of each genotype is
   used in the analysis.The default value for this option is FALSE,
   i.e. the genomes are assumued to be diploid and heterozygous. }
   \item{ancestryfile}{ An optional file name that is used to provide
     subject-specific ancestry information. More Soon...}
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

  \item{h}{ A HAPPY object}
  \item{pkg}{The name of the R delayed data package to be created}
  \item{dir}{ Name of directory to create a delayed data package for a
    single happy object}
  \item{model}{ One of "additive", "full", "genotype"}.
}
\value{
  \code{save.genome()} returns NULL.
  \code{load.genome()} returns a list object which contains information
  about the delayed datapackages loaded, and how the markers are
  distributed between the packages. The list comprises two components,
  named "genome" and "subjects". The former is a datatable with columns
  "marker", "chromosome", "map", "ddp" which acts as a genome-wide
  lookup-table for each marker. The latter lists the subject names
  corresponding to the rows in the design matrices or genotypes. NOTE:
  The software assumes that all the chromosome-specific files used in
  \code{save.genome()} are consistent. i.e. the same subjects in the
  same order occur in each chromosome, and that a marker is only present
  once across the genome.
  \code{load.markers()} returns a list of data (either matrices or genotype
    vectors), each datum being named accoring to the relevant marker
  \code{the.chromosomes()} returns a character vector of chromosome names,
  like \code{ c( "chr1", "chr2" ..., "chrX", "chrY" )}.
}
\seealso{ happy(). Note that the function happy.save() differs from
  save.happy(), in that it saves a single
  happy object for reloading with \code{load()}; it does not use
  delayed data loading.}
\examples{
}
\keyword{models}
\author{Richard Mott}
