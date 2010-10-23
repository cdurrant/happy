\name{happyplot}
\alias{happyplot}
\alias{mergeplot}
\title{Plotting functions for happy model fits}
\description{ happyplot() will plot along the genome the log P-value that a QTL is not found in a series of marker intervals. It accepts as input the results of hfit(), mfit() and mergefit(). mergeplot() is a convenience function for calling happlyplot() after a call to mergefit(), with several parameters set. }
\usage{ 
happyplot( fit, mode='logP', labels=NULL, xlab='cM', ylab=NULL, main=NULL, sub=NULL, type='s', pch=20, 
           vlines.lty=3, vlines.col="lightgray", vlines.lwd=1, labels.col="black", labels.srt=270,
	   labels.adj=0, labels.ps=8, lines.lwd=2, lines.col=c( "black", "red", "blue", "green", "orange"),
	   chrs=NULL, ... ) 
mergeplot( fit, mergedata, mode='logP', xlab='bp', ylab=NULL, main=NULL,  t='p', pch=20, ... ) 
}
\arguments{ 
\item{fit}{ an object returned by a previous call to hfit(), mfit(), or mergefit() }
\item{mode}{ the mode of the plot - either 'logP', when the negative base-10 logarithm of the ANOVA P-value of plotted, or 'SS', when the fitting sums-of-squares is plotted}
\item{labels}{ optional matrix detailing marker labels to be drawn on the plot. The labels are written vertically above the plot, with vertical lines extending down into the plot area. labels is a matrix with two named columns 'marker', containing the marker names, and 'POSITION', containing the x-axis positions of the markers.}  
\item{mergedata}{ (mergeplot() only). an object returned by a previous call to mergeprepare(). This is used to construct labels for plotting} 
The following options are passed to the plotting routines: 
\item{xlab}{ the x-axis label }
\item{ylab}{ the y-axis label }
\item{main}{ the titke of the plot }
\item{t}{ the type of plot - either 'p', 'l', 's' or 'S', with the same meanings as in plot() }
\item{pch}{ the plotting character code, with the same meaning as in
  plot() }
\item{vlines.lty}{ type of line to link labels with X axis, set to 0 for blank line, the default draws dotted lines }
\item{vlines.col}{ color of line linking labels with X axis, set to light grey by default }
\item{vlines.lwd}{ width of line that is linking labels with the X axis }
\item{labels.col}{ colour for labels }
\item{labels.srt}{ degree of rotation for label, default orients them from top downwardss }
\item{labels.adj}{ adjustment for lables (left-justified) }
\item{labels.ps}{ point size of label text }
\item{lines.lwd}{ width of lines of main plot }
\item{lines.col}{ colours in which draw the main plots, which should be an array of colours when multiple plots are expected for a single diagram }
\item{chrs}{ chromosomes to print diagrams for. The default is to print them all. }
\item{...}{ arguments passed down to R plotting functions}
}
\value{ A plot to the current graphics device is produced. For
  happyplot, if fit$permdata is not NULL (i.e. hfit() was run using the
  permutation test option) then the plot contains two step functions, labelled
  global.logp and point.logp.  The global.logp plot shows the empirical
  log-pvalue relative the whole region (i.e. adjusted for the number of
  markers) while point.logp shows the empirical log-pvalue for each
  interval. If fit$permdata is NULL then the plots give the ANOVA logP
  values. If the model used in  hfit() is 'additive' then the logP for
  the additive model vs the null model is plotted; if the model is
  'full' then the curves for the full, additive and partial F-test logP
  values are plotted.  }
\examples{
\dontrun{ h <- happy( 'HS.data', 'HS.alleles' ) }
\dontrun{ fit <- hfit( h, h$markers, model='full' ) }
\dontrun{ happyplot( fit ) }
}

\seealso{ hfit(), mfit(), mergefit(), par() }
\keyword{aplot}
\author{Richard Mott}
