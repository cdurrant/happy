                                        # create plots of the result returned by hfit, mfit

happyplot <- function ( fit, mode='logP', labels=NULL,
			xlab=ifelse(missing(chrs),'cM',paste('cM on chromosome',chrs)),
			ylab=NULL, main=NULL, sub=NULL, type='s',
			vlines.lty=3,		# to link labels with graph use dotted lines
			vlines.col="lightgray", # and those lines should not distract from the main graph too much
			vlines.lwd=1,		# which is why they are rather thin by default
			labels.col="black",	# colour for labels (from top downwards)
			labels.srt=270,		# degree of rotation for labels
			labels.adj=0,		# adjustment for lables (left-justified)
			labels.ps=8,		# point size of label text
			lines.lwd=2,		# width of lines of main plot
			lines.col=c( "black", "red", "blue", "green", "orange"),
						# colours in which draw the main plots
			pch=20,
			chrs=NULL,		# chromosomes to print
			... ) {

  cat("modified version.\n")

  chromosome<-NULL				# single chromosome to work with in this plot
  if (is.null(chrs)) {
	# no chromosomes specifies, print them all
	cat("Not chromosome to sprint specified, preparing to print them all.\n")
	chromosome<-unique(fit$chromosome)
  } else {
	chromosome<-chrs
	cat("Chromosome: "); print(chromosome)
  }

  if (!is.null(chromosome)) {
	  if (length(chromosome)>1) {
		# we can only deal with a single chromosome per plot, call the function
		# recursively
		for (chr in chromosome) {
			happyplot(fit,mode,labels,xlab,ylab,main,sub,type,
			          vlines.lty,vlines.col,vlines.lwd,
				  labels.col,labels.srt,labels.adj,labels.ps,
				  lines.lwd,lines.col,pch, chrs=chr,  ...)
		}
		rm(chr)
		return(NULL)
	  }
  }

  selected.markers<-T				# take all markers by default
  if (!is.null(chromosome)) {
	selected.markers<-which(fit$chromosome==chromosome)
	if (0==length(selected.markers)) {
		warning(paste("No data available for chromosome ",chromosome,".\n",sep=""))
		return(NULL)
	}
	selected.markers<-selected.markers[1:(length(selected.markers)-1)]
	cat("Plotting fit for chromosome ",chromosome,", selected positions ",
		paste(range(selected.markers),sep="-"),".")
  }


  plot.new(...)					# initialisation of graphics, calling it
						# early to allow use of strwidth
  def.par <- par(no.readonly=TRUE)

  model <- fit$model
  lp <- na.omit(fit$table[selected.markers,])
  test <- fit$test
  offset <- fit$offset
  plots <- fit$width

                                        # title and labels
  if ( ! is.null( fit$permdata ) ) {
    mode <- 'permutation'
    main <- 'permutation plot'
    ylab <- 'permutation logp'
    lp <- na.omit(fit$permdata$permutation.pval[selected.markers,])
  }
  else if ( mode == 'logP' ) {
    ylab <- mode
    if ( is.null(main) )
      main <- 'log probability plot'
  }
  else {
    ylab <- 'SS'
    if ( is.null(main) )
      main <- 'Fitting Sums of Squares plot'
    offset <- offset + plots 
  }

						# the y-axis range
  mx <- 0
  rangemax <- offset + plots-1 
  for( i in offset:rangemax ) {
    r <- range( as.numeric(lp[,i]))
    mx <- range(c( mx, r))
  }
  ymax <- mx[2]
  
						# work out how much vertical space to allocate
						# to the marker labels, if present
  if ( ! is.null( labels ) ) {

    if ( is.logical( labels ) & T==labels) {	# allow decent labels with a mere setting to TRUE
      labels<-list(
	text=fit$table[selected.markers,"marker"],
	POSITION=as.numeric(fit$table[selected.markers,"cM"])
      )
    }

    ps <- par('ps')
    par(ps=8)
    lwidth <- strwidth( as.character(labels$text), units='inches' );
    par(ps=ps)
    lr <- range(lwidth)
    H <- lr[2]*1.2
    
    area <- par( 'fin' ) # width and height in inches of figure
    mai <- par( 'mai' ) # margins in inches 
    lambda <- (area[2]-mai[1]-mai[3]-H)/mx[2] # expansion factor

    h <- H/lambda
    mx[2] <- mx[2] + h
  }

						# preparing window and its dimensions to plot in
  
  if (is.null(lines.col)) {
     warning("The system needs a series of colours in which to draw the probability lines. The 'lines.col' argument should thus not be set to 'NULL'.\n")
     lines.col <- c( "black", "red", "blue", "green", "orange")
  }

  cnames = colnames(lp );
  rx = range( as.numeric( lp[,1] ) )
  lx <- rx[2]-rx[1]
  plot.window(xlim=rx,ylim=mx, ...)

  title(main=main,sub=sub,xlab=xlab,ylab=ylab, ...)
  axis(side=1)
  axis(side=2)

						# the labels and connecting lines
						# print early to have them covered by later text
  if ( ! is.null(labels) ) {			
    y <- rep( mx[2]*0.99, length(labels$text) )
    text( labels$POSITION, y, as.character(labels$text), col=labels.col, srt=labels.srt, ps=labels.ps, adj=0 )
    for( m in labels$POSITION) {
      lines( x=c( m,m ), y=c(0,ymax) , lty=vlines.lty, col=vlines.col, lwd=vlines.lwd)
    }
  }

						# the main data, not using 'plot' to avoid
						# empty pages e.g. when printing to PDF
  if ( rangemax >= offset ) {
    tx <- rx[1] + 0.02*(lx)
    #ty <- c( 0.95*mx[2] ) # to be removed
    wd <- strwidth(cnames)
    buff <- strwidth("spa")
    for( i in offset:rangemax ) {
      col=lines.col[i-offset+1]
      text( tx, 0, cnames[i], adj=c(0),ps=12,col=col)
      tx <- tx + wd[i] + buff[1] 
      lines( x=lp[,1], y=lp[,i], type=type, pch=pch, col=col, lwd=lines.lwd)
    }
  }
  par(def.par)

  NULL
}


