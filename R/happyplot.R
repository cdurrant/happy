                                        # create plots of the result returned by hfit, mfit

happyplot <- function ( fit, mode='logP', labels=NULL,
			xlab=ifelse(together,
				ifelse(missing(chrs),'Marker position',paste('Marker positions from chromosome',paste(chrs,collapse=","))),
				ifelse(missing(chrs),'cM',paste('cM on chromosome',paste(chrs,collapse=",")))),
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
			together=TRUE,		# show multiple chromosomes on single plot
			vlines.chr.lty=vlines.lty, # line width of vertical line identifying chromosome
			vlines.chr.lwd=3*vlines.lwd, # type of vertical line identifying chromosome
			vlines.chr.col=vlines.col, # colour of vertical line identifying chromosome
			vlines.peak.lty=vlines.lty, # line width of vertical line identifying chromosome
			vlines.peak.lwd=2*vlines.lwd, # type of vertical line identifying chromosome
			vlines.peak.col=vlines.col, # colour of vertical line identifying chromosome
			verbose=TRUE,
			... ) {

  #cat("modified version.\n")

  chromosome<-NULL				# single chromosome to work with in this plot
  if (is.null(chrs)) {
	# no chromosomes specifies, print them all
	if (verbose) cat("Not chromosome to plot specified, preparing to print them all.\n")
	chromosome<-unique(fit$chromosome)
  } else {
	chromosome<-chrs
  }
  if (verbose) cat("Chromosome: "); print(chromosome)

  if (!is.null(chromosome)) {
	  if (!together && length(chromosome)>1) {
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

  selected.markers<-TRUE			# take all markers by default
  run.lengths=list(lengths=nrow(fit$table),values="any",cumsum=nrow(fit$table))

  if (!is.null(chromosome)) {
	run.lengths<-rle(fit$chromosome)
	run.lengths$cumsum<-cumsum(run.lengths$lengths)
	if (verbose){cat("run.lengths       :");print(run.lengths)}
	if (verbose){cat("run.lengths$cumsum:");print(run.lengths$cumsum)}
	

	selected.markers <- (fit$chromosome %in% chromosome)
	if (0==length(selected.markers)) {
		warning(paste("No data available for chromosome ",chromosome,".\n",sep=""))
		return(NULL)
	}
	selected.markers[run.lengths$cumsum]<-FALSE # may not be the brightest thing on earth to do

	if (length(chromosome)==1) {
		if (verbose) cat("Plotting fit for chromosome ",chromosome,", selected positions ",
			paste(range(which(selected.markers)),sep="-"),".'\n")
	} else {
		if (verbose) cat("Plotting fit for chromosomes ",paste(chromosome,collapse=",",sep=""),"\n",sep="")
	}
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
  ymax <- mx[2]  # max of value, not excluding artificial extension for labels
  
						# work out how much vertical space to allocate
						# to the marker labels, if present
  if ( ! is.null( labels ) ) {

    if ( is.logical( labels ) & TRUE==labels) {	# allow decent labels with a mere setting to TRUE
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

  cnames = colnames(lp);
  if(verbose) {cat("colnames: "); print(cnames)}
  if (together) {
	  plot.window(xlim=c(1,nrow(lp)),ylim=mx+c(-0.2,0), ...)
  }
  else {
	  rx = range( as.numeric( lp[,1] ) )
	  lx <- rx[2]-rx[1]
	  plot.window(xlim=rx,ylim=mx, ...)
  }

  title(main=main,sub=sub,xlab=xlab,ylab=ylab, ...)
  axis(side=1)
  axis(side=2)
						#
						# PLOTTING
						#

  if (together) {				#       Variant 1: all selected chrs in one graph

	if (!is.null(fit$chromosome)) {
						# plotting separator of chromosomes
		run.lengths.selected<-rle(fit$chromosome[selected.markers])
		run.lengths.selected$cumsum<-cumsum(run.lengths.selected$lengths)
		if (verbose) {
			cat("run.lengths.selected: "); print(run.lengths.selected)
			cat("run.lengths.selected$cumsum: "); print(run.lengths.selected$cumsum)
		}

		data.chromosome<-NULL
		if (length(chromosome)==1) {
			data.chromosome <- rbind(
				start=1,
				end=run.lengths.selected$cumsum-1
			)
		} else {
			data.chromosome <- rbind(
				start=c(1,1+run.lengths.selected$cumsum[1:(length(run.lengths.selected$cumsum)-1)]),
				end=run.lengths.selected$cumsum-1
			)
		}
		if (verbose) {
			cat("data.chromosome: "); print(data.chromosome)
			cat("-run.lengths.selected$values"); print(run.lengths.selected$values)
		}
		colnames(data.chromosome)<-run.lengths.selected$values

		cat("These are the number of markers-pairs for every selected chromosome:\n")
		print(run.lengths.selected)
		cat("and here as a matrix the bits of interest for every chromosome:\n")
		print(data.chromosome)

		for(p in 1:ncol(data.chromosome)) {
			from<-data.chromosome["start",p]
			to<-data.chromosome["end",p]
			chr<-colnames(data.chromosome)[p]
			lines( x=c(to+1,to+1), y=c(0,ymax) , lty=vlines.chr.lty, col=vlines.chr.col, lwd=vlines.chr.lwd)
			text(x=floor(mean(data.chromosome[,p])),y=ymax/22,labels=paste("Chr",chr),col=vlines.chr.col)
						# plotting fragments of data per chromosomes to avoind 'wrong' links between chrs
			for( i in offset:rangemax ) {

				col=lines.col[i-offset+1]

				y.local <- lp[from:to,i]

				if (!is.null(labels)) {
						cat("Plotting lables for chromosomal peaks only.\n")
					y.local.max<-max(y.local)
					y.local.max.which<-from+which(y.local>=y.local.max)-1
					cat("Peaks at: ") ; print(y.local.max.which)
					peak.labels<-fit$table[selected.markers,"marker"][y.local.max.which]
					cat("Peaks labels: ") ; print(peak.labels)
					cat("Peaks height: ") ; print(y.local.max)
					for(max.pos in 1:length(y.local.max)) {
						m=y.local.max.which[max.pos]
						cat("m=",m,", y.local.max=",y.local.max,"\n")
						text(x=m,y=mx[2],labels=peak.labels[max.pos],col=col,srt=labels.srt, ps=labels.ps, adj=0 )
						lines(x=(0.5+c(m,m)), y=c(0,y.local.max) , lty=vlines.peak.lty, col=vlines.peak.col, lwd=vlines.peak.lwd)
					}
				}

				lines(x=from:to,y=y.local, type=type, pch=pch, col=col, lwd=lines.lwd)
			}
		}

	} else {
						# the main data printed all together disregarding all centiMorgan positions
		for( i in offset:rangemax ) {
			lines(x=1:nrow(lp),y=lp[,i], type=type, pch=pch, col=col, lwd=lines.lwd)
		}
	}
						# now printing column names, loop separated to make sure they are readable
	tx <- 1
	wd <- strwidth(cnames)
	buff <- strwidth("spa")
	cat("cnames: "); print(cnames)

	for( i in offset:rangemax ) {
		col=lines.col[i-offset+1]
		cat("cnames[i]: "); print(cnames[i])
		text( x=tx, y=-0.2, labels=cnames[i], adj=c(0),ps=12,col=col)
		tx <- tx + wd[i] + buff[1] 
	}

  } else {					#       Variant 2: one graph per chromosome

						# the labels and connecting lines
						# print early to have them covered by later text
	if ( ! is.null(labels) ) {			
	  	y <- rep( mx[2]*0.99, length(labels$text) )
	  	text( labels$POSITION, y, as.character(labels$text), col=labels.col, srt=labels.srt, ps=labels.ps, adj=0 )
	  	for( m in labels$POSITION) {
			lines( x=c( m,m ), y=c(0,ymax) , lty=vlines.lty, col=vlines.col, lwd=vlines.lwd)
	  	}
	}
						# the main data
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
  }
  par(def.par)

  NULL
}


