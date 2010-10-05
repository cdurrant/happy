                                        # create plots of the result returned by hfit, mfit

happyplot <- function ( fit, mode='logP', labels=NULL, xlab='cM', ylab=NULL, main=NULL, sub=NULL, type='s',
			vlines.lty=3,		# to link labels with graph use dotted lines
			vlines.col="lightgray", # and those lines should not distract from the main graph too much
			vlines.lwd=1,		# which is why they are rather thin by default
			labels.col="black",	# colour for labels (from top downwards)
			labels.srt=270,		# degree of rotation for labels
			labels.adj=0,		# adjustment for lables (left-justified)
			labels.ps=8,		# point size of label text
			lines.lwd=2,		# width of lines of main plot
			lines.col="black",	# colour in which draw the main plot
			pch=20,
			 ... ) {

  plot.new(...)					# initialisation of graphics, calling it
						# early to allow use of strwidth

  def.par <- par(no.readonly=TRUE)

  model <- fit$model
  lp <- na.omit(fit$table)
  test <- fit$test
  offset <- fit$offset
  plots <- fit$width

                                        # title and labels
  if ( ! is.null( fit$permdata ) ) {
    mode <- 'permutation'
    main <- 'permutation plot'
    ylab <- 'permutation logp'
    lp <- na.omit(fit$permdata$permutation.pval)
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
	text=fit$table[,"marker"],
	POSITION=as.numeric(fit$table[,"cM"])
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

  colours <- c( "black", "red", "blue", "green", "orange")
  cnames = colnames(lp );
  rx = range( as.numeric( lp[,1] ) )
  lx <- rx[2]-rx[1]
  tx <- c( rx[1] + 0.02*(lx) )
  plot.window(xlim=rx,ylim=mx, ...)

  ty <- c( 0.95*mx[2] ) 
  text( tx, ty, cnames[offset], adj=c(0))
  wd <- strwidth(cnames)
  buff <- strwidth("spa")

  title(main=main,sub=sub,xlab=xlab,ylab=ylab, ...)
  axis(side=1)
  axis(side=2)

						# the main data, not using 'plot' to avoid
						# empty pages e.g. when printing to PDF
  lines( x=lp[,1], y=lp[,offset], type=type, pch=pch, col=lines.col, lwd=lines.lwd,...)

  if ( rangemax > offset ) {
    for( i in (offset+1):rangemax ) {
      col=colours[i-offset+1]
      tx <- c( tx[1] + wd[i-1] + buff[1]) 
      text( tx, ty, cnames[i], adj=c(0),ps=12,col=col)
      lines( x=lp[,1], y=lp[,i], type=type, pch=pch,ps=1,col=col)
    }
  }

						# the labels
  if ( ! is.null(labels) ) {
    y <- rep( mx[2]*0.99, length(labels$text) )
    text( labels$POSITION, y, as.character(labels$text), col=labels.col, srt=labels.srt, ps=labels.ps, adj=0 )
    for( m in labels$POSITION) {
      lines( x=c( m,m ), y=c(0,ymax) , lty=vlines.lty, col=vlines.col, lwd=vlines.lwd)
    }
  }

  par(def.par)


  NULL
}


