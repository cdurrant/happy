                                        # create plots of the result returned by hfit, mfit

happyplot <- function ( fit, mode='logP', labels=NULL, xlab='cM', ylab=NULL, main=NULL, t='s', pch=20, ... ) {

  def.par <- par(no.readonly=TRUE)
  plot.new()

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
  
                                        # work out how much vertical space to allocate to the marker labels, if present
  
  if ( ! is.null( labels ) ) {
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

  
  colours <- c( "black", "red", "blue", "green", "orange")
  cnames = colnames(lp );
  par(col="black")
  par(lwd=2)
  plot( x=lp[,1], y=lp[,offset],  ylim=mx,main=main,xlab=xlab,ylab=ylab, t=t, pch=pch, ...)

  rx = range( as.numeric( lp[,1] ) )
  lx <- rx[2]-rx[1]
  tx <- c( rx[1] + 0.02*(lx) )

  ty <- c( 0.95*mx[2] ) 
  text( tx, ty, cnames[offset], adj=c(0))
  wd <- strwidth(cnames)
  buff <- strwidth("spa")

  if ( rangemax > offset ) 
    for( i in (offset+1):rangemax ) {
      par(col=colours[i-offset+1])
      tx <- c( tx[1] + wd[i-1] + buff[1]) 
      text( tx, ty, cnames[i], adj=c(0))
      par(ps=1)
      lines( x=lp[,1], y=lp[,i], t=t, pch=pch)
      par(ps=12)
    }


  par(col="black")

                                        # the labels
  
  if ( ! is.null(labels) ) {
    par(srt=270)
    par(adj=0)
    par(ps=8)
    y <- rep( mx[2]*0.99, length(labels$text) )
    text( labels$POSITION, y, as.character(labels$text) )
    par(lwd=1)
    par(col='black')
    x <- labels$POSITION
    for( m in x) {
      lines( x=c( m,m ), y=c(0,ymax) )
    }
    par(srt=0)
  }

  par(def.par)
  NULL
}


