 #library(happy.hbrem)

# simulation code

uniroot.func <- function(effect, allele, noise, effect.size) {
  p <- ifelse(allele,-effect,effect) + noise
  a <- anova( lm(p ~ allele))
  return( a[1,2]/sum(a[,2])-effect.size)
}

adjust.effect.size <- function( allele, noise, effect.size ) {

  N <- length(allele)
  maf <- sum(allele)/N
  if ( maf > 0.5 ) maf = 1.0-maf

  if ( effect.size > 0 ) {
    effect0 <- sqrt(effect.size/(2*maf*(1-maf)))
    interval=c(0.1*effect0,10*effect0)
    while ( uniroot.func( interval[1],allele,noise,effect.size) * uniroot.func( interval[2],allele,noise,effect.size) > 0 ) {
      interval[1] = interval[1]*0.1
      interval[2] = interval[2]*10
    }
  
    u <- uniroot( uniroot.func,interval,allele, noise, effect.size)
    return(u$root)
  }
  else {
    return(0)
  }
}

#simulate

simulate.MAGIC.phenotypes <- function( genome.cache, maf=0.22, effect.size=0.05, n.sim=1, chrom=5) {

  cat("simulating n.sim=", n.sim, " effect.size=", effect.size, "\n")
  genome <- genome.cache$additive$genome


  nchr <- 1
  chrom.end <- rep(0,len=chrom)
  chrom.start <- rep(0,len=chrom)
  chrom.start[1] <- genome$bp[nchr]
  prev.chr <- genome$chromosome[nchr]

  for( i in 2:nrow(genome)) {
#    cat(i, genome$chromosome[i], prev.chr, "\n")
    if ( genome$chromosome[i] != prev.chr ) {
      chrom.end[nchr] <- genome$bp[i-1]
      prev.chr <- genome$chromosome[i]
      nchr <- nchr+1
      chrom.start[nchr] <- genome$bp[i]
    }
  }
  chrom.end[nchr] <- genome$bp[nrow(genome)]
  chrom.length <- rep(0, nchr)
  cum.length <- rep(0, nchr)
  total.len <- 0

  
  for ( j in 1:length(chrom.start)) {
    chrom.length[j] <- chrom.end[j]-chrom.start[j]+1;
    cum.length[j] <- total.len + chrom.length[j]
    total.len <- cum.length[j] 
  }


  simulated.phenotypes <- matrix( NA, nrow=length(g$subjects), ncol=n.sim)
  rownames(simulated.phenotypes) <- g$subjects
  simulated.locations <- data.frame( true.chr=rep(NA,n.sim), true.bp=rep(NA,n.sim), true.locus=rep(NA,n.sim), true.snp=rep(NA,n.sim), true.allele.var=rep(0,n.sim), true.allele.logP=rep(0,n.sim), true.happy.var=rep(0,n.sim), true.happy.logP=rep(0,n.sim), merged.happy.var=rep(0,n.sim), merged.happy.logP=rep(0,n.sim) ) #, best.merged.logP=rep(0,n.sim), best.merged.var=rep(0,n.sim))
  simulated.locations$true.bp <- as.numeric(simulated.locations$true.bp)  

#  browser()
  n <- 1
  while( n <= n.sim ) {
    chr <- NA
    pos <- floor(total.len*runif( 1 ) )
    chr <- min(which(pos<=cum.length))
    if ( chr > 1 )
      pos.bp <- pos - cum.length[chr-1]
    else
      pos.bp <- pos
    
    locus.1 <- max(which(genome$chromosome == chr & genome$bp <= pos.bp))
    locus.2 <- min(which(genome$chromosome == chr & genome$bp >= pos.bp))
    n.strains = length(g$strains)
    n.maf = floor(maf*n.strains)
    if ( n.maf == 0 ) n.maf = 1
		
    if ( locus.1 == locus.2 -1 ) {
      sdp <- sample( 1:n.strains )[1:n.maf]
      high.allele <- 1:n.strains %in% sdp
      locus <- locus.1
      
      d <- hdesign(g,locus)
      samp <- apply( d, 1, function(x) { sample(1:n.strains,size=1,prob=x) } )
      allele <- high.allele[samp]
      maf <- sum(allele)/length(allele)
      if ( maf > 0.5 ) maf = 1.0-maf
      effect <- sqrt(effect.size/(2*maf*(1-maf)))
      noise <- rnorm( nrow(d) )
      effect <- adjust.effect.size( allele, noise, effect.size )

      test <- simulated.phenotypes[,n]<-  noise + ifelse(allele,-effect,effect)
      if (!(n%% 100 )) cat( n, " " )

      
      a <- anova(lm(test~allele))
      true.allele.var <- a[1,2]/sum(a[,2])
      true.allele.logP <- -log10(a[1,5])
      
      a <- anova(lm(test~d))
      true.happy.var <- a[1,2]/sum(a[,2])
      true.happy.logP <- -log10(a[1,5])

      merged <- apply(d[,high.allele],1,sum)

      a <- anova(lm(test~merged))
      merged.happy.var <- a[1,2]/sum(a[,2])
      merged.happy.logP <- -log10(a[1,5])

#      b <- best.merge( test, d )
                      
#      simulated.locations[n,] <- c( chr, pos.bp, locus, genome$marker[locus], true.allele.var, true.allele.logP, true.happy.var, true.happy.logP, merged.happy.var, merged.happy.logP, b$best.logP,b$best.var )
      simulated.locations[n,] <- c( chr, pos.bp, locus, genome$marker[locus], true.allele.var, true.allele.logP, true.happy.var, true.happy.logP, merged.happy.var, merged.happy.logP )
      n <- n+1
    }
  }
  cat("\n")

  simulated.locations$true.bp <- as.numeric(simulated.locations$true.bp)  
  return( list( simulated.phenotypes=simulated.phenotypes, simulated.locations=simulated.locations ))
}


SS <- function(x, mean.correct=T)
{
  if (mean.correct)
    {
      x <- x - mean(x)
    }
  x %*% x
}

# use qr decomposition to perform genome scans on a set of simulated data

quick.scan <- function( g, simulated ) {

  cat("scanning genome...\n")
  genome <- g$additive$genome
  phenotypes <- simulated$simulated.phenotypes
  test <- phenotypes[,1]
  N <- nrow(phenotypes)
  nsim <- ncol(phenotypes)
  logp <- matrix( 0, nrow=nrow(genome), ncol=ncol(phenotypes))

  if ( 1 ) {
    for(j in 1:nrow( genome ) ) {
      cat( j, genome$marker[j], "\n" )
      d <- hdesign(g,j)
      f <- lm( test ~ d )
      qr <- f$qr
      logp[j,] <- apply( phenotypes, 2, anova.logP,qr, N )
    }
  }

  logp <- -logp/log(10)
  w <- apply(logp, 2, which.max )
  loc <- sim$simulated.locations
  compare <- data.frame(  max.logP=rep(0,nsim), max.var=rep(0,nsim), max.chr=rep(NA,nsim), max.locus=rep(0,nsim), max.snp=rep(NA,nsim), max.bp=rep(0,nsim), diff.bp=rep(0,nsim))
  for( i in 1:nsim ) {
    lp <- logp[w[i],i]
    chr <- genome$chromosome[w[i]]
    bp <- genome$bp[w[i]]
    diff.bp <- -1
    if ( identical(as.character(chr),as.character(loc$true.chr[i]) ) )
      diff.bp <-  abs(bp-loc$true.bp[i])
    compare$max.logP[i] <- lp
    test <- phenotypes[,i]
    d <- hdesign(g,w[i])
    f <- lm( test ~ d )
    qr <- f$qr
    compare$max.var[i] <- anova.logP( test, qr, N, get.var=TRUE )
    compare$max.chr[i] <- chr
    compare$max.locus[i] <- w[i]
    compare$max.bp[i] <- bp
    compare$diff.bp[i] <- diff.bp
    compare$max.snp[i] <- genome$marker[w[i]]
    
  }

  compare <- cbind( loc, compare )
  return(list( logP=logp, compare=compare))
}


anova.logP <- function( p, qr, N, get.var=FALSE ) {
  tss <- SS(p)
  rss <- SS(qr.resid(qr,p))
  fss <- tss-rss
  effect <- fss/tss
  df <- qr$rank-1
  F <- (fss/df)/(rss/(N-df))
  if ( get.var ) 
    return ( fss/tss )
  else
    return( pf( F, df, N-df, lower=FALSE, log=TRUE))
}

best.merge<- function( phen, d ) {
  fit <- lm( phen ~ d-1 )
  idx <- order(fit$coef)
  x <- d[,idx[1]]
  csum <- t(apply(d[,idx],1,cumsum))
  splits = length(idx)-1
  pval <- rep(0,splits)
  gvar <- rep(0,splits)
  for ( j in 1:splits) {
    f <- lm( phen ~ csum[,j] )
    a <- anova(f)
    pval[j] <- a[1,5]
    gvar[j] <- a[1,2]/sum(a[,2])
  }

  best <- which.min(pval)
  f <- lm( phen ~ csum[,best] )
  a <- anova(f,fit)
  partial.pval <- a[2,6]
  
  merge <- idx[1:best]
  merge.strains = colnames(d)[merge]
  r <- list(logP=-log10(pval), best.logP=-log10(pval[best]),best.var=gvar[best], partial.logP=-log10(partial.pval), merge=merge, strain.names=merge.strains)
  print(r)
  return( r)
}
  
resample.func <- function( X, idx, phen.resampled ) {
        X.resampled <- X[idx,]
        a <- anova(lm( phen.resampled ~ X.resampled ))
        return( a[1,5])
      }

hdesign2 <- function( snp, g ) { hdesign( g, snp ) }

resample <- function( g, sim, window=4.0e6, resamples=100, sample.fraction=0.8, replace=FALSE ) {
  genome <- g$additive$genome
  locus <- sim$simulated.locations
  n.subjects <- length(g$subjects)
  sample.size <- ceiling(n.subjects*sample.fraction)
  
  n.sim <- nrow(locus)
  results <- list()
  for( i in 1:n.sim ) {
    loc <- locus$true.locus[i]
    loci.idx <- genome$chr==locus$true.chr[i] & abs( genome$bp-locus$true.bp[i] ) < window
    loci <- genome[ loci.idx, ]
    mat <- list()
    phen <- sim$simulated.phenotypes[,i]
    mat <- lapply( loci$marker, hdesign2, g )
    names(mat) <- loci$marker
    n.snp = length(names(mat))
    logP <- matrix(0,nrow=n.snp,ncol=resamples)
    rownames(logP) <- names(mat)
    logP <- replicate( resamples, { idx <- sample(1:n.subjects, size=sample.size, replace=replace);  phen.resampled <- phen[idx]; sapply( mat, resample.func,idx, phen.resampled, simplify=TRUE )}, simplify=TRUE )
    true.logP <- sapply( mat, resample.func,1:n.subjects, phen, simplify=TRUE )

    logP <- -log10(logP)
    w <- apply( logP, 2, which.max )
    freq <- tabulate(w)
    rmip <- rep(0,length(mat))
    rmip[1:length(freq)] <- freq
    names(rmip) <- names(mat)
    rmip <- rmip/resamples
    print(locus$true.snp[i])
    print(rmip)
    results[[i]] <- list( snp=locus$true.snp[i], rmip=rmip, true.logP=true.logP, loci=loci )
  }
  return(results)
}

resample.wrapper <- function( g, n.sim=1, resamples=100, effect.size=0.1, file=NULL) {

  if ( ! is.null( file ) ) {
    sim <- simulate.MAGIC.phenotypes(g, effect.size=effect.size, n.sim=n.sim)
    res <-  res <- resample(g, sim, resamples=resamples)
    save( res, file=file ) 
  }
}

resample.stats <- function( q=0.5) {
  files <- list.files( ".", pattern="^resample" )
  res.all <- list()
  for (f in files) {
    load(f)
    cat(f, "\n")
    res.all <- c( res.all, res )
  }

  N <- length(res.all)
  df <- data.frame( chr=rep(0,N), bp=rep(0,N), snp=as.character(rep(NA,N)), logp.bp=rep(0,N), rmip.bp=rep(0,N), low.q=rep(1,N), high.q=rep(1,N), peak.logP=rep(0,N))
  df$snp <- as.character(df$snp)
  i<-1

  for( r in res.all ) {
    idx <- match( r$snp, r$loci$marker )
    df$bp[i] <- r$loci$bp[idx]
    df$logp.bp[i] <- r$loci$bp[which.min(r$true.logP)]
    df$rmip.bp[i] <- r$loci$bp[which.max(r$rmip)]
    df$chr[i] <- r$loci$chromosome[idx]
    df$snp[i] <- r$snp
    cum <- cumsum(r$rmip)
    df$low.q[i] <- min(r$loci$bp[cum>q/2])
    df$high.q[i] <- max(r$loci$bp[cum<1-q/2])
    df$peak.logP[i] = -log10(min(r$true.logP))
    i <- i+1
  }

  df$logp.diff <- abs(df$bp-df$logp.bp)
  df$rmip.diff <- abs(df$bp-df$rmip.bp)

  return(df);
}

background.distribution <- function( q, sim, g ) {
  chr <- as.numeric(sim$simulated.locations$true.chr)
  genome <- g$additive$genome
  n.sim <- nrow(sim$simulated.locations)
  different.chr <- matrix( FALSE, nrow=nrow(genome), ncol=5)
  for( ch in 1:5 ) 
    different.chr[,ch] <- genome$chromosome != ch

  mx <- rep(0,n.sim)
  for( i in 1:n.sim)
    mx[i] <- max( q$logP[different.chr[,chr[i]],i] )


  return(mx);
  
}

mapping.resolution <- function( sim, mapping.threshold=3.0e6, logp.threshold=3.0 ) {


  good <- good.sims( sim, mapping.threshold=3.0e6, logp.threshold=3.0 ) 
  sim <- sim(good,  mapping.threshold, logp.threshold)
                    
  sim$logp.diff <- sim$max.logP - sim$true.happy.logP
  plot(sim$logp.diff, sim$true.happy.logP)
  return(sim)
}


good.sims <- function(sim, mapping.threshold=3.0e6, logp.threshold=3.0 ) {

  centromere.1 <- 1.0e6* c( 14.15,14.90 )
  centromere.2 <- 1.0e6* c( 3.24, 3.60 )
  centromere.3 <- 1.0e6* c(13.56, 14.05)
  centromere.4 <- 1.0e6* c(3.04, 3.90)
  centromere.5 <- 1.0e6* c( 11.05,11.89)

  centromeric.1 <- sim$true.chr == 1 & sim$true.bp > centromere.1[1] & sim$true.bp < centromere.1[2]
  centromeric.2 <- sim$true.chr == 2 & sim$true.bp > centromere.2[1] & sim$true.bp < centromere.2[2]
  centromeric.3 <- sim$true.chr == 3 & sim$true.bp > centromere.3[1] & sim$true.bp < centromere.3[2]
  centromeric.4 <- sim$true.chr == 4 & sim$true.bp > centromere.4[1] & sim$true.bp < centromere.4[2]
  centromeric.5 <- sim$true.chr == 5 & sim$true.bp > centromere.5[1] & sim$true.bp < centromere.5[2]

  good <- sim$diff.bp != -1 & sim$diff.bp < mapping.threshold & sim$max.logP > logp.threshold & ! ( centromeric.1 | centromeric.2 | centromeric.3 | centromeric.4 | centromeric.5 )
  
  return(good)
}

conf.int <- function( q, g, logP.drop, mapping.threshold=3.0e6 ) {

  logP = q$logP
  sim = q$compare
  which.max.logP <- apply( logP, 2, which.max )
  max.logP <- apply( logP, 2, max )
  logP.threshold <- max.logP-logP.drop
  logP.include <- apply( logP, 1, function(X,thresh) { X > thresh  } , logP.threshold )
  genome <- g$additive$genome

  midpt <- midpoints(g)
  
  
  pos <- cbind(as.numeric(sim$true.chr),as.numeric(sim$true.bp))
  qtl.interval <- t(apply( pos, 1, function( X, genome) { genome$chr == X[1] & abs( genome$bp - X[2]) < mapping.threshold }, genome ))
  logP.include <- logP.include & qtl.interval
  w <- apply(logP.include, 1, which )
  
  ci.left <- sapply(w, min)
  ci.right <- sapply(w, max )
  ci.width <- midpt[ci.right] - genome$bp[ci.left]
  return( data.frame( left=ci.left, right=ci.right, width=ci.width )) 
}


mid.points <- function( g ) {

  genome <- g$additive$genome

  nr <- nrow(genome)
  incr <- 1:nr
  incr = incr+1
  incr[nr] = nr
  chr.end = genome$chromosome[1:nr] != genome$chromosome[incr]
  incr[chr.end] = incr[chr.end]-1
  midpt <- 0.5*(genome$bp[1:nr]+genome$bp[incr])
  return(midpt)
}

estimate.ecdf <- function( sim, logP=NA, logP.window=0.25, mapping.threshold=3.0e6, centromere.rm=TRUE ) {


  good <- sim$diff.bp != -1 & sim$diff.bp < mapping.threshold & abs(sim$max.logP -logP ) < logP.window 

  if ( centromere.rm == TRUE ) {
    centromere.1 <- 1.0e6* c( 14.15,14.90 )
    centromere.2 <- 1.0e6* c( 3.24, 3.60 )
    centromere.3 <- 1.0e6* c(13.56, 14.05)
    centromere.4 <- 1.0e6* c(3.04, 3.90)
    centromere.5 <- 1.0e6* c( 11.05,11.89)
    
    centromeric.1 <- sim$true.chr == 1 & sim$true.bp > centromere.1[1] & sim$true.bp < centromere.1[2]
    centromeric.2 <- sim$true.chr == 2 & sim$true.bp > centromere.2[1] & sim$true.bp < centromere.2[2]
    centromeric.3 <- sim$true.chr == 3 & sim$true.bp > centromere.3[1] & sim$true.bp < centromere.3[2]
    centromeric.4 <- sim$true.chr == 4 & sim$true.bp > centromere.4[1] & sim$true.bp < centromere.4[2]
    centromeric.5 <- sim$true.chr == 5 & sim$true.bp > centromere.5[1] & sim$true.bp < centromere.5[2]

    good <- good  & ! ( centromeric.1 | centromeric.2 | centromeric.3 | centromeric.4 | centromeric.5 )
  }

  return( ecdf( sim$diff.bp[good] ) )
}



sample.estimates <- function( g, phenotype, marker, n.sample=100 ) {
  
  genome <- g$additive$genome
  cc <- complete.cases(phenotype)
  phenotype <- phenotype[cc]
  use.subjects <- g$subjects[match(names(phenotype), g$subjects, nomatch=0)]
  p.subjects <- match(use.subjects,names(phenotype),nomatch=0)
  phenotype <- phenotype[p.subjects]
  g.subjects <- match(use.subjects,g$subjects)

  d <- hdesign(g,marker)
  d <- d[g.subjects,]
  f <- lm( phenotype ~ d )
  a <- anova(f)
  
  X <- apply( d, 1, function(x) { sample(colnames(d),size=n.sample,prob=x,replace=TRUE) } )
  estimates <- apply( X, 1, lm.func, phenotype )
  logp <- estimates[1,]
  estimates <- estimates[2:nrow(estimates),]
  estimates[2:nrow(estimates),] <- estimates[1,] + estimates[2:nrow(estimates),]
  
  return(list(logp=logp, estimates=estimates, happy.logp=-log10(a[1,5]), happy.coef=coef(f)))
}

# imputation code


imputed.one.way.anova <- function( g, phenotype, marker, n.sample=100, model="linear" ) {

  genome <- g$additive$genome
  cc <- complete.cases(phenotype)
  phenotype <- phenotype[cc]
  use.subjects <- g$subjects[match(names(phenotype), g$subjects, nomatch=0)]
  p.subjects <- match(use.subjects,names(phenotype),nomatch=0)
  phenotype <- phenotype[p.subjects]
  g.subjects <- match(use.subjects,g$subjects)
  
  cat( marker, "\n" )
  d <- hdesign(g,marker)
  d <- d[g.subjects,]
  nc = ncol(d)
  mat = matrix( 0, nrow=nc, ncol=nc)
  diag(mat) <- rep(1,nc)
  X <- apply( d, 1, function(x, n.sample) { sample(nc,size=n.sample, prob=x,replace=TRUE) }, n.sample)
  if ( model=="linear") {
    results = as.matrix(t(apply( X, 1, function( x, mat, phenotype ) {  one.way.anova( mat[x,], phenotype) }, mat, phenotype )))
    colnames(results) <- c( "logP", "sigma2.hat", "df", "N.df", "p",  paste("N.", g$strains, sep=""), g$strains )
    return(results)
  }
  else if ( model == "binary") {
    results = as.matrix(t(apply( X, 1, function( x, mat, phenotype ) {  one.way.binary( mat[x,], phenotype) }, mat, phenotype )))
    colnames(results) <- c( "logP", "pi.null", "df", "N.df", "p",  paste("N.", g$strains, sep=""), g$strains )
    return(results)
  }
  else {
    warning( "unknown model ", model, "\n")
    return (NULL)
  }
}



imputed.parameter.distributions <- function( results, model="linear", delta=0.01, fact=4 ) {
  mean.results = apply(results,2,mean, na.rm=TRUE)
  results = data.frame(results)
  p <- results$p[1]
  n <- results$N.df[1]
  N.strains = (ncol(results)-5)/2
  p6 = p+6
  p12 = p+12
  N <- as.matrix(results[,seq(6,6+N.strains-1,1)])
  N <- ifelse( N>0, N, NA)
  est <- results[,seq(6+N.strains,ncol(results),1)]
  logP = mean(results$logP)

  if ( model == "linear" ) {
    sigma2.hat <- results$sigma2.hat %o% rep(1,p)
    stderr <- sqrt(sigma2.hat/N)
  
    x.range = seq( -fact, +fact, delta )
    all.results <- NULL
    print(N)
    print(est)
    all.ok = sum(is.na(N))
                                        #  if (all.ok>0 )
                                        #    browser()
    
    for( i in 1:p ) {
      ok = !is.na(N[,i])
      Nok = sum(ok)
      
      if ( Nok > 0 ) {
        mu.mean = mean(est[ok,i]) 
        sigma2.mean = mean(sigma2.hat [ok,i])
        Ni <- N[ok,i]
        N.mean = mean(Ni)
        se = sqrt(sigma2.mean/N.mean)
        x.target = mu.mean + seq( -4, +4, 0.001 ) *se
        s = (x.target %o% rep(1,Nok) - rep(1,length(x.target)) %o% est[ok,i])/(rep(1,length(x.target)) %o% sqrt(sigma2.hat[ok,i]/Ni) )
        
        t.distribution <- pt( s, n )
        avg.t.distribution = apply(t.distribution, 1, mean )
        distribution <- data.frame( x=x.target, d=avg.t.distribution)
                                        #      plot(distribution$x, distribution$d)
        
        density <- diff(distribution$d)
        x.density = distribution$x[1:length(density)]
        
        mean.i <- sum( x.density*density)
        var.i = sum( x.density*x.density*density) - mean.i*mean.i
        pct.50 =  distribution$x[c(max(which( distribution$d < 0.25 )), min(which( distribution$d > 0.75 )))]
        pct.90 =  distribution$x[c( max(which( distribution$d < 0.05 )), min(which( distribution$d > 0.95 )))]
        
        cat( Nok, names(est)[i], mean.i, sqrt(var.i), pct.50, pct.90, "\n")
        all.results = rbind( all.results, c( names(est)[i], logP, mean.i, sqrt(var.i), pct.50, pct.90))
      }
      else
        all.results = rbind( all.results, c( names(est)[i], logP, NA, NA, NA, NA, NA, NA ))
      
    }
    all.results = data.frame(all.results)
    names(all.results) = c( "Accession", "logP", "mean", "stderr", "lower.50", "upper.50", "lower.90", "upper.90")
    return(all.results)
  }
  else if ( model == "binary" ) {
    all.results <- NULL
    for( i in 1:p ) {
      ok = !is.na(N[,i])
      Nok = sum(ok)
      if ( Nok > 0 ) {

        n.max = max(N[ok,i])
        binom.dist = apply( cbind( est[ok,i], N[ok,i]), 1, function( x, n ) { d =rep(0, n); d[1:(x[2]+1)] = dbinom( 0:x[2], x[2], x[1]); return(d) }, n.max+1   )
        total.dist = apply(binom.dist, 1, mean)
        print(total.dist)
        z=seq(0,n.max); cs = cumsum(total.dist); q05 = max(1,which(cs<=0.05))-1; q25 = max(1,which(cs<=0.25))-1; q50=max(1,which(cs<=0.50))-1; q75 =min(which(cs>=0.75))-1; q95=min(which(cs>=0.95))-1;
        print(cs)
        mu = sum(z*total.dist)/n.max
        se = sqrt(sum(z*z*total.dist)/(n.max*n.max)-mu*mu)
        cat( names(est)[i], logP, mu, se, q25, q75, q05,q95  , "\n")
        all.results = rbind( all.results, c( logP, names(est)[i],  mu, se, q25/n.max, q75/n.max, q05/n.max,q95/n.max  ))
      }
      else {
        all.results = rbind( all.results, c( logP, names(est)[i], logP, NA, NA, NA, NA, NA, NA ))
      }
    }
    all.results = data.frame(all.results)
    names(all.results) = c(  "impute.logP", "Accession", "mean", "stderr", "lower.50", "upper.50", "lower.90", "upper.90")
    return(all.results)
  }
  
}

one.way.anova <- function( X, y ) {

# assumes design matrix X correponds to one-way anove - ie exatly one non-zero element per row, equal to 1, and assumes all y's are non-missing

  n = length(y)
  
  N <- apply(X, 2, sum ) # number of replicates of each class
  df = sum(N>0)
  S <- t(X) %*% y  # sums of replicates
  beta.hat = ifelse( N>0 , S/N, 0 ) # estimates
  y.hat = X %*% beta.hat
  y.resid = y - y.hat
  RSS = sum(y.resid*y.resid)
  sigma.hat = RSS/(n-df)
  yy = y-mean(y)
  TSS = sum(yy*yy)
  FSS = TSS-RSS
  F = 0
  logP = NA
  if (df > 1) {
    F = (FSS/(df-1))/(RSS/(n-df))
    logP = -pf( F, df-1, n-df, lower.tail =FALSE, log=TRUE)/log(10)
  }
#  a = anova(lm(y ~ X))
 
  return( c( logP, sigma.hat, df-1, n-df, ncol(X), N, beta.hat ) )
}

one.way.binary <- function( X, y ) {

# assumes design matrix X correponds to one-way design - ie exatly one non-zero element per row, equal to 1, and assumes all y's are non-missing and either 0 or 1 or TRUE or FALSE
# estimates are given as proportions pi
  
  n = length(y)
  
  N <- apply(X, 2, sum ) # number of replicates of each class
  df = sum(N>0)
  S <- t(X) %*% y  # sums of replicates
  pi.hat = ifelse( N>0 , S/N, NA ) # estimates
  pi.null = sum(y)/n
  y.hat = X %*% pi.hat

  if ( df > 1 ) {
    log.LR = 2*sum( ifelse( y > 0, y*log(y.hat/pi.null), 0) + ifelse( 1-y>0, (1-y)*log((1-y.hat)/(1-pi.null)),0))
    logP = -pchisq( log.LR, df-1, log=TRUE, lower.tail=FALSE )/log(10)
  }
  else{
    logP = NA
    log.LR = 0
  }
  return( c( logP, pi.null, df-1, n-df, ncol(X), N, pi.hat ) )
}

read.phenotypes <- function( file="/raid/ARABIDOPSIS/GENOTYPES/11092008/QTL_MAPPING_IGNORE/fitted.ignore.30112008.txt" ) {
  hsril <- read.table(file, sep="\t", h=T)
  ignore <- c(  "HSRIL-137","HSRIL-191","HSRIL-199","HSRIL-20","HSRIL-214","HSRIL-287","HSRIL-292","HSRIL-344","HSRIL-41", "HSRIL-435","HSRIL-457","HSRIL-465", "HSRIL-483","HSRIL-522","HSRIL-526","HSRIL-544","HSRIL-641","HSRIL-65", "HSRIL-248","HSRIL-605","HSRIL-91", "HSRIL-28", "HSRIL-215","HSRIL-409", "HSRIL-401","HSRIL-719","HSRIL-5",  "HSRIL-648","HSRIL-547","HSRIL-67", "HSRIL-486","HSRIL-552","HSRIL-550","HSRIL-632","HSRIL-662","HSRIL-74" )

  print(ignore)
  print(dim(hsril))
  hsril <- hsril[ !(as.character(paste("HSRIL-",hsril$HSRIL,sep="")) %in% ignore),]
  print(dim(hsril))

  hs <-  unique(sort(as.character(as.factor(paste("HSRIL-",hsril$HSRIL,sep="")))))
  hsril$HSRIL2 <- as.factor(paste("HSRIL-",hsril$HSRIL,sep="")) 
  hsril$HSRIL <- as.factor(paste("-",hsril$HSRIL,sep="")) 
  
  return(hsril)
}

imputed.qtl.estimates <- function( g, qtl.file="/raid/ARABIDOPSIS/GENOTYPES/11092008/QTL_MAPPING_IGNORE/qtl.summary.csv", phenotype.dir="/raid/ARABIDOPSIS/GENOTYPES/11092008/QTL_MAPPING_IGNORE/" , menu.file = "/raid/ARABIDOPSIS/GENOTYPES/11092008/QTL_MAPPING_IGNORE/modelmenu.txt", n.sample=100) {

  qtls = read.delim(qtl.file, sep=",")
  menu = read.delim(menu.file)
  for( i in 1:nrow(qtls) ) {
    phen = as.character(qtls$phenotype[i])
    print(phen)
    phenotype.file = menu$PhenotypeFile[match( phen,menu$Phenotype)]
    model = menu$Function[match( phen,menu$Phenotype)]
    if ( ! is.null( phenotype.file )) {
      phenotype.file = paste( phenotype.dir, phenotype.file, sep="/")
      pheno = read.phenotypes( phenotype.file );
      y = pheno[[phen]]
      if ( ! is.null( y ) ) {
        snp = as.character(qtls$central.snp[i])
        d = hdesign(g,snp)
        names(y) = pheno$SUBJECT.NAME
        results = imputed.one.way.anova( g, y, snp, n.sample=n.sample, model=model ) 
        p = data.frame(imputed.parameter.distributions ( results,  model=model, delta=0.01, fact=4 ))
        np = nrow(p)
        parameters = data.frame( Phenotype=rep(phen,np), SNP = rep(snp, np),chr=rep(as.character(qtls$central.chr[i]),np), bp.from=rep(as.character(qtls$central.from.bp[i]),np), bp.to = rep(as.character(qtls$central.to.bp[i]),np) ,RMIP = rep(as.character(qtls$total.RMIP[i]),np), anova.logP = rep(as.character(qtls$logP[i]),np))
        
        parameters = cbind(parameters,p)
        file=paste(phen, snp, "imputed.txt", sep=".")
        print(file)
        write.table( parameters, file=file, quote=F, sep="\t", row=FALSE)
      }
    }
  }
}

qtl.confidence.intervals <- function( g, qtl.file="/raid/ARABIDOPSIS/GENOTYPES/11092008/QTL_MAPPING_IGNORE/qtl.summary.csv", mapping.error.database.file="/raid/ARABIDOPSIS/GENOTYPES/11092008/SIMULATIONS/mapping.error.database.RData", logP.window=0.25, q=c( 0.5, 0.9, 0.95)) {
  
  qtls = read.delim(qtl.file, sep=",")
  load(mapping.error.database.file)
  ci.all = NULL
  for( i in 1:nrow(qtls) ) {
    phen = as.character(qtls$phenotype[i])
    print(phen)
    snp=as.character(qtls$central.snp[i])
    logP = qtls$logP[i]
    rmip = qtls$total.RMIP[i]
    ci = confidence.interval( phen, snp, g, logP, rmip, mapping.error.database, logP.window=logP.window, q=q)
    ci.all = rbind(ci.all, ci)  
  }
  
  return(ci.all)
}

trait.loci <- function( ci.all, g, file.name="trait.locus.02022009.csv", q=0.90, rmip.threshold=0.25 ) {

  ci.all = ci.all[ci.all$q==q & ci.all$rmip >=rmip.threshold,]

  genome=g$additive$genome
  closest = NULL
  nr = nrow(ci.all)
  for( i in 1: nr) {
    x = ci.all[i,]
    d1 = genome$marker[which.min(ifelse( x$chr == genome$chromosome, abs(x$ci.lower-genome$bp), NA))];
    d2 = genome$marker[which.min(ifelse( x$chr == genome$chromosome, abs(x$ci.upper-genome$bp), NA))]
    print(d1, x$chr, x$ci.lower)
    closest = rbind(closest,c(d1,d2))
  }
  
  
  df = data.frame(
    name=paste( ci.all$phenotype, ci.all$snp, ci.all$chr, sep="-"),
    population=rep("MAGIC-30112008",nr),
    "genome_scan" = ci.all$phenotype,
    "subscan_label" = rep("additive",nr),
    phenotype = ci.all$phenotype,
    marker1 = closest[,1],
    marker2 = closest[,2],
    species="Arabidopsis thaliana",
    chromosome=ci.all$chr,
    "start_bp"=ci.all$ci.lower,
    "end_bp"=ci.all$ci.upper,
    threshold=rep("",nr),
    score=ci.all$rmip,
    peak = rep("",nr),
    label=ci.all$rmip,
    comment = rep("",nr),
    url = rep("",nr)
  )

  write.table(df, file=file.name, quote=FALSE, row=FALSE, sep=",")
}

confidence.interval <- function( phenotype, snp, g, logP, rmip, mapping.error.database, logP.window=0.25, q=c( 0.50, 0.90, 0.95 ) ) {
  
  subset = mapping.error.database[abs(mapping.error.database$logP-logP)<=logP.window,] 
  while( nrow(subset) < 2000 ) {
    logP.window = logP.window+0.05
    subset = mapping.error.database[abs(mapping.error.database$logP-logP)<=logP.window,] 
  }
  
  cdf = ecdf( subset$mapping.error )
  interval = c( min(subset$mapping.error), max(subset$mapping.error))
  w = sapply( q, function(qq,cdf,interval) { u = uniroot( function(x, cdf, qq) { cdf(x)-qq } , interval, cdf, qq); return(u$root) }, cdf, interval )
  genome = g$additive$genome
  idx = match( snp, genome $marker )
  if ( idx == nrow(genome) | genome$chr[idx] != genome$chr[idx+1])
    midpt = 0.5*(genome$bp[idx]+genome$bp[idx+1])
  else
    midpt = genome$bp[idx]
  
  ci.lower = floor(midpt - w)
  ci.lower = pmax(0, ci.lower)
  ci.upper = floor(midpt + w)
  
  chr = genome$chr[idx]
  lq = length(q)
  
  df = data.frame( phenotype=rep( phenotype, lq), snp=rep( snp, lq), logP=rep(logP, lq), rmip=rep(rmip,lq), chr=rep( chr, lq),midpt=midpt,q=q, ci.lower=ci.lower, ci.upper=ci.upper,  width=floor(2*w+1) )
  return( df)
         
}
  
good.sims <- function(sim, mapping.threshold=3.0e6, logp.threshold=3.0 ) {

  centromere.1 <- 1.0e6* c( 14.15,14.90 )
  centromere.2 <- 1.0e6* c( 3.24, 3.60 )
  centromere.3 <- 1.0e6* c(13.56, 14.05)
  centromere.4 <- 1.0e6* c(3.04, 3.90)
  centromere.5 <- 1.0e6* c( 11.05,11.89)

  centromeric.1 <- sim$true.chr == 1 & sim$true.bp > centromere.1[1] & sim$true.bp < centromere.1[2]
  centromeric.2 <- sim$true.chr == 2 & sim$true.bp > centromere.2[1] & sim$true.bp < centromere.2[2]
  centromeric.3 <- sim$true.chr == 3 & sim$true.bp > centromere.3[1] & sim$true.bp < centromere.3[2]
  centromeric.4 <- sim$true.chr == 4 & sim$true.bp > centromere.4[1] & sim$true.bp < centromere.4[2]
  centromeric.5 <- sim$true.chr == 5 & sim$true.bp > centromere.5[1] & sim$true.bp < centromere.5[2]

  good <- sim$mapping.error != -1 & sim$mapping.error < mapping.threshold & sim$max.logP > logp.threshold & ! ( centromeric.1 | centromeric.2 | centromeric.3 | centromeric.4 | centromeric.5 )
  
  return(good)
}

prepare.ci.database <- function( g, files = c("simulated-0.075-1.RData"  ,"simulated-0.1-1.RData"	  ,"simulated-0.15-1.RData"  ,"simulated-0.1-6.RData"    ,"simulated-0.1-7.RData" ,"simulated-0.05-1.RData" ,"simulated-0.125-1.RData"  ,"simulated-0.1-5.RData"   ,"simulated-0.175-1.RData"  ,"simulated-0.20-1.RData","simulated-0.25-1.RData", "simulated-0.30-1.RData"), outfile = "mapping.error.database.RData"  ) {

  ci.database = NULL
  for( f in files ) {
    load(f)
    print(f)
    ci.database = rbind(ci.database, q$compare)
  }
  ci.database$true.locus = as.integer(ci.database$true.locus)
  ci.database$max.locus = as.integer(ci.database$max.locus)
  print(dim(ci.database))
  
  genome = g$additive$genome
  s = 1:nrow(genome)
  n = length(s)
  s1 = s+1
  s1[n] = n
  same.chr = genome$chr[s] == genome$chr[s1]
  midpt = ifelse( same.chr, floor(0.5*(genome$bp[s] + genome$bp[s1])), genome$bp[s])
  ci.database$mapping.error = ifelse( ci.database$true.chr==ci.database$max.chr , abs(midpt[ci.database$true.locus]-midpt[ci.database$max.locus]), -1 )

  good = good.sims( ci.database)
  mapping.error.database <- data.frame( logP = ci.database$max.logP[good], mapping.error=ci.database$mapping.error[good])
  mapping.error.database = mapping.error.database[order(mapping.error.database$logP),]
  save( mapping.error.database, file=outfile )
  return(outfile)
}


  
    
    
  
  
midpt <- function(g) {

  genome = g$additive$genome
  snp = genome$marker
  idx = match( snp, genome $marker )
  idx1 = idx+1
  idx1[nrow(genome)] = idx[nrow(genome)]
  idx1 = ifelse(genome$chromosome[idx] == genome$chromosome[idx1], idx1, idx )
  midpt = 0.5*(genome$bp[idx]+genome$bp[idx1])
  names(midpt) = snp;
  return(data.frame(snp=snp, chr=genome$chromosome, midpt=midpt))
}


# cdf code

qtl.error <- function( g, simdata, error.threshold=-1, pdffile="power.pdf" ) {


  genome <- g$additive$genome
  nc <- nrow(genome)
  cum <- rep(0,length=nc)
  x <- 0
  chr <- genome$chromosome
  bp <- genome$bp
  boundary <- c()
  cum[1] <- bp[1]
  for( i in 2:length(chr) ) {
    if ( chr[i] != chr[i-1] ) {
      x = cum[i-1]
      boundary <- c( boundary, x )
    }
    cum[i] = x + bp[i]
  }

  boundary <- boundary/1.0e6
  cum <- cum/1.0e6

  centromere <- c( 14.15,14.90, 3.24+boundary[1], 3.60+boundary[1], 13.56+boundary[2], 14.05+boundary[2], 3.04+boundary[3], 3.90+boundary[3], 11.05+boundary[4],11.89+boundary[4])

  xleft <- c(  14.15,3.24+boundary[1],13.56+boundary[2],3.04+boundary[3],11.05+boundary[4])
  ybottom <- c( -1, -1, -1, -1,-1)
  xright <- c(14.90,3.60+boundary[1], 14.05+boundary[2], 3.90+boundary[3],  11.89+boundary[4])
  ytop <- c( 2,2,2,2,2)


  m <- match( simdata$true.snp, genome$marker )
  simdata$cum = cum[m]
  simdata$bin = floor(cum[m])
  n.bin <- max(simdata$bin)
  if ( error.threshold == -1 ) 
    simdata.good <- simdata[simdata$diff.bp >= 0,]
  else 
    simdata.good <- simdata[simdata$diff.bp >= 0 & simdata$diff.bp < error.threshold,]
  mean.error <- aggregate( simdata.good$diff.bp, list( bin=simdata.good$bin ), mean )
  overall.median <- median( simdata.good$diff.bp )/1.0e6
  median.error <- aggregate( simdata.good$diff.bp, list( bin=simdata.good$bin ), median )
  total <- tabulate( simdata$bin+1 )
  if ( error.threshold == -1 ) 
    power <- aggregate( simdata$diff.bp, list( bin=simdata$bin), function( X, error.threshold ) { sum(X > error.threshold) }, error.threshold=error.threshold )
  else
    power <- aggregate( simdata$diff.bp, list( bin=simdata$bin), function( X, error.threshold ) { sum(X >= 0 & X < error.threshold) }, error.threshold=error.threshold )
    
  prob <- power$x/total
  print(median(prob))
  binned <- data.frame( bin=seq(0,n.bin, 1 ), total=total, power.count=power$x, power=prob, mean.error = mean.error$x, median.error=median.error$x )

  if ( ! is.null(pdffile)) 
    pdf(pdffile)


  par(mfrow=c(4,1), mar=c(2,5,1,2), oma=c(5,1,2,2))

  plot(binned$bin, binned$power,ylim=c(0,1), xlim=c(0,130), ylab="power", xlab="Mb", t="l", panel.first=rect(xleft,ybottom, xright,ytop,col="skyblue", border=NA))
  abline(v=boundary,col="red")
  
  ytop <- c( 4,4,4,4,4)
  plot(binned$bin, binned$median.error/1.0e6, ylim=c(0,3), xlim=c(0,130),  ylab="median error/Mb", xlab="Position/Mb", t="l", panel.first=rect(xleft,ybottom, xright,ytop,col="skyblue", border=NA))
  abline(v=boundary,col="red")
  abline(h=overall.median, col="blue")
  print(overall.median)
  if ( ! is.null(pdffile)) 
    dev.off()

  return( binned)
}

