.packageName <- "happy.hbrem_2.2"

library(MASS)
library(g.data)
library(multicore)
	
                                        # C interface to read in .data and .alleles files, perform DP and create a happy object
                                        #
                                        # a happy object is a list with the following attributes:
                                        # happy$strains array of strain names
                                        # happy$markers array of marker names
                                        # happy$subjects array of subject names
                                        # happy$phenotype array of phenotype values
                                        # happy$map array of map locations (in cM) pf markers
                                        # happy$handle integer handle which maps the R happy object to the corresponding C QTL object


happy <- function( datafile, allelesfile, generations=200, phase="unknown",
                 file.format="happy", missing.code="NA", do.dp=TRUE,
                 min.dist=1.0e-5, mapfile=NULL, ancestryfile=NULL, haploid=FALSE ) {

  gen <- as.numeric(generations)+0
  if ( phase=="estimate" ) file.format  <- "ped"

  h <- .Call( "happy", datafile, allelesfile, gen, phase, file.format, missing.code, do.dp=as.integer(do.dp), min.dist=min.dist, haploid=as.integer(haploid), ancestryfile=ancestryfile, PACKAGE="happy.hbrem" )

  h$phase <- phase
  h$haploid <- haploid

  strain.names    <- make.names(h$strains)
  num.strains     <- length(h$strains)
  h$names.additive <- strain.names

  diplotype.names <- matrix(kronecker(strain.names, strain.names, paste, sep="."), nrow=num.strains)
  h$names.full.symmetric <- c( diag(diplotype.names), diplotype.names[upper.tri(diplotype.names, diag=FALSE)])
  h$names.full.asymmetric <- c( t(diplotype.names) ) # assumes row major order, ie, (row1, row2, etc), in C object
  
  # for backwards compatibility
  h$nam <- strain.names
  h$nam2 <- h$names.full.symmetric
  h$nam3 <- h$names.full.asymmetric
  h$bp = h$map
	
  if ( ! is.null(mapfile)) {
    map <- read.delim(mapfile)
    if ( !is.null(map$marker) && ! is.null(map$bp)) {
      hbp <- data.frame(markers=h$markers,bp=rep(NA,length(h$markers)))
      h$bp <- map$bp[match(hbp$markers,map$marker)]
    }
    else
      stop( "incorrect column names found in mapfile ", mapfile , "\n")
  }

  h$additive = list()
  nm = length(h$markers)-1
  h$additive$genome <- data.frame(
                         marker     = I(as.character(h$markers))[1:nm],
                         map        = as.numeric(h$map)[1:nm],
                         bp         = as.numeric(h$bp)[1:nm],
                         chromosome = I(as.character(h$chromosome))[1:nm])

  h$full = list()
  h$full$genome <- data.frame(
                         marker     = I(as.character(h$markers))[1:nm],
                         map        = as.numeric(h$map)[1:nm],
                         bp         = as.numeric(h$bp)[1:nm],
                         chromosome = I(as.character(h$chromosome))[1:nm])
  
  return(h)
}


happy.matrices <- function( h ) {

  if ( class(h) == "happy" ) {
    if (is.null(h$matrices) ) {
      matrices <- list()
      for( m in h$markers) {
        add <- hdesign( h, m, model='additive' )
        full <- hdesign( h, m, model='full' )
	if ( h$use.pedigrees || h$phase.known )
	 full.asymmetric <- hdesign( h, m, model="full.asymmetric")
	else
	 full.asymmetric <- NULL 
        id <- m
        matrices[[m]] <- list( id=id, additive=add, full=full, full.asymmetric=full.asymmetric )
      }
      h$matrices <- matrices
    }
  }
  return(h)
}


happy.save <- function( h, file ) {

  if ( class(h) == "happy" ) {
    if (is.null(h$matrices) )
      h <- h$matrices
    save(h,file=file,compress=TRUE)
  }
}


                                        # C interface to return the design matrix for a marker interval
                                        # h is a happy object returned by a previous call to happy
                                        # marker is the name of the left-hand marker in the interval, or
					#    the integer index of the marker (starting from 1)
                                        # model can be 'additive' or 'full'
                                        # if mergematrix is non-null then the columns of the design matrix are merged 

hdesign <- function( h, marker, model='additive', mergematrix=NULL ) {

  d<- NULL

  if ( class(h) == "happy" ) {
    if (  ! is.null(h$matrices) ) { # data are in R memory
    
      if ( !is.integer(marker) ) { # integer marker index
        m <- match( marker, h$markers )
        if ( ! is.na(m) )
          return( hdesign( h, m, model=model, mergematrix=mergematrix))
        else
          return(NULL)
      }
      else { # marker name
        obj <- h$matrices[marker]
        if ( model == 'additive' ) 
          return ( obj[[1]]$additive )
        else
          return ( obj[[1]]$full )
      }
    }
    else { # data are in C memory
      handle <- as.numeric(h$handle)+0
      if ( h$haploid  ) {
        d <- .Call( "haploid_happydesign", handle, marker,  PACKAGE="happy.hbrem")
        if ( ! is.null(d) ) 
          if ( model == 'additive' ) colnames(d) <- h$nam
          else colnames(d) <- h$nam2
      }
      else {
        d <- .Call( "happydesign", handle, marker, model, PACKAGE="happy.hbrem")
        if ( ! is.null(d) ) {
          if ( model == 'additive' ) colnames(d) <- h$nam
          else colnames(d) <- h$nam2
        }
      }
    }
  }
  else if ( class(h) == "happy.genome" ) { # delayed data package
    loaded.markers <- load.markers( h, c(marker), model=model )
    d <- loaded.markers[[1]]
  }

  # merge the matrix if required
  
  if ( ! is.null(d) && ! is.null( mergematrix ) ) {
    if ( model == 'additive' ) {
      d <- d %*% mergematrix$amat
      colnames(d) <- colnames(mergematrix$amat)
    }
    else if (model == "full" ) {
      d <- d %*% mergematrix$imat
      colnames(d) <- colnames(mergematrix$imat)
    }
    else if (model == "full.asymmetric") {
      d <- d %*% mergematrix$famat
      colnames(d) <- colnames(mergematrix$famat)
    }
  }
  return(d)
}

hnonrecomb<- function( h, marker=NULL, do.mean=TRUE ) {
   handle <- as.numeric(h$handle)+0
   if ( is.null(marker) ) {
     nm <- length(h$markers)-1
     r <- vector("numeric",length=nm)
     for(m in 1:nm) {
       x <- .Call( "happynonrecomb", handle, m, PACKAGE="happy.hbrem")
       r[m] <- 0.5*mean(x)
     }
   }
   else {
     r <- .Call( "happynonrecomb", handle, marker, PACKAGE="happy.hbrem")
     r <- 0.5*r
     if ( do.mean )
       r <- mean(r)
   }
   return(r)
}
  

hprob <- function( h, marker=NULL ) {

   handle <- as.numeric(h$handle)+0
   p <- .Call( "happyprobs", handle, marker, PACKAGE="happy.hbrem")
   colnames(p) <- h$nam2
   return(p)
}

hprob2 <- function( h, marker=NULL, symmetrize=FALSE ) {
   handle <- as.numeric(h$handle)+0
   if ( symmetrize==TRUE )
     symmetrize=1
   else
     symmetrize=0

   p <- .Call( "happyprobs2", handle, marker, symmetrize, PACKAGE="happy.hbrem")
   if ( symmetrize==1 )
     colnames(p) <- h$nam2
   else
     colnames(p) <- h$nam3
   return(p)
}

h.sum.prob2 <- function( h, marker=NULL ) { # for Caroline - the sum of squares of the probabilities
   handle <- as.numeric(h$handle)+0
   if ( ! is.null( marker ) ) {
     p <- .Call( "happyprobs2", handle, marker, PACKAGE="happy.hbrem")
     p2 <- apply( p*p, 1, sum )
     return(p2)
   }
   else {
     nm <- length(h$markers)-1
     mat <- matrix( nrow=length(h$subjects), ncol=nm )
     for(i in 1:nm ) {
       p <- .Call( "happyprobs2", handle, i, PACKAGE="happy.hbrem")
       mat[,i] <- apply( p*p, 1, sum )
     }
     rownames(mat) <- h$subjects
     colnames(mat) <- h$markers[1:nm]
     return(mat)
   }
}


hgenotype <- function( h, marker=NULL, collapse=FALSE, sep="" ) {

  if ( class(h) == "happy" ) {
    handle <- as.numeric(h$handle)+0
    g <- .Call( "happygenotype", handle, marker, PACKAGE="happy.hbrem")
  }
  else if ( class(h) == "happy.genome" ) { # delayed data package
    loaded.markers <- load.markers( h, c(marker), model="genotype" )
    g <- loaded.markers[[1]]
  }
  if ( collapse ) {
    y <- paste( g[,1], g[,2], sep=sep )
    g <- ifelse ( y == "NANA", NA, y )
  }
  else {
    colnames(g) <- c("allele1", "allele2")
  }
  return(g)
}
	

                                        # fit a QTL to the markers using the specified mode

                                        # happy is a happy object returned by a call to happy
                                        # markers is a array of marker names or marker indices
                                        # mode is one of 'additive', 'full' or 'partial'
                                        # verbose controls the amount of output
                                        # return value is a table giving the logP values for each marker tested

hfit <- function( h, markers=NULL, model='additive', mergematrix=NULL, covariatematrix=NULL, verbose=FALSE, phenotype=NULL, family='gaussian', permute=0 ) {

  if ( class(h) == "happy.genome" ) {
    if ( !is.null( h[[model]] ) )
      map <- h[[model]]$map
    if ( is.null(markers) ) {
      nm <- length(h[[model]]$markers)-1
      markers <- h[[model]]$markers[1:nm]
    }
    if ( is.null(phenotype)) {
      error( "phenotype must be set\n")
    }
  }
  else {
    
    map <- h$map

    if ( is.null(markers) ) 
      markers <- h$markers[1:length(h$markers)-1]

    if ( is.null(phenotype) )
      phenotype = h$phenotype
  }
  
  if ( model == 'partial' || model == 'full' ) {
    lp <- matrix( ncol=8, nrow=length(map)-1)
    colnames(lp) <-c('cM', 'marker', 'additive logP', 'full logP', 'partial logP', 'additive SS', 'full SS', 'partial SS')
    offset <- 3
    width <- 3
    idx <- 3
  }
  else {
    lp <- matrix( ncol=4, nrow=length(map)-1)
    colnames(lp) <-c('cM', 'marker', 'additive logP', 'additive SS')
    offset <- 3
    width <- 1
    idx <- 3
  }    

  permdata <- NULL
  if ( permute > 0 ) { # permutation test
    offset <- 6
    width <- 2
    
    
    if ( verbose ) cat ("performing permutation anaysis on ", permute, " permutations\n")
    hf0 <- hfit( h, markers=markers, model=model, mergematrix=mergematrix, covariatematrix=covariatematrix, verbose=FALSE, family=family, permute=0 )
    emp <- matrix( ncol=7, nrow=length(map)-1)
    colnames(emp) <-c('cM', 'marker', 'anova.logP', 'global.pval',  'point.pval', 'global.logP', 'point.logp')
    maxlogp <- vector(mode='numeric',length=permute)
    logpk <- vector(mode='numeric',length=length(map)-1)
    for(k in 1:permute) {
      shuf <- sample(phenotype)
      pres <- hfit( h, markers=markers, model=model, mergematrix=mergematrix, covariatematrix=covariatematrix, verbose=FALSE, phenotype=shuf, family=family, permute=0 )
      maxlogp[k] <- max( as.numeric(pres$table[,idx]) )
      logpk <- logpk + ifelse(pres$table[,idx] > hf0$table[,idx], 1, 0 )
      if ( verbose ) cat(k, maxlogp[k], "\n")
    }
    maxlogp <- sort( maxlogp, decreasing=TRUE )
    logpk <- logpk / permute
    p01 <- maxlogp[as.integer(permute/100)]
    p05 <- maxlogp[as.integer(permute/20)]
    if ( verbose) cat( 'p01 ', p01, ' p05 ', p05 , "\n")
    mi <- 1

    for ( m in hf0$table[,"marker"]) {
      emp[mi,4] <- NA
      if ( mi <= nrow(hf0$table)) {
        logp <- as.numeric(hf0$table[mi,idx])
        n <- 0
        if ( is.numeric(logp) ) {
          n <- sum(ifelse(maxlogp>logp,1,0),na.rm=TRUE)
          emp[mi,"global.pval"] <- n/permute
        }
      }
      mi <- mi+1
    }
    logemp <- ifelse( emp[,'global.pval'] >= 1/permute , -log10(emp[,'global.pval']), log10(permute))
    emp[,"cM"] <- hf0$table[,"cM"]
    emp[,"marker"] <- hf0$table[,"marker"]
    emp[,"anova.logP"] <-hf0$table[,idx]
    emp[,"global.logP"] <- logemp
    emp[,'point.pval'] <- logpk
    emp[,'point.logp'] <- ifelse( logpk >= 1/permute, -log10(logpk), log10(permute))
    permdata <- list( N=permute, p01=p01, p05=p05, permutation.dist=maxlogp, permutation.pval=emp )
    hf0$width <- width
    hf0$offset <- offset
    hf0$permdata <- permdata
    return(hf0);
  }
  else {
    i <- 1
    maxp <- 0
    maxm <- NA
    maxSS <- NA
    
    for( m in markers ) {
      if ( verbose ) cat( "\n\n****** ", i, "marker interval ",  m, "\n\n" )
      if ( model == 'partial' || model == 'full' ) {
        if ( ! is.null(full<- hdesign( h, m, model='full', mergematrix=mergematrix )) ) {
          additive <- hdesign( h, m, model='additive', mergematrix=mergematrix )
          if ( is.null(covariatematrix) ) {
            cfit <- glmfit( phenotype ~ 1, family=family)
            ffit <- glmfit( phenotype ~ full , family=family)
            afit <- glmfit( phenotype ~ additive , family=family)
          }
          else {
            cfit <- glmfit( phenotype ~ covariatematrix , family=family)
            full <- cbind(covariatematrix, full )
            ffit <- glmfit( phenotype ~ full , family=family)
            additive <- cbind(covariatematrix, additive )
            afit <- glmfit( phenotype ~ additive , family=family)
          }
          
          if ( family == "gaussian") {
            an <- anova( cfit, afit, ffit );
            if ( verbose ) print(an)
            logP <- -log(an[[6]])/log(10)
          
            an2 <- anova( cfit, ffit );
            logP2 <- -log(an2[[6]])/log(10)
          }
          else {
            an <- anova( cfit, afit, ffit, test="Chisq" );
            if ( verbose ) print(an)
            logP <- -log(an[[5]])/log(10)
          
            an2 <- anova( cfit, ffit, test="Chisq" );
            logP2 <- -log(an2[[5]])/log(10)
          }

          if ( ! is.na(logP[2]) && logP[2] > maxp ) {
            maxp <- logP[2]
            maxm <- m
            maxSS <- an2[[5]][2]
          }

          lp[i,1] <- (map[i]+map[i+1])/2
          lp[i,2] = m
          lp[i,3] <- logP[2]
          lp[i,4] <- logP2[2]
          lp[i,5] <- logP[3]
          lp[i,6] <- an[[4]][2]
          lp[i,7] <- an2[[4]][2]
          lp[i,8] <- an[[4]][3]
          i <- i+1
        }
      }
      else {
        if ( ! is.null(d<- hdesign( h, m, model='additive', mergematrix=mergematrix )) ) {
          if ( is.null(covariatematrix) ) {
            cfit <- glmfit( phenotype ~ 1 , family=family)
            afit <- glmfit( phenotype ~ d , family=family)
          }
          else {
            cfit <- glmfit( phenotype ~ covariatematrix )
            d <- cbind( covariatematrix, d )
            afit <- glmfit( phenotype ~ d )
          }

          if ( family=="gaussian") {
            an <- anova( cfit, afit )
            logP <- -log( an[[6]])/log(10)
          }
          else {
            an <- anova( cfit, afit, test="Chisq" )
#            print(an)
            logP <- -log( an[[5]])/log(10)
#            print(logP)
          }
            
          if ( verbose ) {
            print( an )
            strain.effects( h, afit )
          }



          if ( ! is.na(logP[2]) && logP[2] > maxp ) {
            maxp <- logP[2]
            maxm <- m
            maxSS <- an[[4]][2]
          }

          lp[i,1] <- (map[i]+map[i+1])/2
          lp[i,2] <- m
          lp[i,3] <- logP[2]
          lp[i,4] <- an[[4]][2]

          i <- i+1
        }
      }
    }
    return(list( table=lp, model=model, family=family, test='hfit', offset=offset, width=width, maxp=maxp, maxm=maxm, maxSS=maxSS, permdata=NULL ))
  }
}

                                        # Calculate the T-tests for comparing strain effects.
                                        # Called from hfit: only suited to the additive model at present

strain.effects <- function( h, fit, family='gaussian' ){

  if ( class(fit) == "lm" ) {
    c <- coef(fit)
    len <- length(c)
    cov <- vcov(fit)
    dia <- 1/sqrt(diag(cov))
    v <- dia %o% dia
    corr <- cov * v

    cat('\nCorrelation Matrix:\n')
    print (corr)
    cat('\n')
    nam <- names(c)
    nam2 <- c( 'Mean', h$strains ) # the true strain names

    
    se <- c()
    for(m in 1:len ) {
      strain1 <- nam[m]
      if ( ! is.na(c[strain1]) ) 
        se <- c( se, sqrt(cov[strain1,strain1]))
      else
        se <- c( se, NA )
    }
    
    df <- df.residual(fit)
    intercept <- "(Intercept)"

    cat('\nStrain Main Effects, with standard errors:\n' )
    cat(formatC( nam2, width=10, format='s' ), '\n')
    cat(formatC( c, width=10, digits=4, format='f'), '\n' )
    cat(formatC( se, width=10, digits=4, format='f'), '\n' )


    
    e <- matrix( nrow=len*(len+1)/2, ncol=4 )
    en <- matrix( nrow=len*(len+1)/2, ncol=2 )
    colnames(e) <- c(  'diff', 'se', 'T', 'P' )
    colnames(en) <- c( 'strain1', 'strain2' )
    k <- 1
    for(m in 1:len ) {
      strain1 <- nam[m]
      if ( strain1 != intercept && ! is.na(c[strain1]) ) {
        for(n in 1:m) {
          strain2 <- nam[n]
          if ( strain2 != intercept && ! is.na(c[strain2]) ) {
            if ( n < m ) {
              d <- c[m]-c[n]
              v <- cov[strain1,strain1] + cov[strain2,strain2] - 2*cov[strain1,strain2]
              en[k,1] <- nam2[m]
              en[k,2] <- nam2[n]
              e[k,1] <- d
              if ( v > 1.0e-6 ) {
                se <- sqrt(v)
                t <- d/se
                p <- pt(t,df=df,lower.tail=FALSE)
                e[k,2] <- se
                e[k,3] <- t
                e[k,4] <- p
              }
              k <- k+1
            }
            else if ( n == m ) {
              d <- c[m] -c[intercept]
              v <- cov[strain1,strain1] + cov[strain2,intercept] - 2*cov[strain1,intercept]
              en[k,1] <- nam2[m]
              en[k,2] <- '(Mean)'
              e[k,1] <- d
              if ( v > 1.0e-6 ) {
                se <- sqrt(v)
                t <- d/se
                p <- pt(t,df=df,lower.tail=FALSE)
                e[k,2] <- se
                e[k,3] <- t
                e[k,4] <- p
              }
              k <- k+1
            }
          }
        }
      }
    }

    k <- k-1
    cat('\nTests of Strain Differences (note that differences may be hard to interpret when strains are indistinguishable)\n\n')
    cat(
        formatC('strain1',width=12,format='s'),
        formatC('strain2',width=12,format='s'),
        formatC('diff',width=10,format='s'),
        formatC('se', width=10, format='s'),
        formatC('T', width=7, format='s'),
        formatC('P',width=7,format='s'),
        '\n')
    
    for( j in 1:k) {
      cat(formatC( en[j,'strain1'], width=12,format='s'),
          formatC( en[j,'strain2'], width=12,format='s'),
          formatC(e[j,'diff'], width=10, digits=4, format='f'),
          formatC(e[j,'se'], width=10, digits=4, format='f'),
          formatC(e[j,'T'], width=7, digits=3, format='f'),
          formatC(e[j,'P'],digits=4,width=7,format='e'),
          '\n')
    }
    e[1:k,]
  }
  NULL
}


                                        # merge fit a QTL to the markers using the specified mode

                                        # happy is a happy object returned by a call to happy
                                        # markers is a array of marker names or marker indices
                                        # model is one of 'additive', 'full' or 'partial'
                                        # covariatematrix is an optional design matrix to include on all the models.
                                        # (This can be used to include additional markers or covariates)
                                        # verbose controls the amount of output
                                        # merge contains the merge object used for merging strains
                                        # return value is a table giving the logP values for each marker tested
                                        # the test is for merged strains versus unmerged

mfit <- function( happy, markers=NULL, model='additive', mergematrix=NULL, covariatematrix=NULL, verbose=TRUE, family='gaussian', variants=NULL  ) {

  map <- happy$map
  if ( is.null( markers ) )
    markers <- happy$markers
  
  if ( model == 'partial' || model == 'full' ) {
    lp <- matrix( ncol=8, nrow=length(markers))
    if ( is.null(variants) )
      rownames(lp) <- markers
    else
      rownames(lp) <- variants
      
    colnames(lp) <-c('cM', 'marker', 'full-merged', 'full-unmerged', 'partial', 'full-merged-SS', 'full-unmerged-SS', 'partial-SS' )
  }
  else {
    lp <- matrix( ncol=8, nrow=length(markers))
    if ( is.null(variants) )
      rownames(lp) <- markers
    else
      rownames(lp) <- variants
    colnames(lp) <-c('cM', 'marker', 'additive-merged', 'additive-unmerged', 'partial', 'additive-merged-SS', 'additive-unmerged-SS', 'partial-SS' )
  }    

  i <- 1
  maxp <- 0
  maxm <- NA
  for( m in markers ) {
    if ( verbose ) print( paste( " ", i, "marker interval ",  m ) )
    if ( model == 'partial' || model == 'full' ) {
      if ( ! is.null(full<- hdesign( happy, m, model='full', mergematrix=mergematrix )) ) {
        fullu <- hdesign( happy, m, model='full', mergematrix=NULL)
        if ( is.null( covariatematrix ) ) {
          cfit <- glmfit( happy$phenotype ~ 1 , family=family)
          ffit <- glmfit( happy$phenotype ~ full , family=family)
          fufit <- glmfit( happy$phenotype ~ fullu , family=family)
        }
        else {
          cfit <- glmfit( happy$phenotype ~ covariatematrix , family=family)
          full <- cbind( covariatematrix, full )
          ffit <- glmfit( happy$phenotype ~ full , family=family)
          fullu <- cbind( covariatematrix, fullu )
          fufit <- glmfit( happy$phenotype ~ fullu , family=family )
        }          
        an <- anova( cfit, ffit, fufit );
        if ( verbose ) print(an)
        logP <- -log(an[[6]])/log(10)
        if ( ! is.na(logP[2]) && logP[2] > maxp ) {
          maxp <- logP[2]
          maxm <- m
        }

        an2 <- anova( cfit, fufit );
        logP2 <- -log(an2[[6]])/log(10) # p-value for full unmerged

        
        lp[i,1] <- (map[i]+map[i+1])/2
        lp[i,2] <- m
        lp[i,3] <- logP[2] # full merged 
        lp[i,4] <- logP2[2] # full unmerged
        lp[i,5] <- logP[3] # partial test of merged vs unmerged
        if ( is.na(lp[i,5]) ) lp[i,5] <- 0
        lp[i,6] <- an[[4]][2]
        lp[i,7] <- an2[[4]][2]
        lp[i,8] <- an[[4]][3]
                                        #        lp[is.na(lp)] = 0.0
        i <- i+1
      }
    }
    else {
      if ( ! is.null(d<- hdesign( happy, m, model='additive', mergematrix=mergematrix )) ) {
        du <- hdesign( happy, m, model='additive', mergematrix=NULL )

        if ( is.null( covariatematrix ) ) {
          cfit <- glmfit( happy$phenotype ~ 1 , family=family)
          afit <- glmfit( happy$phenotype ~ d , family=family)
          aufit <- glmfit( happy$phenotype ~ du , family=family)
        }
        else {
          cfit <- glmfit( happy$phenotype ~ covariatematrix , family=family)
          d <- cbind( covariatematrix, d )
          afit <- glmfit( happy$phenotype ~ d , family=family)
          du <- cbind( covariatematrix, du )
          aufit <- glmfit( happy$phenotype ~ du , family=family)
        }
        an <- anova( cfit, afit, aufit )
        if ( verbose ) print( an )

        logP <- -log(an[[6]])/log(10)
        if ( ! is.na(logP[2]) && logP[2] > maxp ) {
          maxp <- logP[2]
          maxm <- m
        }

        an2 <- anova( cfit, aufit );
        logP2 <- -log(an2[[6]])/log(10)

        lp[i,1] <- (map[i]+map[i+1])/2
        lp[i,2] <- m
        lp[i,3] <- logP[2]
        lp[i,4] <- logP2[2]
        lp[i,5] <- logP[3]
        if ( is.na(lp[i,5]) ) lp[i,5] <- 0
        lp[i,6] <- an[[4]][2]
        lp[i,7] <- an2[[4]][2]
        lp[i,8] <- an[[4]][3]


        i <- i+1
      }
    }
  }
  list( table=lp, model=model, test='mfit', offset=3, width=3, maxp=maxp, maxm=maxm )
}

                                        # fit two markers simultaneously, assuming additive model

twofit <- function ( happy, marker1, marker2, merge1=NULL, merge2=NULL, model = 'additive', verbose=TRUE, family='gaussian' ) {
  
  if ( ! is.null(a1<- hdesign( happy, marker1, mergematrix=merge1, model=model )) && ! is.null(a2<- hdesign( happy, marker2, mergematrix=merge2, model=model )) ) {
    print ( paste( 'joint test of markers ', marker1, ', ', marker2 ) )
    both <- cbind( a1, a2 )
    f1 <- glmfit( happy$phenotype ~ a1 , family=family)
    f2 <- glmfit( happy$phenotype ~ a2 , family=family)
    b <-  glmfit( happy$phenotype ~ both , family=family)
    cfit <- glmfit( happy$phenotype ~ 1 , family=family)
    s <- summary(glmfit(f1))
    f <- s$fstatistic
    p<-pf( f[1], f[2], f[3], lower.tail=FALSE, log.p=TRUE )
    if ( verbose ) {
      print(anova( cfit, f1, b ))
      print(anova( cfit, f2, b ))
      print(names(s))
      print(-p/log(10))
    }
  }
}



                                        # fit a model to a set of markers conditional upon another marker being in the model

condfit <- function( happy, markers, condmarker, merge=NULL, condmerge=NULL, model='additive', condmodel='additive', epistasis=FALSE, verbose=TRUE,family='gaussian' ) {
  if ( ! is.null( cond <- hdesign( happy, condmarker, model=model, mergematrix=condmerge ))) {
    cfit <- glmfit( happy$phenotype ~ 1  , family=family)
    condfit <- glmfit( happy$phenotype ~ cond , family=family)
    if ( verbose ) print ( paste( 'conditional additive test on marker ', condmarker ) )
    for( m in markers ) {
      additive <- hdesign( happy, m, model=model, mergematrix=merge )
      condadditive <- cbind( additive, cond ); 
      cafit <- glmfit( happy$phenotype ~ condadditive , family=family)
      if ( verbose ) print(anova(cfit,condfit,cafit))
    }
  }
}

mergedpositionmatrix <- function( h, position, prepmerge, model='additive', verbose=FALSE, design=TRUE ) {

                                        # creates the merged design matrix for the specified position
                                        # optionally returns the mergematrix instead
  
  ind <- match( position, prepmerge$testmarkerdata$POSITION ) # the index of the position
  retval <- NULL
  if ( ! is.null(ind) ) {
    alleles <- prepmerge$testmarkerdata[ind,] # the allele distribution amounst the strains
    mlist <- mergelist( h$strains, alleles ) # list of lists representation of the allele distribution 
    mergematrix <- mergematrices( h$strains, mergelist=mlist) # the matrices representing the merge
    if ( design ) { # return the design matrix
      marker <- prepmerge$interval[ind]
      retval <- hdesign( h, marker, merge=mergematrix, model=model )
    }
    else # return the mergematrix
      retval <- mergematrix
  }
  else {
    print(paste('test marker', position, 'not found' ))
  }
  retval
}



                                        # merge strains together.
                                        # creates matries used to multiply with hapy design matrices
                                        # Subsequent analyses of the happy data can use the merged strains.
                                        # happy is a happy object
                                        # mergedata is a list of lists of strains eg c( c( "Balbc", "AKR"), c("C57", "DBA") )
                                        # NOTE: all strain names MUST be in happy$strains
                                        # NOTE: unlisted strains are dropped

mergematrices <- function ( strains, mergelist=NULL, verbose=FALSE ) {

  if ( is.null( mergelist ) ) {
    return(NULL)
  }
  else {
    ls <- length(strains)
    lm <- length(mergelist)
    merge <- strains
    group <- c()

                                        # the matrix for the additive model
    
    amat <- matrix( nrow=ls, ncol=lm, data=0 )
    mapping <- 1:ls
    
    i <- 1
    for ( m in mergelist ) {
      ind <- match( m, strains, nomatch=NA);
      merge[ind] = paste( sep='','group',i)
      amat[ind,i] = 1
      p <- paste(as.list(strains[ind]),collapse=",")
      group <- c( group, p)
      mapping[ind] = i
      i <- i+1
    }

                                        # the matrix for the full model
    
    imat <- matrix( nrow=ls*(ls+1)/2, ncol=lm*(lm+1)/2, data=0 )

    n <- 0
    imap <- matrix(nrow=lm,ncol=lm)
    for(i in 1:lm) {
      for(j in 1:i ) {
        n <- n + 1
        imap[i,j] <- n
        imap[j,i] <- n
      }
    }

    iname <- array( dim=lm*(lm+1)/2 )
    n <- 0

#    for(i in 1:ls) {
#      n <- n+1
#      k<- imap[mapping[i],mapping[i]]
#      imat[n,k] <- 1
#      iname[k] <- paste(sep="", "group", mapping[i], ",", mapping[i] )
#    }
#    for(i in 2:ls) {
#     for(j in 1:(i-1) ) {
#       n <- n+1
#       k<- imap[mapping[i],mapping[j]]
#        imat[n,k] <- 1
#       iname[k] <- paste(sep="", "group", mapping[i], ",", mapping[j] )
#
#      }
#    }


    for(i in 1:ls) {
      for(j in 1:i) {    
        n <- n+1
        k<- imap[mapping[i],mapping[j]]
        imat[n,k] <- 1
        iname[k] <- paste(sep="", "group", mapping[i], ",", mapping[i] )
      }
    }
    

	# the matrix for the full asymmetric model

    famat <- matrix( nrow=ls*ls, ncol=lm*lm, data=0 )

    n <- 0
    famap <- matrix(nrow=lm,ncol=lm)
    for(i in 1:lm) {
      for(j in 1:lm ) {
        n <- n + 1
        famap[i,j] <- n
      }
    }

    faname <- array( dim=lm*lm )
    n <- 0
    for(i in 1:ls) {
      for(j in 1:ls ) {
        n <- n+1
        k<- famap[mapping[i],mapping[j]]
#        print( c(i, j, n, k, mapping[i], mapping[j]))
        famat[n,k] <- 1
        faname[k] <- paste(sep="", "group", mapping[i], ",", mapping[j] )

      }
    }

######
	
    names(merge) <- strains
    colnames(amat) <- group
    colnames(imat) <- iname
    colnames(famat) <- faname
	
    if ( verbose ) print(paste("strains merged into ", length(group), "groups", paste( group, collapse="| " )))
    list(merge=merge, group=group, amat=amat, imat=imat, famat=famat)
  }
}

fastmergefit <- function( datafile, allelesfile, markerposfile, testmarkerfile, generations=200, model='additive', verbose=FALSE ) {

  h <- happy( datafile, allelesfile, generations=generations ) 
  if ( ! is.null( h ) ) {
    prep <- mergeprepare( h, markerposfile, testmarkerfile )
    fit <- mergefit(  h, prep, model=model, verbose=verbose )
    return(fit)
  }
  else {
    return(NULL)
  }
}

mergeprepare <- function( h, markerposfile, testmarkerfile, verbose=FALSE ) {
  
                                        # read in the marker positions from markerposfile
                                        # format comprises two columns header 'marker' and 'POSITION' 
                                        # the file must be sorted by position. Position refers to bp coordinate

  mergedata <- list()
  
  markerpos <- read.table(markerposfile,header=TRUE)
  if ( is.null(markerpos$POSITION) ) {
    print (paste(' ERROR - required column POSITION is missing from ', markerposfile ))
    return (NULL)
  }  
  if ( is.null(markerpos$marker) ) {
    print (paste(' ERROR - required column marker is missing from ', markerposfile ))
    return (NULL)
  }  

                                        # determine those markers that are also in the happy object

  markerpos$marker <- as.character(markerpos$marker)

  if ( verbose ) print(markerpos)

  if ( class(h) == "happy" ) 
    mmatch <- markerpos[match( h$markers, markerpos$marker, nomatch=0 ),]
  else 
    mmatch <- markerpos[match( h$additive$markers, markerpos$marker, nomatch=0 ),]
  nmarkers <- nrow(mmatch)
  markers <- mmatch[order(mmatch$POSITION),]
  if ( verbose ) print(markers)
  print(paste( nrow(markers), ' skeleton markers' ))

                                        # read in a file of markers to test 
                                        # must contain columns titled 'marker' and 'POSITION' and an optional column "Variant"
                                        # together with a column for every strain in happy$strains
  
  testmarkerdata<-read.table(testmarkerfile, header=TRUE)
  
  if ( is.null(testmarkerdata$Variant) )
    testmarkerdata$Variant <- paste( "var", as.character(testmarkerdata$POSITION), sep=".")


                                        # check that the required columns are all present

  required <- c( 'marker', 'POSITION' , make.names(h$strains) )
  s <-  required %in% names(testmarkerdata)
  n <- 1
  m <- 0  
  if ( verbose ) print(s)
  if ( verbose ) print(required)
  
  for ( x in s ) {
    if ( ! x ) {
      print( paste( 'required column', required[n], 'missing from', testmarkerfile ) );
      m <- m+1
    }
    n <- n+1
  }

  if ( m > 0 ) {
    print( 'mergeprepare halted')
    return(NULL)
  }


                                        # identify the happy marker interval for each marker to be tested

  ntest <- nrow(testmarkerdata)
  interval <- matrix(nrow=ntest,ncol=1)
  last <- 1
  good <- 0
  
  for ( t in 1:ntest ) {
    i <- last
    pos <- testmarkerdata$POSITION[t]

    ok <- FALSE
    while( ! ok && i < nmarkers ) { 
      if ( pos >= mmatch$POSITION[i] && pos < mmatch$POSITION[i+1] ) {
        interval[t] <- as.character(mmatch$marker[i])
        if ( verbose ) cat(paste(sep=' ', good, as.character(testmarkerdata$marker[t]), interval[t], mmatch$POSITION[i],pos,mmatch$POSITION[i+1], '\n'))
        last <- i
        ok <- TRUE
        good <- good + 1
      }
      i <- i +1
    }

                                        # out of range markers
    
    if ( is.na(interval[t])) {
      if ( pos <= mmatch$POSITION[1] ){
        interval[t] = as.character(mmatch$marker[1])
        cat(paste('marker ', testmarkerdata$marker[t], 'out of range, assigned to interval ', interval[1], '\n'  ))
      } 
      else if ( pos >= mmatch$POSITION[nmarkers] ) {
        interval[t] = as.character(mmatch$marker[nmarkers-1]) 
        cat(paste('marker ', testmarkerdata$marker[t], 'out of range, assigned to interval ', interval[nmarkers-1], '\n' ))
      }
    }
  }

  print(paste(good, ' test markers placed on skeleton map'))
  
  mergedata$markerpos <- markerpos
  mergedata$interval <- interval
  mergedata$markers <- markers
  mergedata$testmergedata <- testmarkerdata
  mergedata
}

                                        # compute the strain distribution pattern from a list of alleles for the strains.
                                        # results is a string of 0'1 and 1's, in the order dictated by strains. The first character is always 0

sdp <- function( strains, alleles ) {

  s <- list()
  n <- 0
  for (a in alleles) {
    for (x in a) {
      s[x] = n
    }
    n <- n+1
  }


  flip <- FALSE
  if ( s[strains[1]] == 1 )
    flip <- TRUE

  if ( n == 2 && flip ) {
    for ( y in strains) {
      if ( s[y] == 1 )
        s[y] = 0
      else
        s[y] = 1
    }
  }

  sdpvalue <- ""
  for ( y in strains) {
    sdpvalue <- paste( sdpvalue, s[y], sep="" )
  }
  return(sdpvalue)
}


condmergefit <- function( h, mergedata, model='additive', covariatematrix=NULL, verbose=FALSE ) {

  interval <- mergedata$interval
  testmergedata <- mergedata$testmergedata
  ntest <- nrow(testmergedata)

  matrices <- list()
  strains <- h$strains
  logPmax <- matrix(nrow=ntest,ncol=6)
  colnames(logPmax) <- c( "position", "interval", "sdp", "logPself", "logPmax", "logPmaxPosition" )
  logP <- list()
  logPm <- list()
  logPmatrix <- matrix(nrow=ntest,ncol=ntest)
  for ( m in 1:ntest ) {
    if ( !is.null(interval[m]) ) { # this test eliminates markers that could not be placed on the map
      im <- interval[m]
      alleles <- testmergedata[m,] # the allele distribution amounst the strains
      mlist <- mergelist( strains, alleles ) # list of lists representation of the allele distribution
      sdpvalue <- sdp( strains, mlist ) # Strain Distribution pattern
      if ( is.null( matrices[[sdpvalue]] ) )
        matrices[[sdpvalue]] <- mergematrices( strains, mergelist=mlist) # the matrices representing the merge
      d <- hdesign( h, im, model=model, mergematrix=matrices[[sdpvalue]] )
      
      ckey <- paste(im,sdpvalue) # cache key
      if ( is.null(logP[[ckey]]) ) {
        if ( ! is.null(covariatematrix) )
          d <- cbind( d, covariatematrix )

        f <- mfit( h, im, model=model, mergematrix=matrices[[sdpvalue]], covariatematrix=covariatematrix, verbose=verbose )
        mf <- mergefit( h, mergedata, model=model, covariatematrix=d, verbose=verbose )
        mx <- max(as.numeric(mf$table[,3]),na.rm=TRUE)
        mp <- which.max(as.numeric(mf$table[,3]))
        logPm[[ckey]] <- ifelse( is.na(mf$table[,3]), 0, mf$table[,3] )
        logP[[ckey]] <- c( im, sdpvalue, f$maxp, mx, mp)
      }
      logPmatrix[m,] <- logPm[[ckey]]
      logPmax[m,] <- c( testmergedata$POSITION[m], logP[[ckey]])
      cat(paste(logPmax[m,]),"\n")
    }        
  }
  rownames(logPmatrix) <- as.character(testmergedata$POSITION)
  colnames(logPmatrix) <- as.character(testmergedata$POSITION)
  return( list(logPmax=logPmax, logPmatrix=logPmatrix ) )
}



mergefit <- function( h, mergedata,  model='additive', covariatematrix=NULL, verbose=FALSE ) {

  fit <- NULL
  n <- 1
  interval <- mergedata$interval
  testmergedata <- mergedata$testmergedata
  ntest <- nrow(testmergedata)
  
  cache <- list()
  
  f <- NULL
  saved <- 0
  calculated <- 0

  strains <- h$strains

  for ( m in 1:ntest ) {
    if ( !is.null(interval[m]) ) { # this test eliminates markers that could not be placed on the map
      im <- interval[m]
      alleles <- testmergedata[m,] # the allele distribution amounst the strains
      mlist <- mergelist( strains, alleles ) # list of lists representation of the allele distribution
      sdpvalue <- sdp( strains, mlist ) # Strain Distribution pattern
      ckey <- paste(im,sdpvalue) # cache key
                                        #      cat( "ckey ", ckey, '\n'  )
      if ( is.null( cache[[ckey]] )) { # test if this fit has been cached 
        mergematrix <- mergematrices( strains, mergelist=mlist) # the matrices representing the merge
        f <- mfit( h, im, model=model, mergematrix=mergematrix, covariatematrix=covariatematrix, verbose=verbose, variants=testmergedata$Variant[m] )
        calculated <- calculated +1
        cache[[ckey]] <- f
      }
      else {
        saved <- saved +1
        f <- cache[[ckey]]
      }
      f$table[1,1] <- alleles$POSITION
      
      if ( n == 1 ) {
        fit <- f
      }	
      else {
        fit$table <- rbind( fit$table, f$table[1,] )
      } 
      n <- n+1
    }
  }
  
  rownames(fit$table) <- as.character(testmergedata$Variant)
  print(paste(saved, 'fits saved', calculated, 'fits calculated'))
  fit$testmarkerdata <- testmergedata
  fit$interval <- interval
  return(fit)
}


                                        # create a mergelist structure from a list of alleles associated with each strain

mergelist <- function( strains, alleles ) {
  l <- c()
  strains <- make.names(strains)
  for ( s in strains ) {
    l <- c( l, as.character(alleles[[s]] ))	
  }
  l <- unique(sort(as.character(l)))
  mlist <- list()
  for ( a in l ) {
    v <- c()
    for ( s in strains ) {
      if ( alleles[[s]] == a ) {
        v <- c( v, s )
      }
    }
    mlist[[a]] <- v
  }
  mlist
}


mergeplot <- function( fit, mergedata, mode='logP', xlab='bp', ylab=NULL, main=NULL,  t='p', pch=20, ... ) {

  def.par <- par(no.readonly = TRUE)# save default, for resetting...

                                        #  layout( matrix( c(2,1), nrow=2, ncol=1  ), widths = c( 1), heights=c( 5, 1 ) )
                                        #  layout.show(2)

  labels <- list( text=as.character(mergedata$markers$marker), POSITION=mergedata$markers$POSITION )
  
  happyplot( fit, mode=mode, labels=labels, xlab=xlab, ylab=ylab, main=main,  t=t, pch=pch, ... )


  par(def.par)#- reset to default

}

epistasis <- function( h, markers1, markers2=NULL, merge1=NULL, merge2=NULL, model='additive', verbose=FALSE, family='gaussian' ) {

  if ( is.null( markers2 ) ) { # all pairwise interactions 
    nmarkers = length(markers1)
    ninteractions <- nmarkers*(nmarkers-1)/2
    results <- matrix( nrow=ninteractions, ncol=7 )
    colnames(results) <- c( 'marker1', 'marker2', 'main1', 'main2', 'main1+main2', 'main1*main2', 'main1.main2' )
    
    r <- 1
    d <- list()
    length(d) <- nmarkers
    main <-list()
    length(main) <- nmarkers
    logten <- log(10)
    
                                        # precalculate the main effects
    
    for ( m in 1:nmarkers ) {
      d[[m]] <- hdesign( h, markers1[m], mergematrix=merge1, model=model )
      f <- glmfit( h$phenotype ~ d[[m]] , family=family)
      a <- anova(f)
      main[[m]] <- -log(a[1,5])/logten # its log-P value
    }
    
    for ( m1 in 2:nmarkers ) {
      ma1 <- markers1[m1]
      print(ma1)
      for( m2 in 1:(m1-1) ) {
        results[r,] <- epistasispair( h, ma1, markers1[m2], merge1=merge1, merge2=merge2, mode=mode, verbose=verbose, d1=d[[m1]], d2=d[[m2]], main1=main[[m1]], main2=main[[m2]] )
        r <- r+1
      }
    }
  }
  else {
    nmarkers1 = length(markers1)
    nmarkers2 = length(markers2)
    ninteractions <- nmarkers1*nmarkers2
    results <- matrix( nrow=ninteractions, ncol=7 )
    colnames(results) <- c( 'marker1', 'marker2', 'main1', 'main2', 'main1+main2', 'main1*main2', 'main1.main2' )
    
    r <- 1

    for ( m1 in 1:nmarkers1 ) {
      for( m2 in 1:nmarkers2 ) {
        results[r,] <- epistasispair( h, markers1[m1], markers2[m2], merge1=merge1, merge2=merge2, model=model, verbose=verbose )
        r <- r+1
      }
    }
  }

  results
}

epistasispair<- function( h, marker1, marker2, merge1=NULL, merge2=NULL, model='additive', verbose=FALSE, d1=NULL, d2=NULL, main1=0, main2=0, family='gaussian' ) {
  

  logten <- log(10)
  
                                        # fit the first marker
  if ( is.null(d1)) {
    d1 <- hdesign( h, marker1, mergematrix=merge1, model=model )

    m1 <- glmfit( h$phenotype ~ d1 , family=family)
    a1 <- anova(m1)
    main1 <- -log(a1[1,5])/logten # its log-P value
  }
  
                                        # fit the second marker
  if ( is.null(d2)) {

    d2 <- hdesign( h, marker2, mergematrix=merge2, model=model )
    
    m2 <- glmfit( h$phenotype ~ d2 , family=family)
    a2 <- anova(m2)
    main2 <- -log(a2[1,5])/logten
  }
  
  additive <- 0
  epistatic <- 0
  full <- 0
  
  if ( ! is.null(d1) && ! is.null(d2) ) {

    
    const <- glmfit( h$phenotype ~ 1 , family=family)

    
                                        # fit both markers additively
    
    d12 <- cbind( d1, d2 )
    m12 <- glmfit( h$phenotype ~ d12 , family=family)
    a12 <- anova(m12)
    additive <- -log(a12[1,5])/logten
    
                                        # fit both markers epistatically
    D12 <- matrixSquared( d1, d2 )
    M12 <- glmfit( h$phenotype ~ D12 , family=family)
    A12 <- anova(const,m12,M12)
    
    if ( verbose ) print(A12)
    additive <- -log(A12[2,6])/logten # the log P for the additive main effects
    epistatic <- -log(A12[3,6])/logten # the log P for the interaction after removing main effects
    an12 <- anova( M12 );
    full <- -log(an12[1,5])/logten # log p-value for full interaction
  }

                                        #  print(c( marker1, marker2, 	main1, main2, additive, full, epistatic ))
  c( marker1, marker2, 	main1, main2, additive, full, epistatic )
  
}


matrixSquared <- function( matrix1, matrix2 ) {
  dim1 <- dim(matrix1)
  dim2 <- dim(matrix2)
  I <- NULL
  if ( dim1[1] == dim2[1] && dim1[2] == dim2[2] ) {
    I <- matrix( nrow=dim1[1], ncol=dim1[2]*dim2[2])
    c12 <- 1
    for( c1 in 1:dim1[2] )
      for( c2 in 1: dim2[2] ) {
        I[,c12] <- matrix1[,c1]*matrix2[,c2]
        c12 <- c12+1
      }
  }
  I 
}


                                        # support for a seqential fit of marker intervals 
hfit.sequential<- function ( h, threshold=2,  markers=NULL, model='additive', mergematrix=NULL, covariatematrix=NULL, verbose=FALSE, family='gaussian') {


  if ( is.null(covariatematrix) ) 
    covariatematrix <- matrix(1, nrow=length(h$phenotypes),ncol=1)
  
  nullfit <- glmfit( as.formula('h$phenotypes ~ covariatematrix'), family=family )
  lastfit <- nullfit
  cat('\nNull Model\n')
  print(anova(nullfit))


  fit <- hfit( h, markers=markers, model=model, mergematrix=mergematrix, covariatematrix=covariatematrix, verbose=verbose, family=family );
  logP <- fit$maxp;
  m <- fit$maxm
  kmax <- length(h$markers)-1
  k <- 0
  intervals <- c()
  interval <- list()
  oldformula <- 'h$phenotypes ~ covariatematrix'
  submatrix <- covariatematrix
  
  while ( logP > threshold && k < kmax ) {
    k <- k+1

    intervals <- c(intervals, m)
    d <- hdesign( h, m, model=model, mergematrix=mergematrix )
    interval[[m]] <- d
    newformula <- paste(oldformula, ' + interval[[\'', m, '\']]', sep="")
    nextfit <- glmfit( as.formula(newformula), family=family )

    cat(paste( "\n", k, " *** intervals", paste(intervals, collapse=","), "logP", logP), "\n\n" )

    cat('\nPartial Comparison:\n')
    an1 <- anova(lastfit,nextfit)
    print(an1)
    cat('\nFull Comparison:\n')
    an2 <- anova(nullfit,nextfit)
    print(an2)
    print(summary.lm(nextfit))

    lastfit <- nextfit
    oldformula <- newformula
    submatrix <- cbind( submatrix, d )
    
    fit <- hfit( h, markers=markers, model=model, mergematrix=mergematrix, covariatematrix=submatrix, verbose=verbose, family=family );
    logP <- fit$maxp;
    m <- fit$maxm

  }
  NULL
}


                                        # support for multiple phenotypes
                                        # phen should be a data table of numeric phenotypes

pfit <- function( h, phen, markers=NULL, model='additive', mergematrix=NULL, covariatematrix=NULL, verbose=FALSE, family='gaussian') {
  if ( is.data.frame(phen) ) {
    pnames <- names(phen)
    if ( nrow(phen) == length(h$subjects) ) {
      results <- list();
      for ( p in pnames) {
        cat(p, "\n")
        results[[p]] <- hfit(h,phenotype=as.numeric(as.character(phen[[p]])),family=family,covariatematrix=covariatematrix)
      }
      return(results)
    }
    else {
      print(paste('ERROR - number of rows in phen = ', nrow(phen), ', different from number of subjects ', nrow(h$subjects)))
    }
  }
  else {
    print('ERROR phen is not a data.frame')
  }
  NULL
}

glmfit <- function( formula=NA, family='gaussian' ) {
  formula <- as.formula(formula)
  if ( family == 'gaussian' )
    return(lm( formula ))
  else
    return(glm( formula, family=family))
}


                                        # convert a vector of numbers into gaussian deviates

normalise <- function( values=NULL ) {

  n<- length(values)+1;
  r <- rank(values)/n;
  qnorm(r)
}



###################################



gaussian.iterate <- function( d, params ) {

	y <- d$y
	probs <- d$probs # n x p
	beta <- params$beta
	sigma2 <- params$sigma2
	n <- length(y)
        yy <- t(array(y, dim=c(n,length(beta))))
        pp <- array(beta,dim=c(length(beta),n))
        xx <- yy-pp  # p x n
	xxx <- xx*xx # p x n
        e <- t(exp(-xxx/(2*sigma2)))*probs  # n x p

	w <- e / apply( e,1,sum) # divide by col sums  n x p 
	ws <- apply(w,2,sum)
	beta.new <- drop(  y %*% w ) / ws

	rs <- apply( e, 1, sum )
	sigma2.new <- mean (apply( t(w) * xxx, 2, sum ) )

        LogL <- -sum(log(rs)) + n* log(2*pi*sigma2) *0.5

	dbeta <-  ( beta - beta.new )*ws/sigma2
	dsigma2 <- 0.5*n/sigma2 * (1 - sigma2.new/sigma2 )

	return (list( LogL=LogL, beta=beta.new, sigma2=sigma2.new, dbeta=dbeta, dsigma2=dsigma2 ))
}

gaussian.init <- function( d ) {

	dm<- dim(d$probs)
	df <- dm[2]
	n <- length(d$y)
	sigma2 <- var(d$y)*(n-1)/n
	mu <- mean(d$y)
	LogL <- gaussian.null( n, sigma2 )
	return( list( sigma2=sigma2, beta=array(mu, dim=df), LogL=LogL))
}
	
gaussian.null <- function( n, sigma2 ) {

	return( 0.5*n*(1+log(2*pi*sigma2)))
}

gaussian.loop <- function ( d, maxit=100, eps=1.0e-3, df=NULL ) {
	
	i <- 0
	e <- 2*eps 
	Flast <- 0.0
	params.null <- gaussian.init( d )
	params <- params.null

	while( i < maxit && e > eps ) {
		params.new <- gaussian.iterate( d, params ) 	
		i <- i+1
		params <- params.new
		e <- abs(params.new$LogL-Flast);
#		print(c(i,e,params.new$LogL,params.new$sigma))
		Flast <- params.new$LogL
		
	}

	params$it <- i
	params$eps <- e
	params$N <- length(d$y)
	params$Null <- gaussian.null(params$N,params.null$sigma2)
	params$chi <- 2*( params$Null - params.new$LogL )
	if ( is.null(df) ) 
		params$df <- length(params$beta)-1
	else
		params$df <- df

	print( c(params$chi, params$df , params$Null, params$LogL))
	params$Pval <- pchisq( params$chi, params$df , lower.tail=FALSE)
	params$LogPval <- -log10(params$Pval)
	
	if ( params$chi < 0 ) {
	cat (c("error ", params$chi, params$Null , params$LogL, "\n"))
#	print( params$sigma)
#	print( params$beta)
}


	return( params ) 
}

gaussian.fn <- function( p, d=NULL ) {
	
	params <- list( beta=p[2:length(p)], sigma2=p[1] )
	params.new <- gaussian.iterate( d, params )
	res <- params.new$LogL
	attr(res, "gradient") <- c( params.new$dsigma2, params.new$dbeta ) 
	res
}

gaussian.gr <-  function( p, d=NULL ) {
	
	params <- list( beta=p[2:length(p)], sigma2=p[1] )
	params.new <- gaussian.iterate( d, params )
	c( params.new$dsigma2, params.new$dbeta ) 
}



gfit <- function( h, eps=1.0e-4, shuffle=FALSE, method="optim" ) {

	y <- h$phenotypes
	nm <- length(h$markers)-1
	table <- matrix(nrow=nm,ncol=7)
	colnames(table) <- c( "cM", "marker", "LogP", "ChiSq",  "Null", "df", "Pval")

	if ( shuffle ) y <- sample(y)
	for( m in 1:nm ) {
		p <- hprob( h, h$markers[m] ) 
		q <- qr(p)
		df <- q$rank
		d <- list( y=y, probs=p )
		if ( method == "optim" ) {
		        params0 <- gaussian.init( d )
		        p0 <- c( params0$sigma2, params0$beta )
#			res <- nlm( gaussian.fn, p0, print.level=2, check.analyticals=TRUE, d=d)
			res <- optim( p0, gaussian.fn, gaussian.gr, method="BFGS",  d=d)
			chi <- 2*(params0$LogL - res$value)
			pval <- pchisq( chi, df , lower.tail=FALSE)
			LogPval <- -log10(pval)
			table[m,] <- c( h$map[m], h$markers[m], LogPval, chi, params0$LogL, df, pval )
			print(table[m,])
		}
		else {
			res <- gaussian.loop( d, eps=eps, df=df )
			cat( h$markers[m], "chi",  res$chi, "df", res$df, "logPval", res$LogPval, "LogL", res$LogL, "Null", res$Null, "\n" )
			table[m,] <- c( h$map[m], h$markers[m], res$LogPval, res$chi, res$Null, res$df, res$Pval )
		}
	}
	return ( list( table=table, offset=3, width=1, model="mixture",test="gfit", method=method, maxm=max(as.numeric(table[,3]),na.rm=TRUE), maxp=which.max(as.numeric(table[,3])) ))
}


# Genome Cache functions

save.genome <- function ( gdir, sdir, prefix, chrs=NULL,
                         file.format="ped", mapfile=NULL, ancestryfile=NULL, generations=50, phase="unknown", haploid=FALSE, mc.cores=1 ) {

  
  if ( is.null(chrs) ) 
    chrs <- the.chromosomes()

  if ( ! file.exists(sdir))
    dir.create(sdir)
  if ( haploid == FALSE ) {
    full <- paste(sdir, "/full/", sep="")
    dir.create(full)
  }
  additive <- paste(sdir, "/additive/", sep="")
  dir.create(additive)
  genotype <- paste(sdir, "/genotype/", sep="")
  dir.create(genotype)

  if ( mc.cores <=1 ) {
    lapply( chrs, save.happy.internal, gdir, prefix, file.format, ancestryfile, generations, mapfile, phase, haploid, additive, full, genotype) 
  }
  else	{
     mclapply( chrs, save.happy.internal, gdir, prefix, file.format, ancestryfile, generations, mapfile, phase, haploid, additive, full, genotype, mc.cores=mc.cores) 
  }
}

save.happy.internal <- function( chr, gdir, prefix, file.format, ancestryfile, generations, mapfile, phase, haploid, additive, full, genotype) { 
	h <- happy( paste( gdir, chr, prefix, ".data", sep="" ),
        	       paste( gdir, chr, prefix, ".alleles", sep="" ),
               	file.format=file.format,
               	ancestryfile=ancestryfile,
	        generations=generations, do.dp=TRUE, mapfile=mapfile, phase=phase, haploid=haploid )
    save.happy( h, chr, dir=additive, model="additive"  )
    if ( h$haploid == FALSE ) save.happy( h, chr, dir=full, model="full"  )
    save.happy( h, chr, dir=genotype, model="genotype"  )
    delete.happy.cobject(h)
}

delete.happy.cobject <- function(h)
{
  cat("delete.happy() called\n")
}


the.chromosomes <- function( autosomes=19, use.X=FALSE ) {

  if ( use.X) 
    return(paste( "chr", c(1:autosomes,"X"), sep=""))
  else
    return(paste( "chr", c(1:autosomes), sep=""))
}


save.happy <- function( h, pkg, dir, model="additive" ) {

  ddp <- paste( dir, "/", pkg, sep="")
  g.data.attach(ddp)

  if ( model == "genotype" )
    nm <- length(h$markers)
  else
    nm <- length(h$markers) -1
#  markers.safe = as.character(h$markers[1:nm])
  markers.safe = make.names(h$markers[1:nm])


  assign("markers", h$markers[1:nm], 2)
  assign("markers.safe",markers.safe,2)
  assign("map", h$map[1:nm], 2 )
  assign("chromosome", h$chromosome[1:nm], 2 )
  assign("subjects", h$subjects, 2 );
  assign("strains", h$strains, 2 )
  assign("haploid", h$haploid, 2)
  print ("saving strains")
  if ( !is.null(h$bp))
    assign("bp", h$bp[1:nm], 2);
  
  if ( model == "genotype" ) {
    for( m in 1:nm ) {
      assign(markers.safe[m], hgenotype( h, m, collapse=FALSE ), 2)
    }
    g.data.save(ddp)
  }
  else {
    for( m in 1:nm ) {
      assign(markers.safe[m], hdesign( h, m, model=model ), 2)
    }
    g.data.save(ddp)
  }
  return(ddp)
}

load.genome <- function (dir,
                         use.X = TRUE,
                         chr   = the.chromosomes(use.X=use.X),
                         n.chr=NA,
                         models=c("additive", "full", "genotype"))  # CHANGED
{
  g <- list()
  old.subjects <- NULL
  old.strains <- NULL
  if ( is.integer(n.chr) )
    chr = paste("chr", 1:chr, sep="")
  
  for (model in models) {
    if ( file.exists( paste(dir, model, sep = "/") )) {
    pkgs <- paste(dir, model, chr, sep = "/")   # CHANGED
    pkgs = pkgs[file.exists(pkgs)]
    markers <- c()
    chromosome <- c()
    map <- c()
    pkgname <- c()
    bp <- c()
    for (p in pkgs) {
#      chromosome <- c(chromosome, g.data.get("chromosome", p))
      chromosome <- c(chromosome, happy.load.data("chromosome", p))
      m <- happy.load.data("markers", p)
      markers <- c(markers, m)
      map <- c(map, happy.load.data("map", p))
      bp <- c(bp, happy.load.data("bp", p))
      pkgname <- c(pkgname, rep(p, length(m)))
      subjects <- happy.load.data("subjects", p)
      strains <- happy.load.data("strains", p)

      if ( is.null(old.subjects)) {
        old.subjects <- subjects
      }
      if ( !identical(subjects,old.subjects )) {
        cat( "ERROR - subject names are inconsistent for chromosome ", chromosome[1] , "\n", subjects, "\n", old.subjects, "\n")
        stop( "FATAL HAPPY ERROR")
      }

      if ( is.null(old.strains)) {
        old.strains <- strains
      }
      if ( ! identical( strains, old.strains) ) {
        cat( "ERROR - strain names are inconsistent for chromosome ", chromosome[1] , "\n", strains, "\n", old.strains, "\n")
        stop( "FATAL HAPPY ERROR")
      }
    }
    genome <- data.frame(
                         marker     = I(as.character(markers)),
                         map        = as.numeric(map),
                         bp         = as.numeric(bp),
                         ddp        = I(as.character(pkgname)),
                         chromosome = I(as.character(chromosome)))
    g[[model]] <- list(
                       genome     = genome,
                       subjects   = subjects,
                       strains    = strains,
                       markers    = as.character(genome$marker),
                       chromosome = as.character(genome$chromosome),
                       map        = genome$map)
  }
  }
  g$subjects <- g$genotype$subjects
  g$strains  <- g$additive$strains
  g$markers  <- g$genotype$markers
  g$chromosome <- g$genotype$chromosome
  g$map <- g$genotype$map
  class(g)   <- "happy.genome"
  return(g)
} 




load.markers <- function( genome, markers, model="additive", include.models=FALSE ) {

  if ( length(model) == 1 )
    model <- rep( model, length(markers))

  marker.list <- list()
  model.list <- list()

  for(i in 1:length(markers)) {
    
    genome.model <- genome[[model[i]]]
    if ( is.numeric(markers))
      rows <- markers
    else
      rows <- pmatch( as.character(markers[i]),
                     as.character(genome.model$genome[,1]), nomatch=NA )
  
    if ( length(rows) > 0 ) {
      r <- rows[1]
      
      m <- as.character(genome.model$genome[r,"marker"])
#      m.names = make.names(m)
      
      pkg <- as.character(genome.model$genome[r,"ddp"])
      marker.list[[m]] <- happy.load.data( m, pkg) ###
      model.list[[m]] <- model[i]
  
    }
  }
  if ( include.models )
    return( list( marker=marker.list, model=model.list ))
  else
    return( marker.list )
}

happy.load.data <- function (item, dir) # replaces calls to g.data.get, to make things backwards compatible.
{
    env <- new.env()

    # determine which version of g.data was used to save the data
    filename.pre2009  <- file.path(dir, "data", paste(item, "RData", sep = "."))
    if (file.exists(filename.pre2009))
    {
        load(filename.pre2009, env)
        return ( get(item, envir = env) )
    }
    
    # assume 2009 version of g.data was used
    filename.post2009 <- file.path(dir,
            paste(gsub("([[:upper:]])", "@\\1", item), "RData", sep = ".")
            )
    if (file.exists(filename.post2009))
    {       
        load(filename.post2009, env)
        return ( get(item, envir = env ) )
    }
    mm = make.names(item)
    filename.make.names = file.path(dir,paste(gsub("([[:upper:]])", "@\\1", mm), "RData", sep = "."))

    if (file.exists(filename.make.names))
    {       
        load(filename.make.names, env)
        return ( get(mm, envir = env ) )
    }
    stop("Could not find data for ", item, " in package ", dir)
}




