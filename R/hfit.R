                                        # fit a QTL to the markers using the specified mode

                                        # happy is a happy object returned by a call to happy
                                        # markers is a array of marker names or marker indices
                                        # mode is one of 'additive', 'full' or 'partial'
                                        # verbose controls the amount of output
                                        # return value is a table giving the logP values for each marker tested

hfit <- function( h, markers=NULL, model='additive', mergematrix=NULL, covariatematrix=NULL, verbose=FALSE, phenotype=NULL, family='gaussian', permute=0, chromosome=NULL ) {

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
  } else {
    
    map <- h$map

    if ( is.null(markers) ) 
      markers <- h$markers[1:length(h$markers)-1]

    if ( is.null(chromosome) ) 
      chromosome <- h$chromosome[1:length(markers)]

    if ( is.null(phenotype) )
      phenotype = h$phenotype
  }
  
  if ( model == 'partial' || model == 'full' ) {
    lp <- matrix( ncol=8, nrow=length(map)-1)
    colnames(lp) <-c('cM', 'marker', 'additive logP', 'full logP', 'partial logP', 'additive SS', 'full SS', 'partial SS')
    offset <- 3
    width <- 3
    idx <- 3
  } else {
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
    p01<-quantile(maxlogp,probs=(1-0.01))
    p05<-quantile(maxlogp,probs=(1-0.05))
    if ( verbose ) print(quantile(maxlogp,0.01*c(0,50,75,90,95,99,99.5,1)))
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

  } else {

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
          } else {
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
          } else {
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

      } else {

        if ( ! is.null(d<- hdesign( h, m, model='additive', mergematrix=mergematrix )) ) {
          if ( is.null(covariatematrix) ) {
            cfit <- glmfit( phenotype ~ 1 , family=family)
            afit <- glmfit( phenotype ~ d , family=family)
          } else {
            cfit <- glmfit( phenotype ~ covariatematrix )
            d <- cbind( covariatematrix, d )
            afit <- glmfit( phenotype ~ d )
          }

          if ( family=="gaussian") {
            an <- anova( cfit, afit )
            logP <- -log( an[[6]])/log(10)
          } else {
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
    return(list( table=lp, model=model, chromosome=chromosome, family=family, test='hfit', offset=offset, width=width, maxp=maxp, maxm=maxm, maxSS=maxSS, permdata=NULL ))
  }
}
