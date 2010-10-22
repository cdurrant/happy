.packageName <- "happy.hbrem_2.2"


hbrem <- function( RX, HaploidInd, Ndip, Nind, Npost=2000, Nbin, Ry ) {

  brem <- .Call( "hbrem", RX, HaploidInd, Ndip, Nind, Npost, Nbin, Ry,  PACKAGE="happy.hbrem" )

  return(brem)
}

hbrem.locus <- function(m, g, model, Ry, cc, HaploidInd, Ndip, Nind, Npost, Nbin) {

    d <- hdesign(g, m, model=model)
    d.cc <- d[cc,]

    hb <- hbrem(RX=d.cc, HaploidInd=HaploidInd, Ndip=Ndip, Nind=Nind, Npost=Npost, Nbin=Nbin, Ry=Ry)

    reg.full.lm <- lm(Ry ~ d.cc)
    reg.null.lm <- lm(Ry ~ 1)
    Ftest <- anova(reg.null.lm, reg.full.lm)
    pval <- Ftest$"Pr(>F)"[2]

    cat( m, -log10(pval), "\n" )
    return( c(  -log10(pval), hb[[1]], hb[[2]], hb[[3]], hb[[4]] ))
  }


hbrem.region <- function(g, markers, Ry, cc, HaploidInd,  Npost, Nbin, mc.cores=1) {

  
  if ( HaploidInd == 0 ) {
    model = "additive"
    Ndip = length(g$strains)
    Nind = length(Ry)
  }
  else {
    model = "full"
    ns = length(g$strains)
    Ndip = ns*(ns+1)/2
    Nind = length(Ry)
  }

  nmark <- length(markers)

  
  if ( mc.cores == 1 ) 
    res =  t(sapply ( markers, hbrem.locus, g, model,Ry, cc, HaploidInd, Ndip, Nind, Npost, Nbin) )
  else {
    res.list=mclapply ( markers, hbrem.locus, g, model,Ry, cc, HaploidInd, Ndip, Nind, Npost, Nbin, mc.cores=mc.cores)
    res = t(do.call( "cbind", res.list ) )
  }


  mark.pars.df <- data.frame(res[,1:46])
  names(mark.pars.df) = c( "F.logPval", "Hbar", "sd.Ni", "BIC.qtl", "BIC.null", "BF", "logBF", "DIC.qtl", "DIC.null", "DIC.diff", "pd.qtl", "pd.null", "mode.kT", "ga", "gb", "mode.var", "med.kT", "med.mu", "med.var", "mean.kT", "mean.mu", "mean.var", "hpd.kT.50.lower", "hpd.kT.50.upper", "hpd.mu.50.lower", "hpd.mu.50.upper", "hpd.var.50.lower", "hpd.var.50.upper", "hpd.kT.75.lower", "hpd.kT.75.upper", "hpd.mu.75.lower", "hpd.mu.75.upper", "hpd.var.75.lower", "hpd.var.75.upper", "hpd.kT.95.lower", "hpd.kT.95.upper", "hpd.mu.95.lower", "hpd.mu.95.upper", "hpd.var.95.lower", "hpd.var.95.upper", "hpd.kT.99.lower", "hpd.kT.99.upper", "hpd.mu.99.lower", "hpd.mu.99.upper", "hpd.var.99.lower", "hpd.var.99.upper")
  offset = ncol(mark.pars.df)
  mark.pars.df$Name=as.character(markers)
  if ( class(g) == "happy.genome" ) {
    idx = match( markers, g[[model]]$genome$marker )
    mark.pars.df$Chr = g[[model]]$genome$chromosome[idx]
    mark.pars.df$Bp = g[[model]]$genome$bp[idx]

    bp2 = mark.pars.df$Bp[2:length(mark.pars.df$Bp)]
    bidx = which(bp2 < mark.pars.df$Bp[1:(length(mark.pars.df$Bp)-1)])-1

    mark.pars.df$CumBp = rep(0,nrow(mark.pars.df))
    mark.pars.df$CumBp[bidx+2] = mark.pars.df$Bp[bidx+1]
    mark.pars.df$CumBp = cumsum(mark.pars.df$CumBp) + mark.pars.df$Bp
         
  }


  strain.means.df <- data.frame( res[,offset:(offset+Ndip-1)])
  strain.sdevs.df <- data.frame( res[,(offset+Ndip):(offset+2*Ndip-1)])
  strain.avNis.df <- data.frame( res[,(offset+2*Ndip):(offset+3*Ndip-1)])

  if ( HaploidInd ==0 ) {
   names(strain.means.df) <- g$strains
   names(strain.sdevs.df) <- g$strains
   names(strain.avNis.df) <- g$strains
  }
  else {
   strain.names = g$strains
   num.strains = length(g$strains)
   diplotype.names <- matrix(kronecker(strain.names, strain.names, paste, sep="."), nrow=num.strains)
   names.full.symmetric <- c( diag(diplotype.names), diplotype.names[upper.tri(diplotype.names, diag=FALSE)])
   names(strain.means.df) <- names.full.symmetric
   names(strain.sdevs.df) <- names.full.symmetric
   names(strain.avNis.df) <- names.full.symmetric
  }


  hbrem.region.list <- list(Summary.Parameters=mark.pars.df, Group.Means=strain.means.df, Group.StDevs=strain.sdevs.df, Group.ExpCounts=strain.avNis.df)

  return(hbrem.region.list)
}
