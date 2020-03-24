library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
library(tidyr)

# Data prep ---------------------------------------------------------------

# set whether or not to truncate to standard global temperature range
doTrunc <- TRUE

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

source('GSA_custom_ecospat_fcns.R')

# if not running any sections above, load the data:
# TODO rename outDf
if (doTrunc){
  outDf <- read.csv('Data/foram_MAT_occs_latlong_8ka_trunc_20-03-24.csv',
                    stringsAsFactors = FALSE)
  samp <- read.csv('Data/samp_MAT_occs_latlong_8ka_trunc_20-03-24.csv',
                   stringsAsFactors = FALSE)
} else {
  outDf <- read.csv('Data/foram_MAT_occs_latlong_8ka_20-03-24.csv',
                    stringsAsFactors = FALSE)
  samp <- read.csv('Data/samp_MAT_occs_latlong_8ka_20-03-24.csv',
                   stringsAsFactors = FALSE)
}

spp <- unique(outDf$species)
bins <- unique(outDf$bin)
nbins <- length(bins)
envCol <- c('temp_ym_hab','temp_ym_0m')
ncores <- detectCores() - 1

# KDE niche summary -------------------------------------------------------

# Output niche overlap (though time), peak abundance, preferred enviro, & tolerance
nicher <- function(dat, b1, b2, s, env, xmn, xmx,
                   w1 = NULL, w2 = NULL, reflect = FALSE, ...){
  
  sp1rows <- which(dat$species==s & dat$bin==b1)
  sp1 <- dat[sp1rows,env]
  
  sp2rows <- which(dat$species==s & dat$bin==b2)
  sp2 <- dat[sp2rows,env]
  
  noWeight <- any(is.null(w1), is.null(w2))
  if (! noWeight){
    d1 <- tryCatch(
      transformEst(sp1, w = w1, bw='SJ-ste', reflect = reflect, a = xmn, b = xmx, ...),
#      JonesEst(sp1, w = w1, bw='brt', reflect = reflect, a = xmn, b = xmx, ...),
      error = function(err){ list() }
    ) 
    d2 <- tryCatch(
      transformEst(sp2, w = w2, bw='SJ-ste', reflect = reflect, a = xmn, b = xmx, ...),
#      JonesEst(sp2, w = w2, bw='brt', reflect = reflect, a = xmn, b = xmx, ...),
      error = function(err){ list() }
    ) 
  } else {
    # use unweighted sampling, either regular density or reflected:
    if (reflect){
      d1 <- tryCatch(
        density.reflected(sp1, bw='SJ-ste', lower = xmn, upper = xmx, ...),
      )
      d2 <- tryCatch(
        density.reflected(sp2, bw='SJ-ste', lower = xmn, upper = xmx, ...),
      )
    } else {
      d1 <- tryCatch(
        density(sp1, bw='SJ-ste', from = xmn, to = xmx, ...),
      )
      d2 <- tryCatch(
        density(sp2, bw='SJ-ste', from = xmn, to = xmx, ...),
      )
    }
  }
  
  # the species may be absent in one or both bins, in which case d is an empty list
  if (length(d1)==0){
    data.frame(bin=NA, sp=NA, h=NA, pa=NA, pe=NA)
    
  } else {
    stats <- nichStats(d1)
    
    if (length(d2)==0){
      data.frame(bin=b1, sp=s, h=NA, t(stats))
    } else{
      h <- hell(d1, d2, extrap = TRUE) 
      data.frame(bin=b1, sp=s, h=h, t(stats))
    }
  }
}

# the older bin is column 1, the younger is column 2
bPairs <- cbind(bins[-1], bins[-length(bins)])
# for the most recent time bin, it's not possible to calculate overlap
# (because no modern data are included), but include it anyway
# so that the standing niche at the last time bin is calculated
recent <- cbind(4, NA)
bPairs <- rbind(recent, bPairs)

# loop over time bins over species over environmental variables
kdeLoop <- function(e,dat){
  xmn <- min(samp[,e])
  xmx <- max(samp[,e])
  
  bList <- apply(bPairs, 1, function(x){
    # estimate bias function for each time bin
    # based on sampling distribution, with boundary reflection
    sampRows1 <- which(samp$b==x[1])
    samp1 <- samp[sampRows1,e]
    densSamp1 <- density(samp1, bw='SJ-ste', from=xmn, to=xmx)
#    densSamp1 <- density.reflected(samp1, from=xmn, to=xmx) # typo - should be lower/upper not from/to
    w1 <- approxfun(densSamp1$x, densSamp1$y)
    
    # in the most recent time bin, there is no subsequent bin
    if (is.na(x[2])){
      w2 <- NA
    } else {
      sampRows2 <- which(samp$b==x[2])
      samp2 <- samp[sampRows2,e]
      densSamp2 <- density(samp2, bw='SJ-ste', from=xmn, to=xmx)
#      densSamp2 <- density.reflected(samp2, from=xmn, to=xmx) # typo - should be lower/upper not from/to
      w2 <- approxfun(densSamp2$x, densSamp2$y)
    } 
    
    sList <- lapply(spp, function(s){
      nicher(dat = dat, b1 = x[1], b2 = x[2], s = s, env = e, 
             w1 = w1, w2 = w2, xmn = xmn, xmx = xmx)
    })
    do.call(rbind, sList)
  })
  sDf <- do.call(rbind, bList)
  # remove NA rows (if a species is not sampled in the focal bin)
  nas <- is.na(sDf$bin)
  sDf[!nas,]
}

# Note: might be weird/problematic to truncate based on surface temp,
# but then use in-habitat temp for niches. Temps at depth can get colder
# than the standardised lower bound from surface data.
# Perhaps an alternative to surface and in-habitat temp is to use MLD temp.
kdeSum <- kdeLoop('temp_ym_0m', outDf)
  # eList <- lapply(envCol, kdeLoop, dat=outDf)
  # kdeSum <- merge(eList[[1]], eList[[2]], by=c('sp','bin'), no.dups=TRUE, 
  #              suffixes=paste0('_', envCol))

# Non-KDE niche summary ---------------------------------------------------

# Need mean, variance, sample size, and age of trait values (MAT) for each sp & bin
sumup <- function(bin, s, dat, binCol, sCol, traitCol){
  slcRows <- which(dat[,binCol] == bin & dat[,sCol] == s)
  if (length(slcRows)>0){
    slc <- dat[slcRows,]
    x <- slc[,traitCol]
    if (length(traitCol)==1){
      m <- mean(x)
      sd <- sd(x)
      n <- length(x)
      rtrn <- data.frame(bin=bin, sp=s, m=m, sd=sd, n=n)
    } else {
      m <- apply(x, 2, mean)
      sd <- apply(x, 2, sd)
      n <- nrow(x)
      rtrn <- data.frame(cbind(bin, t(m), t(sd), n))
      colnames(rtrn) <- c('bin', paste0('m_', traitCol), paste0('sd_', traitCol), 'n')
      rtrn$sp <- s
    }
    
  } else {
    rtrn <- data.frame()
  }
  return(rtrn)
}

rawSumL <- lapply(spp, function(s){
  temp <- lapply(bins, sumup, s=s, dat=outDf, 
                 binCol='bin', sCol='species', traitCol=envCol)
  tempDf <- do.call(rbind, temp)
})
rawSum <- do.call(rbind, rawSumL)

# combine KDE and raw-scale summary/sample statistics into one output file
fullSum <- merge(kdeSum, rawSum, all.x=TRUE)
if (doTrunc){
  sumNm <- paste0('Data/foram_niche_sumry_metrics_trunc_',day,'.csv') 
} else {
  sumNm <- paste0('Data/foram_niche_sumry_metrics_',day,'.csv')
}
write.csv(fullSum, sumNm, row.names=FALSE)

# Inter-specific overlap --------------------------------------------------

interSppD <- function(b, df, env){
  xmx <- max(df[,env])
  xmn <- min(df[,env])
  
  bRows <- which(df$bin==b)
  bSpp <- unique(df$species[bRows])
  
  bSampRows <- which(samp[,'b']==b)
  bSamp <- samp[bSampRows,env]
  sampDens <- density.reflected(bSamp, from=xmn, to=xmx)
  w <- approxfun(sampDens$x, sampDens$y)
  
  # Construct KDE of all species
  kdeL <- lapply(bSpp, function(s){
    spRows <- which(df$species==s & df$bin==b)
    sp1 <- df[spRows,env]
    d <- JonesEst(sp1, bw='brt', w = w, from = xmn, to = xmx, reflect = TRUE)
  }
  )
  names(kdeL) <- bSpp
  
  # Compute Hellinger's H for all species pairs
  fin <- data.frame()
  for (s1 in bSpp){
    for (s2 in bSpp){
      if (s1==s2) {next} else{
        h <- hell(kdeL[[s1]], kdeL[[s2]])
        pairDat <- data.frame(envVar=env, bin=b, sp1=s1, sp2=s2, h=h)
        fin <- rbind(fin, pairDat)
      }
    }
  }
  fin
}

registerDoParallel(ncores)
interLong <- foreach(bin=bins, .packages='pracma', .combine=rbind, .inorder=FALSE) %dopar% {
    interSppD(b=bin, df=outDf, env='temp_ym_0m')
  }
stopImplicitCluster()

# note: if only using surface data, this could be streamlined; doesn't need to be spread
interWide <- spread(interLong, envVar, h)
colnames(interWide)[ncol(interWide)] <- 'h_temp_ym_0m'
if (doTrunc){
  interSppNm <- paste0('Data/foram_species_pairs_KDE_H_trunc_', day, '.csv')
} else {
  interSppNm <- paste0('Data/foram_species_pairs_KDE_H_', day, '.csv')
}
write.csv(interWide, interSppNm, row.names=FALSE)
