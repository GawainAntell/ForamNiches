library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
library(tidyr)

# Data prep ---------------------------------------------------------------

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

source('GSA_custom_ecospat_fcns.R')

df <- read.csv('Data/foram_uniq_occs_latlong_8ka_20-03-27.csv', 
               stringsAsFactors = FALSE)
samp <- read.csv('Data/samp_uniq_occs_latlong_8ka_20-03-27.csv',
                 stringsAsFactors = FALSE)
spAttr <- read.csv('Data/foram_spp_data_20-03-26.csv',
                   stringsAsFactors = FALSE)

envNm <- 'temp_ym'
envCols <- grep(envNm, colnames(df))
allEnvNm <- colnames(df)[envCols]
df <- df[,c('species','bin',allEnvNm)] # 'cell_number','coreUniq',
colnames(samp)[1:2] <- c('bin','cell_number')
spp <- unique(df$species)
bins <- unique(df$bin)
nCore <- detectCores() - 1

# Truncate to standard temp range -----------------------------------------

# Restrict to the last 700 ka, encompassing 7 glacial/interglacial cycles
trimBool <- df$bin <= 700
df <- df[trimBool,]
sampTrim <- samp$bin <= 700
samp <- samp[sampTrim,]
bins <- bins[bins <= 700]

# Find the temperature thresholds for a given depth zone
minmax <- function(df, b, env){
  bBool <- df[,'bin']==b
  slc <- df[bBool,]
  range(slc[,env])
}

# Subset such that each sp has at least 6 occs per bin.
tooRare <- function(sp, bin, df){
  spRows <- which(df$species==sp & df$bin==bin)
  if (length(spRows)<6){
    spRows
  }
}  

# Convert an environmental variable name to habitat name
e2zone <- function(txt){ 
  switch(txt, 
         temp_ym_0m=c('Surface','Surface.subsurface','Subsurface'), 
         temp_ym_surf='Surface',
         temp_ym_surfsub='Surface.subsurface',
         temp_ym_sub='Subsurface')
}

truncatr <- function(e){
  sampSmry <- sapply(bins, minmax, df=samp, env=e)
  uppr <- min(sampSmry[2,])
  lwr <- max(sampSmry[1,])
  
  # Consider only the species within the focal habitat zone
  zone <- e2zone(e)
  zonePos <- spAttr$DepthHabitat %in% zone
  zoneSp <- spAttr$species[zonePos]
  zoneDfPos <- df$species %in% zoneSp
  zoneDf <- df[zoneDfPos,]
  
  trunc <- truncSamp <- data.frame()
  for (b in bins){
    bBool <- zoneDf$bin==b
    slc <- zoneDf[bBool,]
    tooBig <- which(slc[,e] > uppr)
    tooSmol <- which(slc[,e] < lwr)
    out <- c(tooBig, tooSmol)
    if (length(out) > 0){
      slc <- slc[-out,]
    }
    trunc <- rbind(trunc, slc)
    
    bBoolSamp <- samp$bin==b
    slcSamp <- samp[bBoolSamp,]
    tooBig <- which(slcSamp[,e] > uppr)
    tooSmol <- which(slcSamp[,e] < lwr)
    out <- c(tooBig, tooSmol)
    if (length(out) > 0){
      slcSamp <- slcSamp[-out,]
    }
    truncSamp <- rbind(truncSamp, slcSamp)
  }
  
  # The last steps could introduce species with <6 occs.
  tossRowsL <- 
    sapply(spp, function(x){
      sapply(bins, function(b){
        tooRare(sp=x, bin=b, df=trunc)
      } )
    } )
  tossRows <- unlist(tossRowsL)
  trunc <- trunc[-tossRows,]
  
  # Also re-check for per-species continuity through time 
  # (at least 7 successive steps of 8ka i.e. 6 boundary crossings)
  keepSpp <- character()
  spp <- unique(trunc$species)
  binL <- bins[2] - bins[1]
  enuf <- rep(binL, 6)
  enufTxt <- paste0(enuf, collapse='')
  for (s in spp){
    spBool <- trunc$species==s
    spDf <- trunc[spBool,]
    spB <- sort(unique(spDf$bin))
    bDiff <- diff(spB)
    diffTxt <- paste0(bDiff, collapse='')
    srch <- grep(enufTxt,diffTxt)
    if (length(srch) > 0){
      keepSpp <- c(keepSpp, s)
    }
  }
  keepBool <- trunc$species %in% keepSpp
  trunc <- trunc[keepBool,]
  
  # retain the columns necessary and sufficient for KDE
  trunc <- trunc[,c('species','bin',e)]
  colnames(trunc)[3] <- envNm
  truncSamp <- truncSamp[,c('bin',e)]
  colnames(truncSamp)[2] <- envNm
  
  list(sp=trunc, samp=truncSamp)
  # output is a list of two dataframes (species-level data & sampling data)
}

truncEnv <- lapply(allEnvNm, truncatr)
names(truncEnv) <- allEnvNm

# KDE niche summary -------------------------------------------------------

# Output 1 sp's niche overlap (though time), peak abundance, & preferred enviro
nicher <- function(dat, b1, b2, s, env, xmn, xmx,
                   w1 = NULL, w2 = NULL, reflect = FALSE, ...){
  
  sp1rows <- which(dat$sp$species==s & dat$sp$bin==b1)
  sp1 <- dat$sp[sp1rows,envNm]
  
  sp2rows <- which(dat$sp$species==s & dat$sp$bin==b2)
  sp2 <- dat$sp[sp2rows,envNm]
  
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

kde <- function(dat, bPair){
  b1 <- bPair[1]
  b2 <- bPair[2]
  xmn <- min(dat$samp[,envNm])
  xmx <- max(dat$samp[,envNm])
  
  # estimate bias function for each time bin based on sampling distribution
  sampRows1 <- which(dat$samp$bin==b1)
  samp1 <- dat$samp[sampRows1,envNm]
  # Reflecting the sample curve doesn't change the sp KDE much
  # except that the ends turn down a bit more (more convexity).
  # Since it's more complicated and throws warnings, don't do it.
  densSamp1 <- density(samp1, bw='SJ-ste', from=xmn, to=xmx)
  #  densSamp1 <- density.reflected(samp1, lower=xmn, upper=xmx) 
  w1 <- approxfun(densSamp1$x, densSamp1$y)
  
  # in the most recent time bin, there is no subsequent bin
  if (is.na(b2)){
    w2 <- NA
    
    # redefine axis limits - discretization of samp KDE can shrink them a bit
    xmnNew <- min(densSamp1$x)
    xmxNew <- max(densSamp1$x)
  } else {
    sampRows2 <- which(dat$samp$bin==b2)
    samp2 <- dat$samp[sampRows2,envNm]
    densSamp2 <- density(samp2, bw='SJ-ste', from=xmn, to=xmx)
    #  densSamp2 <- density.reflected(samp2, lower=xmn, upper=xmx) 
    w2 <- approxfun(densSamp2$x, densSamp2$y)
    
    # redefine axis limits - discretization of samp KDE can shrink them a bit
    lim1 <- range(densSamp1$x)
    lim2 <- range(densSamp2$x)
    xmnNew <- max(lim1[1], lim2[1])
    xmxNew <- min(lim1[2], lim2[2])
  } 
  if (sum(identical(xmn, xmnNew), identical(xmx, xmxNew))!=2){
    error('code needed after all')
    }
  zoneSp <- unique(dat$sp$species)
  sList <- lapply(zoneSp, function(s){
    nicher(dat = dat, b1 = b1, b2 = b2, s = s,
           w1 = w1, w2 = w2, xmn = xmn, xmx = xmx)
  })
  do.call(rbind, sList)
}

# for the most recent time bin, it's not possible to calculate overlap
# (because no modern data are included), but include it anyway
# so that the standing niche at the last time bin is calculated
pairL <- list()
bPairs <- cbind(bins, c(NA, bins[-length(bins)]) )
for (i in 1:length(bins)){
  entry <- bPairs[i,]
  pairL <- append(pairL, list(entry))
}

# warning - this could take 1 - 2 hours
pkg <- c('pracma','GoFKernel')
pt1 <- proc.time()
registerDoParallel(nCore)
kdeSum <- foreach(dat=truncEnv[2:4], .combine=rbind, .inorder=FALSE, .packages=pkg) %:% 
  foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar% 
   kde(dat, bPair)
kdeSumSS <- foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar% 
  kde(truncEnv[[1]], bPair)
stopImplicitCluster()
pt2 <- proc.time()
pt2-pt1

# remove NA rows (if a species is not sampled in the focal bin)
nas <- is.na(kdeSum$bin)
kdeSum <- kdeSum[!nas,]
naSS <- is.na(kdeSumSS$bin)
kdesumSS <- kdeSumSS[!naSS,]

# Non-KDE niche summary ---------------------------------------------------
 
# Need mean, variance, sample size, & age of trait values for each sp & bin
sumup <- function(bin, s, dat, binCol, sCol, traitCol){
  slcRows <- which(dat$sp[,binCol] == bin & dat$sp[,sCol] == s)
  if (length(slcRows)>0){
    x <- dat$sp[slcRows,envNm]
    m <- mean(x)
    sd <- sd(x)
    n <- length(x)
    rtrn <- data.frame(bin=bin, sp=s, m=m, sd=sd, n=n)
  } else {
    rtrn <- data.frame()
  }
  return(rtrn)
}

# iterate over bins over species over depth habitats
envL <- lapply(truncEnv[2:4], function(dat){
  spL <- lapply(spp, function(s){
    binL <- lapply(bins, sumup, s=s, dat=dat, 
                   binCol='bin', sCol='species', traitCol=envNm)
    binDf <- do.call(rbind, binL)
  })
  spDf <- do.call(rbind, spL)
})
rawSum <- do.call(rbind, envL)

# iterate over bins over species, surface values only
envLss <- lapply(spp, function(s){
  binL <- lapply(bins, sumup, s=s, dat=truncEnv[[1]], 
                 binCol='bin', sCol='species', traitCol=envNm)
  binDf <- do.call(rbind, binL)
})
rawSumSS <- do.call(rbind, envLss)

# combine KDE and finite sample statistics into same output

fullSum <- merge(kdeSum, rawSum, all.x=TRUE)
sumNm <- paste0('Data/foram_niche_sumry_metrics_',day,'.csv')
write.csv(fullSum, sumNm, row.names=FALSE)

fullSumSS <- merge(kdeSumSS, rawSumSS, all.x=TRUE)
sumSSnm <- paste0('Data/foram_niche_sumry_metrics_0m_',day,'.csv')
write.csv(fullSumSS, sumSSnm, row.names=FALSE)
