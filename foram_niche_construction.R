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

source('species_kde_buildr.R')

df <- read.csv('Data/foram_uniq_occs_latlong_8ka_20-07-16.csv', 
               stringsAsFactors = FALSE)
samp <- read.csv('Data/samp_uniq_occs_latlong_8ka_20-07-16.csv',
                 stringsAsFactors = FALSE)
spAttr <- read.csv('Data/foram_spp_data_20-07-16.csv',
                   stringsAsFactors = FALSE)

envNm <- 'temp_ym'
envCols <- grep(envNm, colnames(df))
allEnvNm <- colnames(df)[envCols]
df <- df[,c('species', 'bin', allEnvNm)] # 'cell_number','coreUniq',
colnames(samp)[1:2] <- c('bin','cell_number')
spp <- unique(df$species)
bins <- unique(df$bin)
nCore <- detectCores() - 1

# Truncate to standard temp range -----------------------------------------

# Restrict to the last 700 ka, to trim edge effects of range-through cores
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
  prnt <- c(round(uppr,1), round(lwr, 1))
  print(prnt)
  
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
    toss <- c(tooBig, tooSmol)
    if (length(toss) > 0){
      slc <- slc[-toss,]
    }
    trunc <- rbind(trunc, slc)
    
    bBoolSamp <- samp$bin==b
    slcSamp <- samp[bBoolSamp,]
    tooBigSamp <- which(slcSamp[,e] > uppr)
    tooSmolSamp <- which(slcSamp[,e] < lwr)
    sampToss <- c(tooBigSamp, tooSmolSamp)
    if (length(sampToss) > 0){
      slcSamp <- slcSamp[-sampToss,]
    }
    truncSamp <- rbind(truncSamp, slcSamp)
  }
  
  # The last steps could introduce species with <6 occs
  tossRowsL <- 
    sapply(spp, function(x){
      sapply(bins, function(b){
        tooRare(sp=x, bin=b, df=trunc)
      } )
    } )
  tossRows <- unlist(tossRowsL)
  trunc <- trunc[-tossRows,]
  
  # Check that the species is sampled consecutively at least once 
  keepSpp <- character()
  spp <- unique(trunc$species)
  binL <- bins[2] - bins[1]
  enufTxt <- paste(binL)
  # Alternatively, as for time series analysis, ensure samplng
  # in at least 7 successive steps i.e. 6 boundary crossings
    #  enuf <- rep(binL, 6)
    # enufTxt <- paste0(enuf, collapse='')
  for (s in spp){
    spBool <- trunc$species==s
    spDf <- trunc[spBool,]
    spB <- sort(unique(spDf$bin))
    bDiff <- diff(spB)
    diffTxt <- paste0(bDiff, collapse='')
    srch <- grep(enufTxt, diffTxt)
    if (length(srch) > 0){
      keepSpp <- c(keepSpp, s)
    }
  }
  keepBool <- trunc$species %in% keepSpp
  trunc <- trunc[keepBool,]
  
  # retain the columns necessary and sufficient for KDE
  trunc <- trunc[, c('species', 'bin', e)]
  colnames(trunc)[3] <- envNm
  truncSamp <- truncSamp[, c('bin', e)]
  colnames(truncSamp)[2] <- envNm
  
  list(sp = trunc, samp = truncSamp)
  # output is a list of two dataframes (species-level data & sampling data)
}

truncEnv <- lapply(allEnvNm, truncatr)
names(truncEnv) <- allEnvNm
nrow(truncEnv$temp_ym_0m$sp)
truncNm <- paste0('Data/sampled_',envNm,'_truncated_by_depth_',day,'.rds')
saveRDS(truncEnv, truncNm)

# maximum temperature cutoff for each depth
# > [1] 27.1 -1.4
# > [1] 26.9 -1.3
# > [1] 25.9 -0.6
# > [1] 23.8  0.8

# Evaluate degree of truncation -------------------------------------------

df$trunc <- 'in range'
tooBig <- which(df$temp_ym_0m > 27.1)
tooSmol <- which(df$temp_ym_0m < -1.4)
df$trunc[tooBig] <- 'high'
df$trunc[tooSmol] <- 'low'

# omit the most recent 2 time bins because the data are too numerous
# too plot on the same scale
mdrnBool <- df$bin %in% bins[1:2]
mdrn <- df[mdrnBool,]
table(mdrn$bin)
old <- df[!mdrnBool,]
old$trunc <- factor(old$trunc, levels = c('high','low','in range'))
bars <- ggplot(data = old, aes(fill = trunc, x= - bin)) + 
  scale_x_continuous(name='Time (Ka)', labels = c(600, 400, 200),
                     breaks = c(-600, -400, -200), expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 580)) +
  theme_bw() +
  geom_bar(position="stack", width = 5) +
  scale_fill_manual(name='MAT value in relation to cutoffs', 
                    values=c('plum','gold','grey20')) +
  theme(legend.position = 'top')

barNm <- paste0('Figs/truncated-data-sample-size_bars_',day,'.pdf')
pdf(barNm, width=6, height=4)
print(bars)
dev.off()

# KDE niche summary -------------------------------------------------------

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
pkg <- c('pracma','GoFKernel','kerneval')
pt1 <- proc.time()
registerDoParallel(nCore)
kdeSum <- foreach(dat=truncEnv[2:4], .combine=rbind, .inorder=FALSE, .packages=pkg) %:% 
  foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar% 
   kde(dat, bPair, envNm)
kdeSumSS <- foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar% 
  kde(truncEnv[[1]], bPair, envNm)
stopImplicitCluster()
pt2 <- proc.time()
pt2-pt1

 # remove NA rows (if a species is not sampled in the focal bin)
nas <- is.na(kdeSum$bin)
kdeSum <- kdeSum[!nas,]
naSS <- is.na(kdeSumSS$bin)
kdeSumSS <- kdeSumSS[!naSS,]

# Non-KDE niche summary ---------------------------------------------------

# Need sample size, mean, and variance of trait values for each sp & bin
sumup <- function(bin, s, dat, binCol, sCol, traitCol){
  slcRows <- which(dat$sp[,binCol] == bin & dat$sp[,sCol] == s)
  if (length(slcRows) > 0){
    x <- dat$sp[slcRows,envNm]
    m <- mean(x)
    sd <- sd(x)
    n <- length(x)
    rtrn <- data.frame(bin=bin, sp=s, m=m, sd=sd, n=n)
  } else {
    # some species do not occur in every bin
    rtrn <- data.frame()
  }
  return(rtrn)
}

# iterate over bins over species over depth habitats
envL <- lapply(truncEnv[2:4], function(dat){
  spL <- lapply(spp, function(s){
    binL <- lapply(bins, sumup, s = s, dat = dat, 
                   binCol = 'bin', sCol = 'species', traitCol = envNm)
    binDf <- do.call(rbind, binL)
  })
  spDf <- do.call(rbind, spL)
})
rawSum <- do.call(rbind, envL)

# iterate over bins over species, sea surface values only
envLss <- lapply(spp, function(s){
  binL <- lapply(bins, sumup, s = s, dat = truncEnv[[1]], 
                 binCol = 'bin', sCol = 'species', traitCol = envNm)
  binDf <- do.call(rbind, binL)
})
rawSumSS <- do.call(rbind, envLss)

# combine KDE and finite sample statistics into same output
fullSum <- merge(kdeSum, rawSum, all.x=TRUE)
sumNm <- paste0('Data/foram_niche_sumry_metrics_hab',day,'.csv')
write.csv(fullSum, sumNm, row.names=FALSE)
fullSumSS <- merge(kdeSumSS, rawSumSS, all.x=TRUE)
sumSSnm <- paste0('Data/foram_niche_sumry_metrics_0m_',day,'.csv')
write.csv(fullSumSS, sumSSnm, row.names=FALSE)

# export summary for sampling environment too (all bins and depths)
bsum <- function(d){
  bMat <- sapply(bins, function(b){
    bBool <- d$samp$bin==b
    slc <- d$samp[bBool,]
    m <- mean(slc$temp_ym)
    sdv <- sd(slc$temp_ym)
    nsite <- nrow(slc)
    cbind(b, m, sdv, nsite)
  })
  bDf <- data.frame(t(bMat))
  colnames(bDf) <- c('bin','m','sd','nsite')
  bDf
}
sampSum <- lapply(truncEnv, bsum)
sampSumNm <- paste0('Data/sampled_',envNm,'_bin_summary_by_depth_',day,'.rds')
saveRDS(sampSum, sampSumNm)
