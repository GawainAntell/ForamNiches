library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
library(tidyr)

# set whether to do KDE based on sea surface values or habitat depth zones
ss <- TRUE

# Data prep ---------------------------------------------------------------

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

source('GSA_custom_ecospat_fcns.R')

df <- read.csv('Data/foram_uniq_occs_latlong_8ka_20-03-27.csv', 
               stringsAsFactors = FALSE)
samp <- read.csv('Data/samp_uniq_occs_latlong_8ka_20-03-27.csv',
                 stringsAsFactors = FALSE)

envNm <- 'temp_ym'
envCols <- grep(envNm, colnames(df))
allEnvNm <- colnames(df)[envCols]
df <- df[,c('species','bin','cell_number','coreUniq',allEnvNm)]
# colnames(samp)[2] <- 'cell'

spp <- unique(df$species)
bins <- unique(df$bin)
nCore <- detectCores() - 1

# Truncate to standard temp range -----------------------------------------

# Truncate series to the last 700 ka, encompassing 7 glacial/interglacial cycles
trimBool <- df$bin <= 700
df <- df[trimBool,]
sampTrim <- samp[,'b'] <= 700
samp <- samp[sampTrim,]
bins <- bins[bins <= 700]

minmax <- function(df, b, env){
  bBool <- df[,'b']==b
  slc <- df[bBool,]
  rng <- range(slc[,env])
  c(b, rng)
}

sampSmryM <- sapply(bins, minmax, df=samp, env='temp_ym_0m')
sampSmry <- data.frame(t(sampSmryM))
colnames(sampSmry) <- c('bin','min','max')
uppr <- min(sampSmry$max)
lwr <- max(sampSmry$min)

p <- ggplot(data=sampSmry) + theme_bw() +
  scale_x_continuous(name='Time (ka)', expand=c(0.01,0)) +
  scale_y_continuous(name = 'MAT (degrees C)') +
  geom_linerange(aes(x=-bin, ymin=min, ymax=max), colour='red') +
  geom_linerange(aes(x=-bin, ymin=lwr, ymax=uppr), colour='black') +
  geom_hline(yintercept=uppr, colour='grey', lwd=1) +
  geom_hline(yintercept=lwr, colour='grey', lwd=1)

pNm <- paste0('Figs/standardised_MAT_max_min_', day, '.pdf')
pdf(pNm, width = 6, height=4)
print(p)
dev.off()

trunc <- data.frame()
for (b in bins){
  bBool <- df$bin==b
  slc <- df[bBool,]
  tooBig <- which(slc$temp_ym_0m > uppr)
  tooSmol <- which(slc$temp_ym_0m < lwr)
  out <- c(tooBig, tooSmol)
  if (length(out) > 0){
    slc <- slc[-out,]
  }
  trunc <- rbind(trunc, slc)
}

# apply truncation for sampled site data too
truncSamp <- data.frame()
for (b in bins){
  bBool <- samp[,'b']==b
  slc <- samp[bBool,]
  tooBig <- which(slc[,'temp_ym_0m'] > uppr)
  tooSmol <- which(slc[,'temp_ym_0m'] < lwr)
  out <- c(tooBig, tooSmol)
  if (length(out) > 0){
    slc <- slc[-out,]
  }
  truncSamp <- rbind(truncSamp, slc)
}

# * Evaluate degree of truncation -----------------------------------------

df$trunc <- 'in range'
tooBig <- which(df$temp_ym_0m > uppr)
tooSmol <- which(df$temp_ym_0m < lwr)
df$trunc[tooBig] <- 'high'
df$trunc[tooSmol] <- 'low'

mdrnBool <- df$bin %in% bins[1:2]
mdrn <- df[mdrnBool,]
old <- df[!mdrnBool,]
old$trunc <- factor(old$trunc, levels=c('high','low','in range'))
bars <- ggplot(data=old, aes(fill=trunc, x= - bin)) + 
  scale_x_continuous(name='time (excluding last 16 ka)', expand=c(0.01,0)) +
  scale_y_continuous(expand=c(0,0)) +
  #  theme_bw() +
  geom_bar(position="stack", width = 5) +
  scale_fill_manual(name='MAT value in relation to cutoffs', 
                    values=c('plum','gold','grey20')) +
  theme(legend.position = 'top')

barNm <- paste0('Figs/truncated_data_sample_size_',day,'.pdf')
pdf(barNm, width=6, height=4)
print(bars)
dev.off()

# inspect the proportion of observations remaining
nrow(trunc)/nrow(df) # all data
table(old$trunc)['in range']/nrow(old) # excluding most recent 16 ka

# Clean -------------------------------------------------------------------

# The last steps could introduce species with <6 occs.
# Subset again such that each sp has >5 occs per bin.
tooRare <- function(sp, bin, df){
  spRows <- which(df$species==sp & df$bin==bin)
  if (length(spRows)<6){
    spRows
  }
}  
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

# check for any gaps in the time series
binsObs <- sort(unique(trunc$bin))
if (any(diff(binsObs) != binL)) warning('discontinuous time bins')

# KDE niche summary -------------------------------------------------------

# Output 1 sp's niche overlap (though time), peak abundance, & preferred enviro
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

kde <- function(e, dat, bPair, xmn, xmx){
  b1 <- bPair[1]
  b2 <- bPair[2]
  
  # estimate bias function for each time bin based on sampling distribution
  sampRows1 <- which(samp$b==b1)
  samp1 <- samp[sampRows1,e]
  # Reflecting the sample curve doesn't change the sp KDE much
  # except that the ends turn down a bit more (more convexity).
  # Since it's more complicated and throws warnings, don't do it.
  densSamp1 <- density(samp1, bw='SJ-ste', from=xmn, to=xmx)
  #  densSamp1 <- density.reflected(samp1, lower=xmn, upper=xmx) 
  w1 <- approxfun(densSamp1$x, densSamp1$y)
  
  # in the most recent time bin, there is no subsequent bin
  if (is.na(b2)){
    w2 <- NA
  } else {
    sampRows2 <- which(samp$b==b2)
    samp2 <- samp[sampRows2,e]
    densSamp2 <- density(samp2, bw='SJ-ste', from=xmn, to=xmx)
    #  densSamp2 <- density.reflected(samp2, lower=xmn, upper=xmx) 
    w2 <- approxfun(densSamp2$x, densSamp2$y)
  } 
  
  # redefine axis limits - discretization of samp KDE will shrink them a bit
  lim1 <- range(densSamp1$x)
  lim2 <- range(densSamp2$x)
  xmn <- max(lim1[1], lim2[1])
  xmx <- min(lim1[2], lim2[2])
  
  sList <- lapply(spp, function(s){
    nicher(dat = dat, b1 = b1, b2 = b2, s = s, env = e, 
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

# Note: might be weird/problematic to truncate based on surface temp,
# but then use in-habitat temp for niches. Temps at depth can get colder
# than the standardised lower bound from surface data.
# Perhaps an alternative to surface and in-habitat temp is to use MLD temp.
pt1 <- proc.time()
registerDoParallel(nCore)
kdeSum <- foreach(bPair=pairL) %dopar% 
   kde('temp_ym_0m', df, xmn, xmx)
stopImplicitCluster()
pt2 <- proc.time()
pt2-pt1

# remove NA rows (if a species is not sampled in the focal bin)
nas <- is.na(kdeSum$bin)
kdeSum <- kdeSum[!nas,]

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
  temp <- lapply(bins, sumup, s=s, dat=df, 
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
