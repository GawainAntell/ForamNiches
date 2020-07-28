library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

# Data prep ---------------------------------------------------------------

source('species_kde_buildr.R')
day <- as.Date(date(), format = "%a %b %d %H:%M:%S %Y")
spAttr <- read.csv('Data/foram-spp-data_2020-07-21.csv',
                   stringsAsFactors = FALSE)
dList <- readRDS('Data/spp-and-sampling-data_list-by-depth_2020-07-21.rds')
envNm <- 'temp_ym'
spp <- unique(dList$temp_ym_0m$sp$species)
bins <- unique(dList$temp_ym_0m$sp$bin)
nCore <- detectCores() - 1

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

# * WARNING * - this could take a couple hours
bw <- 'nrd0'
pkg <- c('pracma','GoFKernel','kerneval')
pt1 <- proc.time()
registerDoParallel(nCore)
kdeSum <- foreach(dat=dList[2:4], .combine=rbind, .inorder=FALSE, .packages=pkg) %:%
  foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar%
   kde(dat, bPair, envNm, bw = bw)
kdeSumSS <- foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar%
  kde(dList[[1]], bPair, envNm, bw = bw)
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
envL <- lapply(dList[2:4], function(dat){
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
  binL <- lapply(bins, sumup, s = s, dat = dList[[1]], 
                 binCol = 'bin', sCol = 'species', traitCol = envNm)
  binDf <- do.call(rbind, binL)
})
rawSumSS <- do.call(rbind, envLss)

# combine KDE and finite sample statistics into same output
fullSum <- merge(kdeSum, rawSum, all.x = TRUE)
sumNm <-   paste0('Data/niche-sumry-metrics_', bw, '_hab_', day, '.csv')
write.csv(fullSum, sumNm, row.names=FALSE)
fullSumSS <- merge(kdeSumSS, rawSumSS, all.x = TRUE)
sumSSnm <- paste0('Data/niche-sumry-metrics_', bw, '_SS_', day, '.csv')
write.csv(fullSumSS, sumSSnm, row.names = FALSE)
