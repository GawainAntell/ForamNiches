library(phytools)
library(paleoPhylo)
library(ape)
library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(PBSmapping)
library(ggplot2)

# set whether or not to truncate to standard global temperature range
doTrunc <- FALSE

# save names to put packages on all cores later
pkgs <- c('sp','raster') 

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

# Data import -------------------------------------------------------------

modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)

# Read in occurrence data
occ <- read.csv('Data/foram_uniq_occs_latlong_8ka_200108.csv', stringsAsFactors=FALSE)
bins <- unique(occ$bin)

# Limit analysis to species included in tree from Aze et al. 2011
trRaw <- read.csv('Data/Aze_et_al_2011_bifurcating_tree_data.csv', stringsAsFactors=FALSE)
trRaw$Species.name <- gsub('Globigerinoides sacculifer', 'Trilobatus sacculifer', trRaw$Species.name)
trRaw$Species.name <- gsub('Globigerinoides trilobus', 'Trilobatus trilobus', trRaw$Species.name)
tr_paleo <- with(trRaw, 
                 as.paleoPhylo(Species.code, Ancestor.code, Start.date, End.date)
)
tr <- buildApe(tr_paleo)
sppCodes <- tr$tip.label
rowOrdr <- match(sppCodes, trRaw$Species.code)
tr$tip.label <- trRaw$Species.name[rowOrdr]

# 9 species are not present in the phylogeny, mostly because microporiferate
sppAll <- unique(occ$species)
lostSpp <- setdiff(sppAll, tr$tip.label)

# Beella megastoma is arguably the same species as B. digitata,
# and Truncorotalia crassula is arguably senior synonym to crassaformis
# (Schiebel and Hemleben 2017). The depth ranges for both are unknown.
lostSpp <- c(lostSpp, 'Beella megastoma', 'Truncorotalia crassula')

spp <- setdiff(sppAll, lostSpp)
rows2toss <- ! occ$species %in% spp
occ <- occ[!rows2toss,]
row.names(occ) <- as.character(1:nrow(occ))

# Combine enviro and spp data ---------------------------------------------

envNm <- c('ann_temp_ym_dpth'
           #'month_temp_range',
           #'month_temp_max',
           #'month_temp_min',
           #'ann_otracer14_ym_dpth',
           #'ann_mixLyrDpth_ym_uo',
           #'ann_salinity_ym_dpth',
           #'ann_W_ym_dpth'
)
envNmShort <- sapply(envNm, function(txt){
  paste(strsplit(txt,'_')[[1]][2:3], collapse='_')
}) 

# Note: code below can deal with envNm that's a vector

llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Extract at each of 4 depths: 0, 40, 78, 164 m
# The lower 3 correspond to surface, surface-subsurface, and subsurface species.
dpths <- c(1,4,6,8)
# Compare with the Average Living Depth of species in the dataset:
sppDat <- read.csv('Data/foram_spp_data_200108.csv', stringsAsFactors=FALSE)
sameSpp <- sppDat$species %in% spp
sppDat <- sppDat[sameSpp,]
zones <- unique(sppDat$DepthHabitat)
for (z in zones){
 zBool <- sppDat$DepthHabitat==z
 avg <- mean(sppDat$ALD[zBool])
 avg <- round(avg)
 print(paste(z, avg))
}
# [1] "Subsurface 164"
# [1] "Surface.subsurface 93"
# [1] "Surface 49"

source('raster_brick_import_fcn.R')

addEnv <- function(bin, dat, binCol, cellCol, prj, envNm){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCells <- slc[,cellCol]
  slcEnv <- getBrik(bin=bin, envNm=envNm, mods=modId)
  
  for (i in 1:length(envNm)){
    env <- envNm[i]
    envVals <- raster::extract(slcEnv[[env]], slcCells)
    # Rows = points of extraction, columns = depth layers  
    envVals <- envVals[,dpths]
    nmOld <- envNmShort[i]
    nmNew <- paste(nmOld, c('0m','surf','surfsub','sub'), sep='_')
    slc[,nmNew] <- envVals
  } 
  slc
}

# Fast enough (1 min) this could be done in a loop/lapply rather than parallel
ncores <- detectCores() - 1
registerDoParallel(ncores)
sppEnv <- foreach(bin=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  addEnv(bin=bin,envNm=envNm,dat=occ,binCol='bin',cellCol='cell_number',prj=llPrj)
stopImplicitCluster()

# * Env at preferred depth ------------------------------------------------

habitatCol <- paste0(envNmShort, '_hab')
sppEnv[,habitatCol] <- NA
for (s in spp){
  sRow <- sppDat$species==s
  habitat <- sppDat$DepthHabitat[sRow]
  if (habitat=='Surface'){
    habNm <- paste0(envNmShort, '_surf')
  }
  if (habitat=='Surface.subsurface'){
    habNm <- paste0(envNmShort, '_surfsub')
  }
  if (habitat=='Subsurface'){
    habNm <- paste0(envNmShort, '_sub')
  }
  
  sBool <- sppEnv$species==s
  sppEnv[sBool, habitatCol] <- sppEnv[sBool, habNm]
}

# * Clean -----------------------------------------------------------------

# remove records where environment is unknown
envCol <- c(habitatCol, paste0(envNmShort, '_0m'))
if (length(envCol)==1){
  naRows <- is.na(sppEnv[,envCol])
} else {
  naRows <- apply(sppEnv[,envCol], 1, function(r)
    any(is.na(r))
  )
}
sppEnv <- sppEnv[!naRows,]

df <- sppEnv[,c('species','bin','cell_number','centroid_long','centroid_lat',envCol)]

# inspect correlation between temperature at surface (0m) and near preferred depth
cor(sppEnv$temp_ym_0m, sppEnv$temp_ym_hab)

# Truncate to standard temp range -----------------------------------------

if (doTrunc){
  
minmax <- function(df, b, env){
  bBool <- df[,'bin']==b
  slc <- df[bBool,]
  mn <- min(slc[,env])
  max <- max(slc[,env])
  c(b, mn, max)
}
sampSmryM <- sapply(bins, minmax, df=df, env='temp_ym_0m')
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
p
dev.off()

trunc <- data.frame()
outRows <- numeric()
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
  outRows <- c(outRows, out)
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
bars
dev.off()

# inspect the proportion of observations remaining
nrow(trunc)/nrow(df) # all data
table(old$trunc)['in range']/nrow(old) # excluding most recent 16 ka

# end case where data are truncated to standard temperature range
} else {
  trunc <- df
} 

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

# Also check for per-species continuity through time (at least 6 successive steps of 8ka)
# This precludes 7 species and ~5% of occurrences from analysis
keepSpp <- character()
binL <- bins[2] - bins[1]
enuf <- rep(binL, 6)
enufTxt <- paste0(enuf, collapse='')
spp <- unique(trunc$species)
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
any(diff(binsObs) != binL)

# Calculate 'sampled' environment from all occs in every bin
# analogous to the calculations for each species.
# Note: some species may be at same cell in a bin, so omit duplicates.
getSamp <- function(bin, dat, binCol, cellCol){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCells <- slc[,cellCol]
  dupes <- duplicated(slc[,cellCol])
  smpld <- slc[!dupes,]
  smpld$species <- 'sampled'
  rtrn <- rbind(slc, smpld)
  rtrn
}
outL <- lapply(bins, getSamp, dat=trunc, binCol='bin', cellCol='cell_number')
outDf <- do.call(rbind, outL)

if (doTrunc){
  outNm <- paste0('Data/foram_MAT_occs_latlong_8ka_trunc_',day,'.csv')
} else {
  outNm <- paste0('Data/foram_MAT_occs_latlong_8ka_',day,'.csv')
}
write.csv(outDf, outNm, row.names = FALSE)

# Simulate uniform niche spp ----------------------------------------------

df <- outDf 

# if starting from top of script, run the following lines to jump in from here:
  # df <- read.csv('Data/foram_MAT_occs_latlong_8ka_200127.csv',stringsAsFactors=FALSE)
  # bins <- unique(df$bin)
  # day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

spp <- unique(df$species)

sampBool <- df$species=='sampled'
samp <- df[sampBool,]

simNiche <- function(s, b, dat, sCol, bCol, samp){
  spBinBool <- dat[,sCol]==s & dat[,bCol]==b
  n <- sum(spBinBool)
  if (n>0){
    sampBinBool <- samp[,bCol]==b
    sampBin <- samp[sampBinBool,]
    rows <- 1:nrow(sampBin)
    rows4sim <- sample(rows, n)
    simDat <- sampBin[rows4sim,]
    simDat[,sCol] <- s
    return(simDat)
  }
}

simOverBins <- function(b){
  simL <- lapply(spp, b=b, simNiche, dat=df, sCol='species', bCol='bin', samp=samp)
  simDf <- do.call(rbind, simL)
}
simL <- lapply(bins, simOverBins)
simDf <- do.call(rbind, simL)

if (doTrunc){
  simNm <- paste0('Data/uniform_niche_sim_MAT_occs_latlong_8ka_trunc_',day,'.csv')
} else {
  simNm <- paste0('Data/uniform_niche_sim_MAT_occs_latlong_8ka_',day,'.csv')
}
write.csv(simDf, simNm, row.names=FALSE)

# Niche summary -----------------------------------------------------------

source('GSA_custom_ecospat_fcns.R')

nbins <- length(bins)
env <- 'temp_ym_hab'
h.method <- "nrd0" # "SJ-ste" # "ucv"
R <- 2^8

# * KDE niche summary -----------------------------------------------------

# Output niche overlap (though time), peak abundance, preferred enviro, & tolerance
nicher <- function(dat, b1, b2, sp, env, R, h.method){
  xmin <- min(dat[,env])
  xmax <- max(dat[,env])
  
  sp1rows <- which(dat$species==sp & dat$bin==b1)
  sp1 <- dat[sp1rows,env]
  
  sp2rows <- which(dat$species==sp & dat$bin==b2)
  sp2 <- dat[sp2rows,env]
  
  # for each species at time i and i+1
  z1 <- tryCatch(
    kdeNiche(sp=sp1, xmax=xmax, xmin=xmin, R=R, h.method=h.method),
    error = function(err){ list() }
  ) 
  z2 <- tryCatch(
    kdeNiche(sp=sp2, xmax=xmax, xmin=xmin, R=R, h.method=h.method),
    error = function(err){ list() }
  ) 
  
  # the species may be absent in one or both bins, in which case z is an empty list
  if (length(z1)==0){
    data.frame(bin=NA, sp=NA, n=NA, h=NA,
  #             errBaseLin=NA, errBaseSpl=NA, errSimpLin=NA, errSimpSpl=NA, errTrapLin=NA,
               pa=NA, pe=NA, tol=NA)
    
  } else {
    n <- length(sp1rows)
    if (length(z2)==0){
      data.frame(bin=b1, sp=sp, n=n, h=NA,
  #               errBaseLin=NA, errBaseSpl=NA, errSimpLin=NA, errSimpSpl=NA, errTrapLin=NA,
                 pa=z1$pa, pe=z1$pe, tol=z1$t)
    } else{
   #   err <- data.frame(sumErr(z1, z2))
      h <- hell(z1, z2) 
      data.frame(bin=b1, sp=sp, n=n, h=h,
             #    err,
                 pa=z1$pa, pe=z1$pe, tol=z1$t)
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

nichL <- lapply(spp, function(s){
  l <- apply(bPairs, 1, function(x){
    nicher(dat=df, b1=x[1], b2=x[2], sp=s, env=env, R=R, h.method=h.method)
  })
  do.call(rbind, l)
})
nich <- do.call(rbind, nichL)
# remove NA rows (if a species is not sampled in the focal bin)
nas <- is.na(nich$bin)
nich <- nich[!nas,]

simNichL <- lapply(spp, function(s){
  l <- apply(bPairs, 1, function(x){
    nicher(dat=simDf, b1=x[1], b2=x[2], sp=s, env=env, R=R, h.method=h.method)
  })
  do.call(rbind, l)
})
simNich <- do.call(rbind, simNichL)
simNas <- is.na(simNich$bin)
simNich <- simNich[!simNas,]

# * Non-KDE niche summary -------------------------------------------------

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
  temp <- lapply(bins, sumup, s=s, dat=df, binCol='bin', sCol='species', traitCol=envCol)
  tempDf <- do.call(rbind, temp)
})
rawSum <- do.call(rbind, rawSumL)

# combine KDE and raw-scale summary/sample statistics into one output file
fullSum <- merge(nich, rawSum)
if (doTrunc){
  sumNm <- paste0('Data/foram_niche_sumry_metrics_trunc_',day,'.csv') 
} else {
  sumNm <- paste0('Data/foram_niche_sumry_metrics_',day,'.csv')
}
write.csv(fullSum, sumNm, row.names=FALSE)

rawSimL <- lapply(spp, function(s){
  temp <- lapply(bins, sumup, s=s, dat=simDf, binCol='bin', sCol='species', traitCol=envCol)
  tempDf <- do.call(rbind, temp)
})
rawSim <- do.call(rbind, rawSimL)
fullSim <- merge(simNich, rawSim)
if (doTrunc){
  simNm <- paste0('Data/uniform_niche_sim_sumry_metrics_trunc_',day,'.csv')
} else {
  simNm <- paste0('Data/uniform_niche_sim_sumry_metrics_',day,'.csv')
}
write.csv(fullSim, simNm, row.names=FALSE)

# Inter-specific overlap --------------------------------------------------

interSppD <- function(b, df, env, R, h.method){
  xmax <- max(df[,env])
  xmin <- min(df[,env])
  
  bSppRows <- which(df$species!='sampled' & df$bin==b)
  bSpp <- unique(df$species[bSppRows])
  
  # Construct KDE of all species
  kdeL <- lapply(bSpp, function(s){
    spRows <- which(df$species==s & df$bin==b)
    sp1 <- df[spRows,env]
    z <- kdeNiche(sp1, xmax=xmax, xmin=xmin, R=R, h.method=h.method)
  }
  )
  names(kdeL) <- bSpp
  
  # Compute Hellinger's H for all species pairs
  fin <- data.frame()
  for (s1 in bSpp){
    for (s2 in bSpp){
      if (s1==s2) {next} else{
        h <- hell(kdeL[[s1]], kdeL[[s2]])
        pairDat <- data.frame(bin=b, sp1=s1, sp2=s2, h=h)
        fin <- rbind(fin, pairDat)
      }
    }
  }
  fin
}

registerDoParallel(ncores)
interSppDf <- foreach(bin=bins, .packages='pracma', .combine=rbind, .inorder=FALSE) %dopar%
  interSppD(b=bin, df=df, env=env, R=R, h.method=h.method)
stopImplicitCluster()
interSppNm <- paste0('Data/foram_species_pairs_KDE_H_', day, '.csv')
write.csv(interSppDf, interSppNm, row.names=FALSE)

registerDoParallel(ncores)
interSim <- foreach(bin=bins, .packages='pracma', .combine=rbind, .inorder=FALSE) %dopar%
  interSppD(b=bin, df=simDf, env=env, R=R, h.method=h.method)
stopImplicitCluster()
interSimNm <- paste0('Data/uniform_niche_sim_pairs_KDE_H_', day, '.csv')
write.csv(interSim, interSimNm, row.names=FALSE)
