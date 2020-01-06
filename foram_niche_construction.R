library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(PBSmapping)
library(ggplot2)

# Data import -------------------------------------------------------------

# save names to put packages on all cores later
pkgs <- c('sp','raster') 

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

# Read in occurrence data
source('read_foram_data.R')

modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
envNm <- c('ann_temp_ym_dpth'
            #'month_temp_range', 
            #'month_temp_max',
            #'month_temp_min',
            #'ann_otracer14_ym_dpth',
            #'ann_mixLyrDpth_ym_uo',
            #'ann_salinity_ym_dpth',
            #'ann_W_ym_dpth'
             )
# Note: envNm can be a vector
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Combine enviro and spp data ---------------------------------------------

getBrik <- function(bin, envNm){
  modRow <- modId$age_1000ka == bin
  id <- modId$id[modRow]
  
  # Load the rasters for only the desired env variables and time step
  allFls <- list.files('Data/', recursive = TRUE)
  txt <- paste0(id,'.*tif')
  modFls <- grep(txt, allFls)
  flNms <- paste0('Data/', allFls[modFls])
  envFlPos <- sapply(envNm, grep, flNms)
  envFlNms <- flNms[envFlPos]
  
  # Temperature raster files have 19 layers, 
  # but if using mix layer depth or BVF then modify code for 1 layer
  r <- lapply(envFlNms, brick)
  names(r) <- envNm
  r
}

addEnv <- function(bin, dat, binCol, cellCol, prj, envNm){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCells <- slc[,cellCol]
  slcEnv <- getBrik(bin=bin, envNm=envNm)
  
  for (env in envNm){
    envVals <- raster::extract(slcEnv[[env]], slcCells)
    # Rows = points of extraction, columns = depth layers  
    envVals <- envVals[,1]
    env <- paste0(env,'_surface')
    slc[,env] <- envVals
  } 
  
  slc
}

# Fast enough (1 min) this could be done in a loop/lapply rather than parallel
ncores <- detectCores() - 1
registerDoParallel(ncores)
sppEnv <- foreach(bin=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  addEnv(bin=bin,envNm=envNm,dat=occ,binCol='bin',cellCol='cell_number',prj=llPrj)
stopImplicitCluster()

# * Clean -----------------------------------------------------------------

# remove records where environment is unknown
envCol <- grep(envNm, colnames(sppEnv))
if (length(envCol)==1){
  naRows <- is.na(sppEnv[,envCol])
} else {
  test <- apply(sppEnv[,envCol], 1, function(r)
    any(is.na(r))
  )
}
sppEnv <- sppEnv[!naRows,]

df <- sppEnv[,c('species','bin','cell_number','centroid_long','centroid_lat','ann_temp_ym_dpth_surface')]
colnames(df)[ncol(df)] <- 'mat'

# Truncate to standard temp range -----------------------------------------

minmax <- function(df, b, env){
  bBool <- df[,'bin']==b
  slc <- df[bBool,]
  mn <- min(slc[,env])
  max <- max(slc[,env])
  c(b, mn, max)
}
sampSmryM <- sapply(bins, minmax, df=df, env='mat')
sampSmry <- data.frame(t(sampSmryM))
colnames(sampSmry) <- c('bin','min','max')
uppr <- min(sampSmry$max)
lwr <- max(sampSmry$min)

p <- ggplot(data=sampSmry) + theme_bw() +
  scale_x_continuous(name='Time (ka)', expand=c(0.01,0)) +
  scale_y_continuous(name = 'MAT (degrees C)') +
  geom_linerange(aes(x=-bin, ymin=min, ymax=max)) +
  geom_hline(yintercept=uppr, colour='blue', lwd=1) +
  geom_hline(yintercept=lwr, colour='blue', lwd=1)

pNm <- paste0('Figs/standardised_MAT_max_min_', day, '.pdf')
pdf(pNm, width = 6, height=4)
p
dev.off()

trunc <- data.frame()
for (b in bins){
  bBool <- df$bin==b
  slc <- df[bBool,]
  tooBig <- which(slc$mat > uppr)
  tooSmol <- which(slc$mat < lwr)
  out <- c(tooBig, tooSmol)
  slc <- slc[-out,]
  trunc <- rbind(trunc, slc)
}

# inspect the proportion of observations remaining
nrow(trunc)/nrow(df)

# * Clean -----------------------------------------------------------------

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

# note: no remaining observations at 740 ka!
binsObs <- sort(unique(trunc$bin))
any(diff(binsObs) != binL)
trim <- trunc$bin < 740
trimmd <- trunc[trim,]
binsNew <- sort(unique(trimmd$bin))

# Calculate 'sampled' environment from all occs in every bin
# analogous to the calculations for each species
# Note: some species may be at same cell in a bin, so omit duplicates
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
outL <- lapply(binsNew, getSamp, dat=trimmd, binCol='bin', cellCol='cell_number')
outDf <- do.call(rbind, outL)

outNm <- paste0('Data/foram_MAT_occs_latlong_8ka_',day,'.csv')
write.csv(outDf, outNm, row.names = FALSE)

# KDE niche summary -------------------------------------------------------

# * Data prep -------------------------------------------------------------

source('ecospat.grid.clim.dyn.GSA.fcn.R')

df <- outDf # if starting from top of script
# df <- read.csv('Data/foram_MAT_occs_latlong_8ka_200106.csv',stringsAsFactors=FALSE)
bins <- unique(df$bin)
nbins <- length(bins)
spp <- unique(df$species)

env <- 'mat'
h.method <- "nrd0" # "SJ-ste" # "ucv"
# Resolution of the gridding of the climate space. Ecospat default is 100.
R <- 1000

# Calculate niche overlap (Schoener's D), peak abundance, preferred enviro, & tolerance
nicher <- function(b1, b2, sp, env, h.method){
  
  globBool <- df$species=='sampled'
  glob <- df[globBool,env]
  
  glob1rows <- which(df$species=='sampled' & df$bin==b1)
  glob1 <- df[glob1rows,env]
  
  glob2rows <- which(df$species=='sampled' & df$bin==b2)
  glob2 <- df[glob2rows,env]
  
  sp1rows <- which(df$species==sp & df$bin==b1)
  sp1 <- df[sp1rows,env]
  
  sp2rows <- which(df$species==sp & df$bin==b2)
  sp2 <- df[sp2rows,env]
  
  # for each species at time i and i+1
  z1 <- tryCatch(
    GSA.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0, th.env=0, h.method=h.method),
    error = function(err){ list() }
  ) 
  z2 <- tryCatch(
    GSA.grid.clim.dyn(glob, glob2, sp2, R, th.sp=0, th.env=0, h.method=h.method),
    error = function(err){ list() }
  ) 
  
  # the species may be absent in one or both bins, in which z1|2 is an empty vector
  if (length(z1)==0){
    data.frame(bin=NA, sp=NA, n=NA, d=NA, pa=NA, pe=NA, tol=NA)
  } else {
    n <- length(sp1rows)
    if (length(z2)==0){
      data.frame(bin=b1, sp=sp, n=n, d=NA, pa=z1$pa, pe=z1$pe, tol=z1$t)
    } else{
      ovrlp <- GSA.ecospat.niche.overlap(z1, z2, cor=FALSE)
      data.frame(bin=b1, sp=sp, n=n, d=ovrlp$D, pa=z1$pa, pe=z1$pe, tol=z1$t)
    }
  }
}

# the older bin is column 1, the younger is column 2
bPairs <- cbind(bins[-1], bins[-length(bins)])
nichL <- lapply(spp, function(s){
  l <- apply(bPairs, 1, function(x){
    nicher(b1=x[1], b2=x[2], sp=s, env=env, h.method=h.method)
  })
  do.call(rbind, l)
})
nich <- do.call(rbind, nichL)
# remove NA rows
nas <- is.na(nich$bin)
nich <- nich[!nas,]

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')
dfNm <- paste0('Data/foram_niche_sumry_metrics_KDE_',day,'.csv')
write.csv(nich, dfNm, row.names=FALSE)

# Raw value niche summary -------------------------------------------------
# Note that using the KDE approach allows calculation of Schoener's D overlap.
# Can't calculate overlap from the last time step because no modern data included.
# So omit data from the most recent time step here, for comparability.

recent <- df$bin %in% min(bins)
df <- df[!recent,]

# Need mean, variance, sample size, and age of trait values (MAT) for each sp & bin
sumup <- function(bin, s, dat, binCol, sCol, traitCol){
  slcRows <- which(dat[,binCol] == bin & dat[,sCol] == s)
  if (length(slcRows)>0){
    slc <- dat[slcRows,]
    x <- slc[,traitCol]
    m <- mean(x)
    v <- sd(x)^2
    n <- length(x)
    rtrn <- data.frame(bin=bin, sp=s, m=m, v=v, n=n)
  } else {
    rtrn <- data.frame()
  }
  return(rtrn)
}
rawSumL <- lapply(spp, function(s){
  temp <- lapply(bins[-1], sumup, s=s, dat=df, binCol='bin', sCol='species', traitCol='mat')
  tempDf <- do.call(rbind, temp)
})
rawSum <- do.call(rbind, rawSumL)

rawSumNm <- paste0('Data/foram_niche_sumry_metrics_raw_values_',day,'.csv')
write.csv(rawSum, rawSumNm, row.names=FALSE)
