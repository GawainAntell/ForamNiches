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
library(adehabitatMA)
library(adehabitatHR)
library(ecospat)
library(ggplot2)

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
# Note: code below can deal with envNm that's a vector

llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

source('raster_brick_import_fcn.R')

addEnv <- function(bin, dat, binCol, cellCol, prj, envNm){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCells <- slc[,cellCol]
  slcEnv <- getBrik(bin=bin, envNm=envNm, mods=modId)
  
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
  tooBig <- which(slc$mat > uppr)
  tooSmol <- which(slc$mat < lwr)
  out <- c(tooBig, tooSmol)
  if (length(out) > 0){
    slc <- slc[-out,]
  }
  trunc <- rbind(trunc, slc)
  outRows <- c(outRows, out)
}

# * Evaluate degree of truncation -----------------------------------------

df$trunc <- 'in range'
tooBig <- which(df$mat > uppr)
tooSmol <- which(df$mat < lwr)
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

outNm <- paste0('Data/foram_MAT_occs_latlong_8ka_',day,'.csv')
write.csv(outDf, outNm, row.names = FALSE)

# KDE niche summary -------------------------------------------------------

# * Data prep -------------------------------------------------------------

source('GSA_custom_ecospat_fcns.R')

df <- outDf 

# if starting from top of script, run the following lines to jump in from here:
# df <- read.csv('Data/foram_MAT_occs_latlong_8ka_200116.csv',stringsAsFactors=FALSE)
# bins <- unique(df$bin)

nbins <- length(bins)
spp <- unique(df$species)
env <- 'mat'
h.method <- "nrd0" # "SJ-ste" # "ucv"
# Resolution of the gridding of the climate space. Ecospat default is 100.
# Note that when the value is higher, the density sum will be increasingly
# larger than unity. Possible reason (?): R interpolates between discrete points
# using a trapezoid. This overestimates the integral in concave-up curves,
# which are the tails of the niche distributions.
# (https://commons.wikimedia.org/w/index.php?curid=8541370)
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
  
  # the species may be absent in one or both bins, in which case z is an empty list
  if (length(z1)==0){
    data.frame(bin=NA, sp=NA, n=NA, d=NA, pa=NA, pe=NA, tol=NA)
  } else {
    n <- length(sp1rows)
    if (length(z2)==0){
      data.frame(bin=b1, sp=sp, n=n, d=NA, pa=z1$pa, pe=z1$pe, tol=z1$t)
    } else{
      ovrlp <- GSA.ecospat.niche.overlap(z1, z2, cor=FALSE)
      data.frame(bin=b1, sp=sp, n=n, d=ovrlp, pa=z1$pa, pe=z1$pe, tol=z1$t)
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
# remove NA rows (if a species is not sampled in the focal bin)
nas <- is.na(nich$bin)
nich <- nich[!nas,]

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')
dfNm <- paste0('Data/foram_niche_sumry_metrics_KDE_',day,'.csv')
write.csv(nich, dfNm, row.names=FALSE)

# * Example KDE plots -----------------------------------------------------

s <- 'Neogloboquadrina pachyderma'
b1 <- 204
b2 <- 404

globBool <- df$species=='sampled'
glob <- df[globBool,env]
glob <- as.matrix(glob)

xmax <- max(glob[, 1])
xmin <- min(glob[, 1])
x <- seq(from = min(glob[, 1]), to = max(glob[, 1]), 
         length.out = R)

sp1rows <- which(df$species==s & df$bin==b1)
sp1 <- df[sp1rows,env]
sp1 <- as.matrix(sp1)

sp2rows <- which(df$species==s & df$bin==b2)
sp2 <- df[sp2rows,env]
sp2 <- as.matrix(sp2)

sp1dens <- density(sp1[, 1], 
                   kernel = "gaussian", 
                   bw=h.method,
                   n = R, 
                   cut = 3
)
sp2dens <- density(sp2[, 1], 
                   kernel = "gaussian", 
                   bw=h.method,
                   n = R, 
                   cut = 3
)

p1nm <- paste0('Figs/KDE_Npachyderma_', b1, 'ka_', day, '.pdf')
pdf(p1nm, width=5, height=5)
  plot(sp1dens, xlim=c(xmin-6, xmax+6), main=paste(b1, 'ka,', s))
  rug(sp1[,1])
  polygon(sp1dens$x, sp1dens$y, col='orange')
dev.off()

p2nm <- paste0('Figs/KDE_Npachyderma_', b2, 'ka_', day, '.pdf')
pdf(p2nm, width=5, height=5)
  plot(sp2dens, xlim=c(xmin-6, xmax+6), main=paste(b2, 'ka,', s))
  rug(sp2[,1])
  polygon(sp2dens$x, sp2dens$y, col='orange')
dev.off()

p1 <- as.matrix(sp1dens$y)/sum(as.matrix(sp1dens$y))
p2 <- as.matrix(sp2dens$y)/sum(as.matrix(sp2dens$y))
absDiff <- abs(p1 - p2)
D <- 1 - (0.5 * (sum(absDiff)))
D


# * Inter-specific overlap ------------------------------------------------

interSppD <- function(b, df, R, h.method){
  globBool <- df$species=='sampled'
  glob <- df[globBool,env]
  
  glob1rows <- which(df$species=='sampled' & df$bin==b)
  glob1 <- df[glob1rows,env]
  
  bSppRows <- which(df$species!='sampled' & df$bin==b)
  bSpp <- unique(df$species[bSppRows])
  
  # Construct KDE of all species
  kdeL <- lapply(bSpp, function(s){
    spRows <- which(df$species==s & df$bin==b)
    sp1 <- df[spRows,env]
    z <- GSA.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0, th.env=0, h.method=h.method)
  }
  )
  names(kdeL) <- bSpp
  
  # Compute Schoener's D for all species pairs
  fin <- data.frame()
  for (s1 in bSpp){
    for (s2 in bSpp){
      if (s1==s2) {next} else{
        ovrlpL <- GSA.ecospat.niche.overlap(kdeL[[s1]], kdeL[[s2]], cor=FALSE)
        d <- ovrlpL$D
        pairDat <- data.frame(bin=b, sp1=s1, sp2=s2, d=d)
        fin <- rbind(fin, pairDat)
      }
    }
  }
  fin
}
interSppL <- lapply(bins, interSppD, df=df, R=R, h.method=h.method)
interSppDf <- do.call(rbind, interSppL)
interSppNm <- paste0('Data/foram_species_pairs_KDE_D_', day, '.csv')
write.csv(interSppDf, interSppNm, row.names=FALSE)

# Raw value niche summary -------------------------------------------------

# Need mean, variance, sample size, and age of trait values (MAT) for each sp & bin
sumup <- function(bin, s, dat, binCol, sCol, traitCol){
  slcRows <- which(dat[,binCol] == bin & dat[,sCol] == s)
  if (length(slcRows)>0){
    slc <- dat[slcRows,]
    x <- slc[,traitCol]
    m <- mean(x)
    sd <- sd(x)
    n <- length(x)
    rtrn <- data.frame(bin=bin, sp=s, m=m, sd=sd, n=n)
  } else {
    rtrn <- data.frame()
  }
  return(rtrn)
}
rawSumL <- lapply(spp, function(s){
  temp <- lapply(bins, sumup, s=s, dat=df, binCol='bin', sCol='species', traitCol='mat')
  tempDf <- do.call(rbind, temp)
})
rawSum <- do.call(rbind, rawSumL)

rawSumNm <- paste0('Data/foram_niche_sumry_metrics_raw_values_',day,'.csv')
write.csv(rawSum, rawSumNm, row.names=FALSE)
