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
  
  # Calculate 'sampled' environment from all occs in every bin
  # analogous to the calculations for each species
  # Note: some species may be at same cell in a bin, so omit duplicates
  dupes <- duplicated(slc[,cellCol])
  smpld <- slc[!dupes,]
  smpld$species <- 'sampled'
  slc <- rbind(slc, smpld)
  
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
diff(sort(bins))
trim <- trunc$bin < 740
trimmd <- trunc[trim,]
outNm <- paste0('Data/foram_MAT_occs_latlong_8ka_',day,'.csv')
write.csv(trimmd, outNm, row.names = FALSE)
