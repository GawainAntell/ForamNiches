library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(PBSmapping)

# save names to put packages on all cores later
pkgs <- c('sp','raster') 

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

# Read in occurrence data
source('read_foram_data.R')

modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
envNms <- c('ann_otracer14_ym_dpth','ann_temp_ym_dpth',
             'ann_mixLyrDpth_ym_uo','month_temp_range'
            # 'ann_salinity_ym_dpth','ann_W_ym_dpth','max_month_temp'
             )
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

getR <- function(bin, envNms){
  modRow <- modId$age_1000ka == bin
  id <- modId$id[modRow]
  
  # Annual temp and water age files have 19 layers, 
  # but temp seasonality and mix layer depth have only 1
  rAll <- list.files('Data/', recursive = TRUE)
  txt <- paste0(id,'.*tif')
  rMod <- grep(txt, rAll)
  thickPos <- intersect(rMod, grep('otracer|ann_temp', rAll))
  flatPos <- intersect(rMod, grep('mixLyrDpth|temp_range', rAll))
  
  brickFls <- paste0('Data/', rAll[thickPos])
  flatFls <- paste0('Data/', rAll[flatPos])
  r1 <- lapply(brickFls, brick)
  r2 <- lapply(flatFls, raster)
  r <- append(r1, r2)
  names(r) <- envNms
  r
}

addEnv <- function(bin, dat, binCol, coordCols, cellCol, prj, envNms){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCoords <- slc[,coordCols]
  slcPts <- SpatialPoints(slcCoords, proj4string = CRS(prj))
  
  slcR <- getR(bin=bin, envNms=envNms)
  
  for (env in envNms){
    envVals <- raster::extract(slcR[[env]], slcPts)
    
    # If extracted from Brick, values are put in a matrix
    # Rows = points of extraction, columns = depth layers  
    if (length(envVals)!=length(slcPts)){
      envVals <- envVals[,1]
      env <- paste0(env,'_surface')
    }
    slc[,env] <- envVals
  } 
  
  # Calculate 'sampled' environment from all occs in every bin
  # Calculated the same as for each species, so treat as if a species
  # Note: some species may be at same cell in a bin
  dupes <- duplicated(slc[,cellCol])
  smpld <- slc[!dupes,]
  smpld$species <- 'sampled'
  slc <- rbind(slc, smpld)
  
  slc
}

# Fast enough (1 min) this could be done in a loop/lapply rather than parallel
ncores <- detectCores() - 1
pt1 <- proc.time()
registerDoParallel(ncores)

sppEnv <- foreach(bin=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  addEnv(bin=bin, envNms=envNms, dat=occ, binCol='bin', cellCol='cell_number',
           coordCols=c('pal.long','pal.lat'), prj=llPrj
           )

stopImplicitCluster()
pt2 <- proc.time()
(pt2-pt1)/60

flNm <- paste0('Data/foram_uniq_occs_latlong_8ka_wEnv_',day,'.csv')
# write.csv(sppEnv, flNm)

# Check for collinearity - all good
# But if values taken at surface, water age is all 0
smpld <- sppEnv[sppEnv$species=='sampled',]
envCols <- c('ann_temp_ym_dpth_surface','ann_mixLyrDpth_ym_uo','month_temp_range')
  # ncol(smpld) - (length(envNms)-1):0
smpldEnv <- smpld[,envCols]
smpldEnv <- na.omit(smpldEnv)
cor(smpldEnv)

# Retrieve eigenvalues for PC rotation of sampled environment.
pc <- prcomp(smpldEnv, center=TRUE, scale=TRUE) 
pcRot <- pc$rotation
rotNm <- paste0('Data/all_sampled_env_PC_rotation_',day,'.csv')
write.csv(pcRot, rotNm)

# Print importance of each component
# Axis 1 is half of variance
summary(pc)

pcAll <- stats::predict(pc, sppEnv[envCols])
pcAxes <- paste0('pc', 1:length(envCols))
sppEnv[,pcAxes] <- pcAll
naRows <- is.na(pcAll[,1])
sppEnv <- sppEnv[!naRows,]

# Calculate center of mass and 90% hypervolume of each species' hypervolume
# But only if > 5 occurrence for a given time bin
getNiche <- function(sp, bin){
  spRows <- which(sppEnv$species==sp & sppEnv$bin==bin)
  if (length(spRows)<6){
    return(rep(NA,7))
  }
  spDat <- sppEnv[spRows,]
  n <- length(spRows)
  cm <- colMeans(spDat[,pcAxes])
  xy <- data.frame(X=spDat$pc1, Y=spDat$pc2)
  chull <- calcConvexHull(xy)
  cha <- calcArea(chull)$area
  c(sp,bin,n,cha,cm)
}

nich <- matrix(nrow=0,ncol=7)
for (s in spp){
  for (b in bins){
    spBinNich <- getNiche(sp=s, bin=b)
    nich <- rbind(nich, spBinNich)
  }
}
nich <- na.omit(nich)
colnames(nich) <- c('species','bin','n_sites','convex_hull_2axes',
                    'pc1_centroid','pc2_centroid','pc3_centroid')
rownames(nich) <- 1:nrow(nich)
nichFlNm <- paste0('Data/foram_niche_data_8ky_',day,'.csv')
write.csv(nich, nichFlNm)
