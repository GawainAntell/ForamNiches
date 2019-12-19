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
envNms <- c('ann_temp_ym_dpth', 
            'month_temp_range', 
            'month_temp_max',
            'month_temp_min'
            #'ann_otracer14_ym_dpth',
            # 'ann_mixLyrDpth_ym_uo',
            # 'ann_salinity_ym_dpth',
            #'ann_W_ym_dpth'
             )
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

getBrik <- function(bin, envNms){
  modRow <- modId$age_1000ka == bin
  id <- modId$id[modRow]
  
  # Load the rasters for only the desired env variables and time step
  allFls <- list.files('Data/', recursive = TRUE)
  txt <- paste0(id,'.*tif')
  modFls <- grep(txt, allFls)
  flNms <- paste0('Data/', allFls[modFls])
  envFlPos <- sapply(envNms, grep, flNms)
  envFlNms <- flNms[envFlPos]
  
  # Temperature raster files have 19 layers, 
  # but if using mix layer depth or BVF then modify code for 1 layer
  r <- lapply(envFlNms, brick)
  names(r) <- envNms
  r
}

addEnv <- function(bin, dat, binCol, cellCol, prj, envNms){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCells <- slc[,cellCol]
#  slcCoords <- slc[,coordCols]
#  slcPts <- SpatialPoints(slcCoords, proj4string = CRS(prj))
  
  slcEnv <- getBrik(bin=bin, envNms=envNms)
  
  for (env in envNms){
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
  addEnv(bin=bin,envNms=envNms,dat=occ,binCol='bin',cellCol='cell_number',prj=llPrj)
stopImplicitCluster()

# Check for collinearity
# Beware if using water age: at surface all values = 0
smpld <- sppEnv[sppEnv$species=='sampled',]
envCols <- paste0(envNms,'_surface')
smpldEnv <- smpld[,envCols]
smpldEnv <- na.omit(smpldEnv)
cor(smpldEnv)
# mean annual temp and min and max monthly temp are nearly perfectly colinear

# Retrieve eigenvalues for PC rotation of sampled environment.
# Axis 1 ~ mean/max temperature. Axis 2 ~ seasonal temperature range
pc <- prcomp(smpldEnv, center=TRUE, scale=TRUE) 
pcRot <- pc$rotation
rotNm <- paste0('Data/all_sampled_env_PC_rotation_',day,'.csv')
write.csv(pcRot, rotNm)

# How much variance is explained by axies 1?
pcSmry <- summary(pc)
pcSmry$importance['Proportion of Variance','PC1']
#biplot(pc)

pcAll <- stats::predict(pc, sppEnv[envCols])
pcAxes <- paste0('pc', 1:length(envCols))
sppEnv[,pcAxes] <- pcAll
naRows <- is.na(pcAll[,1])
sppEnv <- sppEnv[!naRows,]

# again remove species for which <6 occs are found
# Subset occurrences such that eas sp has >5 occs per bin
tooRare <- function(sp, bin, df){
  spRows <- which(df$species==sp & df$bin==bin)
  if (length(spRows)<6){
    spRows
  }
}  
tossRowsL <- 
  sapply(spp, function(x){
    sapply(bins, function(b){
      tooRare(sp=x, bin=b, df=sppEnv)
    } )
  } )
tossRows <- unlist(tossRowsL)
sppEnv <- sppEnv[-tossRows,]

outNm <- paste0('Data/foram_uniq_occs_latlong_8ka_wEnv_',day,'.csv')
write.csv(sppEnv, outNm, row.names = FALSE)
