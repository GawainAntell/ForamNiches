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

# The last step could introduce more species with <6 occs
# Subset again such that eas sp has >5 occs per bin
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

outNm <- paste0('Data/foram_uniq_occs_latlong_8ka_MeanAnnT_',day,'.csv')
write.csv(sppEnv, outNm, row.names = FALSE)
