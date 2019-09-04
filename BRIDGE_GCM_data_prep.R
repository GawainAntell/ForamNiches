library(sp)
library(raster)
library(ncdf4)
library(maptools)

# model age increment is 1ky for 0-24 ka, then 4 ky up to 800,000 ka
age <- c(0:23, seq(24, 800, by=4))

# model IDs start with 'tei', plus one of [n, N, o] in every combination w a suffix
dgt123 <- 'tei'
dgt4 <- c('n','N','o','O')
dgt5 <- c(letters, LETTERS, 0:9)
mods <- vector()
for (dgt4i in dgt4){
  for (dgt5i in dgt5){
    iMod <- paste0(dgt123, dgt4i, dgt5i)
    mods <- c(mods, iMod)
  }
}
mods <- head(mods, length(age))
modMeta <- data.frame(id=mods, age_1000ka=age)
# write.csv(modMeta, 'Data/gcm_model_codes.csv', row.names=FALSE)

# spp data binned at 4 ky resolution, so use same for GCMs
ageSteps <- seq(4, 800, by=4)
modSbset <- modMeta$age_1000ka %in% ageSteps
idSbset <- as.character(modMeta$id[modSbset])

# direct download each simulation .nc file from BRIDGE website
for (i in 1:length(idSbset)){
  age <- ageSteps[i]
  id <- idSbset[i]
  rt <- 'https://www.paleo.bristol.ac.uk/ummodel/data/'
  adrs <- paste0(rt, id, '/climate/', id, 'o.pgclann.nc') 
  dest <- paste0('Data/gcm_annual_mean/', id, 'o.pgclann_', age, '.nc')
  # Windows needs 'wb' argument to know that a binary transfer is necessary:
  temp <- download.file(adrs, dest, quiet=TRUE, mode='wb') 
}

# Extract data from each model and convert to raster brick 

# 'ym'= yearly mean
vars <- c(# 'W_ym_dpth', # vertical motion
        #  'ucurrTot_ym_dpth', # W-E current
        #  'vcurrTot_ym_dpth', # N-S current
        #  'salinity_ym_dpth',
          'temp_ym_dpth', # potential temperature
          'otracer14_ym_dpth', # age of ocean water
          'mixLyrDpth_ym_uo' # 1-D mixed layer depth
          )

# set up template raster, 1.25 degree increments (orignal scale)
rEmpt <- raster(ncols=288, nrows= 144, xmn=0, xmx=360, ymn=-90, ymx=90)
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
lat <- seq(-89.375, 89.375, length=144)
lon <- seq(0, 358.75, length=288)
llGrid <- expand.grid(lon, lat)
colnames(llGrid) <- c('lon','lat')

getRast <- function(iDp, vArry, coords, rEmpt, prj){ 
  vDp <- vArry[,,iDp]
  vDf <- data.frame('var'=as.vector(vDp))
  pts <- SpatialPointsDataFrame(coords=coords, data=vDf)
  r <- rasterize(pts, rEmpt, 'var', fun=mean)
  r <- rotate(r)
  projection(r) <- prj
  return(r)
}

# 25 minutes for 50 models, 2 variables
pt1 <- proc.time()
for (i in 1:length(idSbset)){
  age <- ageSteps[i]
  id <- idSbset[i]
  ann <- nc_open(paste0('Data/gcm_annual_mean/', id, 'o.pgclann_', age, '.nc'))

  dpths <- ann$dim$depth$vals
  nDpths <- length(dpths)

# Extract at successive depths (except for mixed layer depth, 1D) and export as brick
# Include age in file name or models that differ only in captalization will overwrite
for (v in vars) {
  vAnn <- ncvar_get(ann, varid=v)
  expNm <- paste0('Data/gcm_annual_mean/', id, age, '_ann_', v, '.tif')
  if (v=='mixLyrDpth_ym_uo'){
    r <- rEmpt
    values(r) <- t(vAnn)
    writeRaster(r, nl=1, filename=expNm, format="GTiff", bylayer=FALSE, overwrite=TRUE)
  } else {
    rastL <- lapply(1:nDpths, getRast, vArry=vAnn, coords=llGrid, rEmpt=rEmpt, prj=llPrj)
    vBrk <- brick(rastL)
    writeRaster(vBrk, nl=nDps, filename=expNm, format="GTiff", bylayer=FALSE, overwrite=TRUE)
  }
}

} # loop through models 
pt2 <- proc.time()
(pt2-pt1)/60


#################################################################
# Build rasters of most extreme monthly temperatures in each cell
# Note that only surface-level data are available

# Download monthly average .nc files
moVect <- c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

# 10 min for max values of 1 var in 50 models * 12 months
pt1 <- proc.time()
for (i in 1:length(idSbset)){
  age <- ageSteps[i]
  id <- idSbset[i]
  for (mo in moVect){
    rt <- 'https://www.paleo.bristol.ac.uk/ummodel/data/'
    adrs <- paste0(rt, id, '/climate/', id, 'o.pfcl', mo, '.nc') 
    dest <- paste0('Data/gcm_monthly_mean/', id, 'o.pfcl', mo, age, '.nc')
    # Windows needs 'wb' argument to know that a binary transfer is necessary:
    temp <- download.file(adrs, dest, quiet=TRUE, mode='wb') 
  }
  
  # Convert all months to (temporary) rasters
  # nc file to matrix to 1-column dataframe to spatialpoints to raster
  moFiles <- paste0('Data/gcm_monthly_mean/', id, 'o.pfcl', moVect, age, '.nc')
  moDat <- lapply(moFiles, nc_open)
  moM <- lapply(moDat, ncvar_get, varid='temp_mm_dpth')
  moDf <- lapply(moM, function(x) data.frame('var'=as.vector(x)))
  moPts <- lapply(moDf, SpatialPointsDataFrame, coords=llGrid)
  moR <- lapply(moPts, rasterize, y=rEmpt, field='var', fun=mean)
  moBrk <- brick(moR)

  # Find the hotest/coldest values per cell
  moMax <- max(moBrk)
  moMin <- min(moBrk)
  seasonal <- moMax - moMin
  seasonal <- rotate(seasonal)
  projection(seasonal) <- llPrj
  ssnlNm <- paste0('Data/gcm_monthly_mean/', id, age, '_month_temp_range.tif')
  writeRaster(seasonal, filename=ssnlNm, format='GTiff', overwrite=TRUE)
}
pt2 <- proc.time()
(pt2-pt1)/60
