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

# spp data binned at 16 ky resolution, so use same for GCMs
ageSteps <- seq(8, 782, by=16)
modSbset <- modMeta$age_1000ka %in% ageSteps
idSbset <- as.character(modMeta$id[modSbset])

# direct download each simulation .nc file from BRIDGE website
for (id in idSbset){
  rt <- 'https://www.paleo.bristol.ac.uk/ummodel/data/'
  adrs <- paste0(rt, id, '/climate/', id, 'o.pgclann.nc') 
  dest <- paste0('Data/gcm_annual_mean/', id, 'o.pgclann.nc')
  # Windows needs 'wb' argument to know that a binary transfer is necessary:
  temp <- download.file(adrs, dest, quiet=TRUE, mode='wb') 
}

# Extract data from each model and convert to raster brick 

# 'ym'= yearly mean
vars <- c(# 'W_ym_dpth', # vertical motion
        #  'ucurrTot_ym_dpth', # W-E current
        #  'vcurrTot_ym_dpth', # N-S current
          'temp_ym_dpth', # potential temperature
          'salinity_ym_dpth',
          'otracer14_ym_dpth' # age of ocean water
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

# 38 minutes for 49 models, 3 variables
pt1 <- proc.time()
for (id in idSbset){
ann <- nc_open(paste0('Data/gcm_annual_mean/', id, 'o.pgclann.nc'))

dpths <- ann$dim$depth$vals
nDpths <- length(dpths)

# extract at successive depths and export as brick
for (v in vars) {
  vAnn <- ncvar_get(ann, varid=v)
  rastL <- lapply(1:nDpths, getRast, vArry=vAnn, coords=llGrid, rEmpt=rEmpt, prj=llPrj)
  vBrk <- brick(rastL) 
  expNm <- paste0('Data/gcm_annual_mean/', id, '_ann_', v, '.tif')
  writeRaster(vBrk, nl=nDps, filename=expNm, format="GTiff", bylayer=FALSE, overwrite=TRUE)
}

} # loop through models 
pt2 <- proc.time()
pt2-pt1

# NOTES for next iteration of analysis:
# Average models within a time interval, e.g. 0ka, 4ka, 8ka, and 12ka models 
# averaged to match with 0-16ka foram data.
# Transform lat-long grid to equal-area projection.

#################################################################
# Build rasters of most extreme monthly temperatures in each cell
# Note that only surface-level data are available

# Download monthly average .nc files
moVect <- c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

# 10 min for max values of 1 var in 49 models * 12 months
pt1 <- proc.time()
for (id in idSbset){
  for (mo in moVect){
  rt <- 'https://www.paleo.bristol.ac.uk/ummodel/data/'
  adrs <- paste0(rt, id, '/climate/', id, 'o.pfcl', mo, '.nc') 
  dest <- paste0('Data/gcm_monthly_mean/', id, 'o.pfcl', mo, '.nc')
  # Windows needs 'wb' argument to know that a binary transfer is necessary:
  temp <- download.file(adrs, dest, quiet=TRUE, mode='wb') 
  }
  
  # Convert all months to (temporary) rasters
  # nc file to matrix to 1-column dataframe to spatialpoints to raster
  moFiles <- paste0('Data/gcm_monthly_mean/', id, 'o.pfcl', moVect, '.nc')
  moDat <- lapply(moFiles, nc_open)
  moM <- lapply(moDat, ncvar_get, varid='temp_mm_dpth')
  moDf <- lapply(moM, function(x) data.frame('var'=as.vector(x)))
  moPts <- lapply(moDf, SpatialPointsDataFrame, coords=llGrid)
  moR <- lapply(moPts, rasterize, y=rEmpt, field='var', fun=mean)
  moBrk <- brick(moR)

  # Find the hotest/coldest values per cell
  moMax <- max(moBrk)
  moMax <- rotate(moMax)
  projection(moMax) <- llPrj
  maxNm <- paste0('Data/gcm_monthly_mean/', id, '_max_month_temp.tif')
  writeRaster(moMax, filename=maxNm, format='GTiff', overwrite=TRUE)
}
pt2 <- proc.time()
(pt2-pt1)/60
