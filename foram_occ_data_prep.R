library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

occ <- read.csv('Data/Foram_data_800ka_IF190709.csv', stringsAsFactors = FALSE)

# easier to work on scale of data: ka not ma
occ$age <- occ$age*1000

# restrict to range-through species
tooYng <- which(occ$earliest < .8)
occ <- occ[-tooYng,]
# species with 'latest' values older than occurrences are suspect
tooOld <- which(occ$latest > .016) 
occ <- occ[-tooOld,]

# excise sp traits into a separate, sp-level file to call later in analysis
traits <- c('morphDes','spinose','structure','aperture','latest','earliest',
            'mph','mphRef','eco','ecoRef','geo','geoRef',
            'Average.Area.microm2','log.area','Approximate.diameter')
sppDat <- occ[,c('species',traits)]
dupes <- duplicated(sppDat$species)
sppDat <- sppDat[!dupes,]
# write.csv(sppDat, 'Data/foram_spp_data.csv', row.names = FALSE)

# cut down file size
irrel <- c('site','hole','core','section','sample.top','sample.type',
           'preservation','processing','Macro.micro')
cols2toss <- colnames(occ) %in% c(traits, irrel)
occ <- occ[,!cols2toss]

# bin to 16 ky intervals to match GCM data

brk <- seq(0, 784, by=16)
bins <- data.frame(t=brk, b=brk+16, mid=brk+8)

# age 'zero' = 1950, and some observations are more recent (i.e. negative age)
bins$t[1] <- -0.1

for (r in 1:nrow(occ)){
  a <- occ$age[r]
  bin <- which(bins$t < a & bins$b >= a)
  aBin <- bins$mid[bin]
  occ$bin[r] <- aBin
}

# Rasterize to the resolution of GCM data

rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
  #eck.proj <- "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
  #rEmpt <- projectRaster(rEmpt, crs = eck.proj)
  #rEmpt <- setValues(rEmpt, c(1:ncell(rEmpt))) 
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
palCoords <- occ[,c('pal.long','pal.lat')]
pts <- SpatialPointsDataFrame(palCoords, data = occ, proj4string = CRS(llPrj))
  #pts_eck6 <- spTransform(pts, eck.proj)

occ$cell_number <- cellFromXY(rEmpt, pts) 
occ$centroid_lat <- occ$centroid_long <- NA
occ[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, occ$cell_number)

# remove records with imprecise age estimates
unconstrnd <- which(occ$age.model %in% c('Berggren1977','Ericson1968','GTSBlow1969'))
occ <- occ[-unconstrnd,]

# remove occurences that are singular in space-time
nmFreq <- table(occ$species)
singles <- names( which(nmFreq == 1) )
occ <- occ[! occ$species %in% singles,]

spp <- unique(occ$species)

# subset by stage, then species, then find unique species-cells combinations

# omit duplicate cell-species combinations
st_sp_summary <- function(dat, sp){
  sp_st <- dat[dat$species==sp,] 
  uniq_record <- function(x){ 
    first <- which(sp_st$cell_number==x)[1]
    sp_st[first,] 
  }
  temp_list <- lapply(unique(sp_st$cell_number), uniq_record)
  do.call('rbind', temp_list)
}

# apply st_sp_summary over species in a stage
st_summary <- function(dat, bin){
  slc <- dat[dat$bin==bin,]
  temp_list <- lapply(unique(slc$species), st_sp_summary, dat=slc)
  do.call('rbind', temp_list)
}

nCore <- detectCores() - 1
pt1 <- proc.time()
registerDoParallel(nCore)
stage_dat_list <- foreach(bin=bins$mid) %dopar% st_summary(bin=bin, dat=occ)
stopImplicitCluster()
pt2 <- proc.time()
pt2-pt1

fin <- do.call('rbind', stage_dat_list)


day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')
df_name <- paste('Data/foram_uniq_occs_latlong', day, '.csv', sep='')
write.csv(fin, file=df_name, row.names=FALSE)

