library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

# Read in occurrence data
occ <- read.csv('Data/Foram_data_800ka_IF190709.csv', stringsAsFactors = FALSE)

# easier to work on scale of data: ka not ma
occ$age <- occ$age*1000

# remove T. cavernula and excelsa, which originated <800 ka
tooYng <- which(occ$earliest < .8)
occ <- occ[-tooYng,]

# excise sp traits into a separate, sp-level file to call later in analysis
traits <- c('morphDes','spinose','structure','aperture','latest','earliest',
            'mph','mphRef','eco','ecoRef','geo','geoRef',
            'Average.Area.microm2','log.area','Approximate.diameter')
sppDat <- occ[,c('species',traits)]
dupes <- duplicated(sppDat$species)
sppDat <- sppDat[!dupes,]

##########################################
# add modern depth ranges to sp-level file, for species in both datasets

# Read in modern depth ranges: Rebotim et al 2017, table 4
dpthTbl <- read.table('Data/Rebotim_et_al_depth_ranges.txt',
                      header=TRUE, stringsAsFactors=FALSE)
abund <- dpthTbl$MaxAbund
penult <- nchar(abund) - 1
dpthTbl$MaxAbund <- substr(abund, 1, penult)
dpthTbl$MaxAbund <- as.numeric(dpthTbl$MaxAbund)

# Add G. trilobus data from Loncaric et al. 2006
trilobus <- c('Globigerinoides.trilobus',
              NA,208,17.4,NA,7.0,NA,'Surface',NA)
dpthTbl <- rbind(dpthTbl,trilobus)

# Add depth ranges from Watkins et al. 1996

# G menardii ranges from 20-60m, ALD of 25 read off fig 6c
menardii <- grep('Globorotalia.menardii',dpthTbl$Species)
dpthTbl$ALD[menardii] <- 25
dpthTbl$VD[menardii] <- 10
dpthTbl$MaxAbund[menardii] <- 3
dpthTbl$DepthHabitat[menardii] <- 'Surface'

# G tumida ranges lives in upper 80m, ALD of 50 read off fig 4b
tumida <- c('Globorotalia.tumida',
            NA,30,50,NA,15,NA,'Surface',NA)

# G conglobatus also lives in upper 80m, ALD of 50 (fig 5c)
conglob <- c('Globigerinoides.conglobatus',
             NA,3,50,NA,15,NA,'Surface',NA)

# G hexagonus is sub-thermocline, 100-200m
hexag <- c('Globorotaloides.hexagonus',
           NA,2,150,NA,75,NA,'Subsurface',NA)
dpthTbl <- rbind(dpthTbl,tumida,conglob,hexag)

# Synonymize according to Microtax for congruence with occurrence database
# Both synonyms of digitata (Beella and Hastigerinella) are present in depth dataset
# use only Beella to match occurrences of the species
dpthTbl$Species <- gsub('Globorotalia.scitula','Hirsutella.scitula',dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.inflata','Globoconella.inflata',dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.hirsuta','Hirsutella.hirsuta',dpthTbl$Species)
dpthTbl$Species <- gsub('Berggrenia.pumillio','Berggrenia.pumilio',dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.crassaformis','Truncorotalia.crassaformis',dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.menardii','Menardella.menardii',dpthTbl$Species)
dpthTbl$Species <- gsub('Globigerinoides.trilobus','Trilobatus.trilobus',dpthTbl$Species)
dpthTbl$Species <- gsub('Globigerinoides.ruber.white','Globigerinoides.ruber',dpthTbl$Species)
dpthTbl$Species <- gsub('Globigerinoides.tenellus','Globoturborotalita.tenella',dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.truncatulinoides','Truncorotalia.truncatulinoides',
                        dpthTbl$Species)

# Add depth data to sp-level dataset (except for 1st column, which would duplicate sp names)
spp <- sppDat$species
spp <- gsub(' ','.',spp)
wOccs <- dpthTbl$Species %in% spp
dpthTbl <- dpthTbl[wOccs,]
newCols <- colnames(dpthTbl)[-1]
sppDat[,newCols] <- NA
spPos <- match(dpthTbl$Species, spp)
sppDat[spPos,newCols] <- dpthTbl[,-1]
# write.csv(sppDat, 'Data/foram_spp_data.csv', row.names = FALSE)

# cut down file size
irrel <- c('site','hole','core','section','sample.top','sample.type',
           'preservation','processing','Macro.micro')
cols2toss <- colnames(occ) %in% c(traits, irrel)
occ <- occ[,!cols2toss]

##########################################################
# bin to 4 ky intervals to match GCM data

brk <- seq(0, 796, by=4)
bins <- data.frame(t=brk, b=brk+4, mid=brk+2)

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
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
palCoords <- occ[,c('pal.long','pal.lat')]
pts <- SpatialPointsDataFrame(palCoords, data = occ, proj4string = CRS(llPrj))

occ$cell_number <- cellFromXY(rEmpt, pts) 
occ$centroid_lat <- occ$centroid_long <- NA
occ[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, occ$cell_number)

# remove records with imprecise age estimates
unconstrnd <- which(occ$age.model %in% c('Berggren1977','Ericson1968','GTSBlow1969'))
occ <- occ[-unconstrnd,]

# remove species too rare to reconstruct niches
nmFreq <- table(occ$species)
rare <- names( which(nmFreq <6) )
occ <- occ[! occ$species %in% rare,]

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

