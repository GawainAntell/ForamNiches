library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)

# Data import -------------------------------------------------------------

tRes <- 8
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

# raster at the resolution of GCM data
rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Read in occurrence data
occ <- read.csv('Data/foram_data_800ka_IF191204.csv', stringsAsFactors=FALSE)

# easier to work on scale of Ka not Ma
occ$age <- occ$age*1000

# typo
occ$species <- gsub('Globototalia anfracta','Globorotalia anfracta', occ$species)

# save names of palaeocoordinate columns to call later
coordCols <- c('pal.long','pal.lat')

# No breeding populations at shallow depths, so any occurrences there are suspect
shall <- which(occ$water.depth < 100)
occ <- occ[-shall,]

# Clean core IDs ----------------------------------------------------------

# There are sometimes 2 or 3 coreID's per hole, with .a/.b/.c suffix,
# seemingly because of differences in data entry rules.
# There are also typos of dashes and underscores mixed up.
occ$coreUniq <- occ$coreID
occ$coreUniq <- gsub('-','_',occ$coreUniq)
occ$coreUniq <- gsub('181_1119B', '181_1119', occ$coreUniq)
rmSuffix <- function(txt){
  strsplit(txt, '[.]')[[1]][1]
}
occ$coreUniq <- sapply(occ$coreUniq, rmSuffix)

# Note that 'core' is not a unique identifier, and many 'hole' values are NA.
# The coreID is nearly unique, but there are mistakes.
# Make sure that the coordinates match or are very close (< 100m)
# and if not, then create a new coreID.

# export a spreadsheet of problem IDs to fix in original database
# corsUniq <- unique(occ$coreUniq)
# problm <- data.frame()
# for (cr in corsUniq){
#   crRows <- which(occ$coreUniq==cr)
#   crDf <- occ[crRows,]
#   mdrnCoords <- c('longitude','latitude')
#   dupes <- duplicated(crDf[,mdrnCoords])
#   uniqCrs <- crDf[!dupes,]
#   if (nrow(uniqCrs) > 1){
#     crPts <- SpatialPoints(uniqCrs[,mdrnCoords])
#     gcdists <- spDists(crPts)
#     far <- any(gcdists > 0.1)
#     if (far){
#       problm <- rbind(problm, uniqCrs)
#     }
#   }
# }
# write.csv(problm, paste0('Data/problem_core_IDs_for_IF_',day,'.csv'))
# badID <- unique(problm$coreUniq)

# Cores RC8C80, RC8C81, RC8C82, RC8C83, RC9C98 and RC9C99
# each contain duplicate sets of entries, with different coordinates 
# for each set. It's unclear which is correct.
problm1 <- occ$coreUniq %in% 
  c('RC8C80','RC8C81','RC8C82','RC8C83','RC9C98','RC9C99')
occ <- occ[!problm1,]

# Cores NA87_22 and NA87_22.a are the same, but the coordinates 
# differ by 0.1 degrees. Use the same average position for both.
problm2 <- which(occ$coreID=='NA87-22')
problm3 <- which(occ$coreID=='NA87-22.a')
examp1 <- occ[problm2[1], coordCols]
examp2 <- occ[problm3[1], coordCols]
newLong <- mean(c(examp1$pal.long, examp2$pal.long))
newLat <- mean(c(examp1$pal.lat, examp2$pal.lat))
occ$pal.long[c(problm2, problm3)] <- newLong
occ$pal.lat[ c(problm2, problm3)] <- newLat

# Extract core range-through data -----------------------------------------

corsUniqCln <- unique(occ$coreUniq)
corInfo <- function(cr){
  crRows <- which(occ$coreUniq==cr)
  crDf <- occ[crRows,]
  ageRng <- range(crDf$age)
  c(crDf$longitude[1], crDf$latitude[1], ageRng)
}
corMat <- sapply(corsUniqCln, corInfo)
coreAtts <- data.frame(t(corMat))
colnames(coreAtts) <- c('long','lat','lad','fad')

corNm <- paste0('Data/core_rangethrough_data_',day,'.csv')
write.csv(coreAtts, corNm)

# Absence data --------------------------------------------------------------

# Any original abundance of NA is suspect. Surprisingly, the abundance
# for some of these records is > 0; those records aren't reliable.
naAbund <- which(is.na(occ$orig.abundance) & occ$abundance > 0)
# summary(occ$abundance[naAbund])
occ <- occ[-naAbund,]

# Some original abundance codings are ambiguous
uniqAbund <- unique(occ$orig.abundance)
numAbund <- suppressWarnings(as.numeric(uniqAbund))
nonNum <- uniqAbund[is.na(numAbund)]
nonNum
# Data in these categories should be omitted:
bad <- c('-', '/', '*', 'rw?', '#0.0', '?1', 'p?')
badBool <- occ$orig.abundance %in% bad
occ <- occ[!badBool,]

# Some original abundance codings imply a species presence, but abundance is 0.
# Modify the abundance so that these records are not discarded.
ok <- setdiff(nonNum, c(bad, NA))
not0 <- which(occ$orig.abundance %in% ok & occ$abundance==0)
occ$abundance[not0] <- 1

# Lastly, there is the case that abundance is recorded as 0; these are presumably absences.
# Absences were already used in the calculation of the sampling universe
# (core range-through data), so omit them now. Analyses will be on presence data.
zeroBool <- occ$abundance==0
occ <- occ[!zeroBool,]

# Species-level data ------------------------------------------------------

# remove T. cavernula and excelsa, which originated <800 ka
tooYng <- which(occ$earliest < .8)
occ <- occ[-tooYng,]

# excise sp traits into a separate, sp-level file to call later in analysis
traits <- c('morphDes','spinose','Macro.micro','structure','aperture',
            'latest','earliest','mph','mphRef','eco','ecoRef','geo','geoRef',
            'log.area','Approximate.diameter') # 'Average.Area.microm2',
sppDat <- occ[,c('species',traits)]
dupes <- duplicated(sppDat$species)
sppDat <- sppDat[!dupes,]

# add modern depth ranges to sp-level file

# Read in modern depth ranges: Rebotim et al 2017, table 4
dpthTbl <- read.table('Data/Rebotim_et_al_depth_ranges.txt',
                      header=TRUE, stringsAsFactors=FALSE)
abund <- dpthTbl$MaxAbund
# remove single-letter code that follows the numeric abund value
penult <- nchar(abund) - 1
dpthTbl$MaxAbund <- substr(abund, 1, penult)
dpthTbl$MaxAbund <- as.numeric(dpthTbl$MaxAbund)
# add a category for depth reference
dpthTbl$ref <- 'Rebotim et al. 2017'

# Schiebel et al. 2004 also confirm that minuta is within top 60m
minut <- which(dpthTbl$Species=='Globigerinita.minuta')
dpthTbl$DepthHabitat[minut] <- 'Surface'

# Specify depth range based on Retailleau et al. 2011
uvula <- which(dpthTbl$Species=='Globigerinita.uvula')
dpthTbl$DepthHabitat[uvula] <- 'Surface'
dpthTbl$ALD[uvula] <- 26
dpthTbl$ref[uvula] <- 'Retailleau et al. 2011'

# Hull et al. did surveys of H digitata; 280-358 m
hdig <- which(dpthTbl$Species=='Hastigerinella.digitata')
dpthTbl$DepthHabitat[hdig] <- 'Subsurface'
dpthTbl$ref[hdig] <- 'Hull et al. 2011'

# Schiebel et al. sampled ALD of 60.5 m, but quite variable
# Schmuker and Schiebel 2002 also got depth of 60-80 m
# but Schiebel and Hemleben 2017 considered the species shallow
anfrac <- which(dpthTbl$Species=='Dentigloborotalia.anfracta')
dpthTbl$DepthHabitat[anfrac] <- 'Surface'
dpthTbl$ref[anfrac] <- c('Schiebel et al. 2004, Rebotim et al. 2017')

# Add G trilobus data from Loncaric et al. 2006
trilobus <- c('Globigerinoides.trilobus',
              NA,208,17.4,NA,7.0,NA,'Surface',NA,'Loncaric et al. 2006')
dpthTbl <- rbind(dpthTbl,trilobus)

# Add G (H) theyeri data from Schiebel et al. 2004
# found only from 40-60 m, very rare
theyeri <- c('Hirsutella.theyeri',
             NA,0.2,50,NA,15,NA,'Surface',NA,'Schiebel et al. 2004')
dpthTbl <- rbind(dpthTbl,theyeri)

# Add depth ranges from Schiebel and Hemleben 2017

# B digitata is subsurface
digitata <- grep('Beella.digitata',dpthTbl$Species)
dpthTbl$DepthHabitat[digitata] <- 'Subsurface'
dpthTbl$ALD[digitata] <- 150
dpthTbl$VD[digitata] <- 75
dpthTbl$ref[digitata] <- 'Schiebel and Hemleben 2017'

# S dehiscens is sub-thermocline, subsurface
dehiscens <- c('Sphaeroidinella.dehiscens',NA,NA,150,NA,75,NA,
               'Subsurface',NA,'Schiebel and Hemleben 2017')

# G ungulata has same distribution as G menardii; 
# density up to 2 per cubic m
ungulata <- c('Globorotalia.ungulata',NA,NA,46.2,NA,25.2,NA,
              'Surface',NA,'Schiebel and Hemleben 2017')
dpthTbl <- rbind(dpthTbl, dehiscens, ungulata)

# Add depth ranges from Watkins et al. 1996

# G menardii ranges from 20-60m, ALD of 25 read off fig 6c
# Schiebel et al. 2004 called menardii an upwelling indicator sp
# and report densities that give ALD=58.6 and Vd=21.9 for one leg,
# ALD=46.2 and VD=25.2 for another leg (with higher concentrations)
menardii <- grep('Globorotalia.menardii',dpthTbl$Species)
dpthTbl$ALD[menardii] <- 46.2
dpthTbl$VD[menardii] <- 25.2
dpthTbl$MaxAbund[menardii] <- 51
dpthTbl$DepthHabitat[menardii] <- 'Surface'
dpthTbl$ref[menardii] <- 'Watkins et al. 1996'

# G tumida ranges lives in upper 80m, ALD of 50 read off fig 4b
tumida <- c('Globorotalia.tumida',NA,30,50,NA,15,NA,
            'Surface',NA,'Watkins et al. 1996')

# G conglobatus also lives in upper 80m, ALD of 50 (fig 5c)
conglob <- c('Globigerinoides.conglobatus',NA,3,50,NA,15,NA,
             'Surface',NA,'Watkins et al. 1996')

# G hexagonus is sub-thermocline, 100-200m
hexag <- c('Globorotaloides.hexagonus',NA,2,150,NA,75,NA,
           'Subsurface',NA,'Watkins et al. 1996')

# G adamsi is sub-thermocline, 100-200 m
# Schiebel et al. 2004 sampled adamsi mainly 40-100m,
# but didn't sample below 100 m
adamsi <- c('Globigerinella.adamsi',NA,2,125,NA,75,NA,
            'Surface.subsurface',NA,
            'Watkins et al. 1996, Schiebel et al. 2004')

# G conglomerata is subsurface, 60-100m, ALD ca 80 (fig 4a)
# Schiebel et al. 2004: G conglomerata is surface, ALD=42.5, VD=30.8
glom <- c('Globoquadrina.conglomerata',
          NA,2.4,42.5,NA,30.8,NA,'Surface.subsurface',NA,
          'Watkins et al. 1996, Schiebel et al. 2004')

dpthTbl <- rbind(dpthTbl,tumida,conglob,hexag,adamsi,glom)

# Re-assign crassaformis as mid-depth as a compromise among conflicting studies. 
# Meilland et al (2019) and Ezard et al (2015) consider the species sub-thermocline,
# but Rebotim et al (2017) and Schiebel (2017) consider it shallow.
crass <- which(dpthTbl$Species=='Globorotalia.crassaformis')
dpthTbl$DepthHabitat[crass] <- 'Surface.subsurface'
dpthTbl$ref[crass] <- 
  c('Ezard et al. 2015, Rebotim et al. 2017, Schiebel 2017, Meilland et al. 2019')

# Synonymize according to Microtax for congruence with occurrence database
dpthTbl$Species <- gsub('Dentigloborotalia.anfracta','Globorotalia.anfracta',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.scitula','Hirsutella.scitula',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.inflata','Globoconella.inflata',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.hirsuta','Hirsutella.hirsuta',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Berggrenia.pumillio','Berggrenia.pumilio',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.crassaformis','Truncorotalia.crassaformis',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.menardii','Menardella.menardii',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globigerinoides.trilobus','Trilobatus.trilobus',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globigerinoides.ruber.white','Globigerinoides.ruber',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globigerinoides.tenellus','Globoturborotalita.tenella',
                        dpthTbl$Species)
dpthTbl$Species <- gsub('Globorotalia.truncatulinoides','Truncorotalia.truncatulinoides',
                        dpthTbl$Species)

# Add depth data to sp-level dataset (except for 1st column, which would duplicate sp names).
# There are species with occurrences for which the habitat is still NA in this object, 
# but these species are omitted in subsequent scripts due to insufficient sample size.
sppDot <- sppDat$species
sppDot <- gsub(' ','.',sppDot)
wOccs <- dpthTbl$Species %in% sppDot
dpthTbl <- dpthTbl[wOccs,]
newCols <- colnames(dpthTbl)[-1]
sppDat[,newCols] <- NA
spPos <- match(dpthTbl$Species, sppDot)
sppDat[spPos,newCols] <- dpthTbl[,-1]
spDatNm <- paste0('Data/foram_spp_data_', day, '.csv')
write.csv(sppDat, spDatNm, row.names = FALSE)

# cut down file size
irrel <- c('orig.species','abundance','orig.abundance','abun.units',
           'age.st','age.en','AM.type','rel.abun','round.age',
           'site','hole','core','section','sample.top','sample.type',
           'total.IDd','preservation','processing','Comments',
           'Average.Area..µm2.','log.area','Approximate.diameter')
cols2toss <- colnames(occ) %in% c(traits, irrel)
occ <- occ[,!cols2toss]

# Bin and rasterize -------------------------------------------------------

# remove records with very imprecise age estimates
ageMods <- c('Berggren1977','Ericson1968','GTSBlow1969','Raffi2006')
unconstrnd <- union( which(occ$age.model %in% ageMods), 
                     which(occ$rng.age > tRes/1000)
                     )
occ <- occ[-unconstrnd,]
# all remaining records have high enough precision for binning
  # unique(occ$rng.age)*1000

brk <- seq(0, 800-tRes, by=tRes)
bins <- data.frame(t=brk, b=brk+tRes, mid=brk+tRes/2)

# age 'zero' = 1950, and some observations are more recent (i.e. negative age)
bins$t[1] <- -0.1

for (r in 1:nrow(occ)){
  a <- occ$age[r]
  bin <- which(bins$t < a & bins$b >= a)
  occ$bin[r] <- bins$mid[bin]
}

# Rasterize to the resolution of GCM data
pts <- SpatialPointsDataFrame(occ[, coordCols], 
                              data = occ, proj4string = CRS(llPrj))
occ$cell_number <- cellFromXY(rEmpt, pts) 
occ[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, occ$cell_number)

# Shorten to unique occs --------------------------------------------------
# subset by bin, then species, then find unique species-cells combinations

# omit duplicate cell-species combinations
slcSpSmry <- function(dat, sp){
  spSlc <- dat[dat$species==sp,] 
  dupes <- duplicated(spSlc$cell_number)
  spSlc[!dupes,]
}

# apply slcSpSmry over species in a stage
slcSmry <- function(dat, bin){
  slc <- dat[dat$bin==bin,]
  tempList <- lapply(unique(slc$species), slcSpSmry, dat=slc)
  do.call('rbind', tempList)
}

nCore <- detectCores() - 1
registerDoParallel(nCore)
stageDfList <- foreach(bin=bins$mid) %dopar% 
  slcSmry(bin=bin, dat=occ)
stopImplicitCluster()

df <- do.call('rbind', stageDfList)
row.names(df) <- as.character(1:nrow(df))

# Combine enviro and spp data ---------------------------------------------

# Calculate the Average Living Depth of species in the dataset:
spp <- unique(df$species) 
row.names(sppDat) <- sppDat$species
sppDat <- sppDat[spp,]
sppDat$ALD <- as.numeric(sppDat$ALD)
zones <- unique(sppDat$DepthHabitat)
for (z in zones){
  zRows <- which(sppDat$DepthHabitat==z)
  avg <- mean(sppDat$ALD[zRows],na.rm=T)
  avg <- round(avg)
  print(paste(z, avg))
}
# [1] "Subsurface 174"
# [1] "Surface.subsurface 88"
# [1] "Surface 42"

# Find equivalent AOGCM depth levels from which to extract temp: 0, 40, 78, 164 m
# The lower 3 correspond to surface, surface-subsurface, and subsurface species
# and the values are based on the mean Average Living Depth of species in the dataset.
dpths <- c(1, 4, 6, 8)

source('raster_brick_import_fcn.R')

modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
bins <- unique(df$bin)
envNm <- 'ann_temp_ym_dpth'
envNmShort <-  paste(strsplit(envNm,'_')[[1]][2:3], collapse='_')
allEnvNm <- paste(envNmShort, c('0m','surf','surfsub','sub'), sep='_')

addEnv <- function(bin, dat, mods, binCol, cellCol, prj, 
                   env, envCols, dpths){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCells <- slc[,cellCol]
  slcEnv <- getBrik(bin = bin, envNm = env, mods = mods)
  
  envVals <- raster::extract(slcEnv[[env]], slcCells)
  # Rows = points of extraction, columns = depth layers  
  envVals <- envVals[,dpths]
  slc[,envCols] <- envVals
  
  # Infer environment if it's missing and some of the adjacent 9 cells have values
  # This reduces the number of occs without enviro from >800 to 200
  naRows <- apply(envVals, 1, function(x) any(is.na(x)) )
  if (sum(naRows)>0){
    naCoords <- slc[naRows,c('centroid_long','centroid_lat')]
    # distance to corner cell is 196 km for 1.25-degree resolution (~111 km/degree)
    fillr <- raster::extract(slcEnv[[env]], naCoords, 
                             buffer = 200*1000, fun = mean)
    slc[naRows, envCols] <- fillr[,dpths]
  }
  slc
}

pkgs <- c('sp','raster') 
registerDoParallel(nCore)
sppEnv <- foreach(bin=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  addEnv(bin = bin, dat = df, mods = modId, binCol = 'bin', cellCol = 'cell_number', 
         prj = llPrj, env = envNm, envCols = allEnvNm, dpths = dpths)
stopImplicitCluster()

# remove records where environment is unknown
naRows <- apply(sppEnv[,allEnvNm], 1, 
                function(r) any(is.na(r))
                )
sppEnv <- sppEnv[!naRows,]

obsNm <- paste0('Data/foram_uniq_occs_latlong_8ka_', day, '.csv')
write.csv(sppEnv, obsNm, row.names = FALSE)

# Enviro at unique sampled sites ------------------------------------------

# read in modified function to reconstruct paleo-coordinates
source('paleocoords_fcn.R')

# get the sampling universe (env at range-through core sites) per bin
sampEnv <- function(b){
  # for each bin, find out which cores range through the 8ky interval
  inBin <- which(coreAtts$fad > (b-4) & coreAtts$lad <= (b+4))
  slc <- coreAtts[inBin,]
  
  # rotate lat-long coordinates to paleo-locations
  palCoords <- pal.coord(slc, 'MATTHEWS2016', Ma = b/1000)
  slc$pal.long <- palCoords$paleolng
  slc$pal.lat  <- palCoords$paleolat
  # use modern coordinates for the few points that can't be reconstructed
  naPts <- apply(slc[, coordCols], 1, 
                 function(r) any(is.na(r))
  )
  slc[naPts, coordCols] <- slc[naPts, c('long','lat')]
  
  # convert (paleo) lat-long to cell numbers, to shorten data to unique cells
  pts <- SpatialPointsDataFrame(slc[, coordCols], 
                                data = slc, proj4string = CRS(llPrj))
  slc$cell <- cellFromXY(rEmpt, pts) 
  
  dupeCell <- duplicated(slc$cell)
  slc <- slc[!dupeCell,]
  
  # get the env values at the core locations
  slcEnv <- getBrik(bin = b, envNm = envNm, mods = modId)
  envVals <- raster::extract(slcEnv[[envNm]], slc$cell)
  # Rows = points of extraction, columns = depth layers  
  envVals <- envVals[,dpths]
  colnames(envVals) <- allEnvNm
  sampled <- cbind(b, cell = slc$cell, envVals)
  
  # interpolate missing values from the adjacent 9 cells
  naRows <- apply(envVals, 1, function(x) any(is.na(x)) )
  if (sum(naRows) > 0){
    naCoords <- slc[naRows, coordCols]
    
    # distance to corner cell is 196 km for 1.25-degree resolution (~111 km/degree)
    fillr <- raster::extract(slcEnv[[envNm]], naCoords, 
                             buffer = 200*1000, fun = mean)
    sampled[naRows, allEnvNm] <- fillr[,dpths]
  }
  sampled
}

# NB: The average distance between modern and reconstructed coordinates for the
# oldest data (800 Ka) is only 26 km - far less than the length of a grid cell.
  # b <- 800
  # inBin <- which(coreAtts$fad > (b-4) & coreAtts$lad <= (b+4))
  # slc <- coreAtts[inBin,]
  # palCoords <- pal.coord(slc, 'MATTHEWS2016', Ma = b/1000)
  # slc$pal.long <- palCoords$paleolng
  # slc$pal.lat  <- palCoords$paleolat
  # mdrn <- SpatialPointsDataFrame(slc[, c('long','lat')], 
  #                               data = slc, proj4string = CRS(llPrj))
  # old  <- SpatialPointsDataFrame(slc[, coordCols], 
  #                               data = slc, proj4string = CRS(llPrj))
  # dists <- pointDistance(mdrn, old, longlat=TRUE)
  # mean(dists)/1000
  # > [1] 26.28166

registerDoParallel(nCore)
samp <- foreach(b=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  sampEnv(b=b)
stopImplicitCluster()

# remove records where environment is unknown
naSampRows <- apply(samp[,allEnvNm], 1, 
                    function(r) any(is.na(r))
                    )
samp <- samp[!naSampRows,]

sampNm <- paste0('Data/samp_uniq_occs_latlong_8ka_', day, '.csv')
write.csv(samp, sampNm, row.names = FALSE)
