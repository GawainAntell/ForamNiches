library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

# Data import -------------------------------------------------------------

tRes <- 8
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

# Read in occurrence data
occ <- read.csv('Data/foram_data_800ka_IF191204.csv', stringsAsFactors=FALSE)

# easier to work on scale of ka not ma
occ$age <- occ$age*1000

# Extract core range-through data -----------------------------------------

# There are sometimes 2 or 3 coreID's per hole, with .a/.b/.c suffix,
# seemingly because of differences in data entry rules.
# There are also typos of dashes and underscores mixed up.
occ$coreID <- gsub('-','_',occ$coreID)
occ$coreUniq <- gsub('181_1119B','181_1119',occ$coreID)
rmSuffix <- function(txt){
  strsplit(txt, '[.]')[[1]][1]
}
occ$coreUniq <- sapply(occ$coreUniq, rmSuffix)
corsUniq <- unique(occ$coreUniq)

# Note that 'core' is not a unique identifier, and many 'hole' values are NA.
# The coreID is nearly unique, but there are mistakes.
# Make sure that the coordinates match or are very close (< 100m)
# and if not, then create a new coreID.
problm <- data.frame()
for (cr in corsUniq){
  crRows <- which(occ$coreUniq==cr)
  crDf <- occ[crRows,]
  mdrnCoords <- c('longitude','latitude')
  dupes <- duplicated(crDf[,mdrnCoords])
  uniqCrs <- crDf[!dupes,]
  if (nrow(uniqCrs) > 1){
    crPts <- SpatialPoints(uniqCrs[,mdrnCoords])
    gcdists <- spDists(crPts)
    far <- any(gcdists > 0.1)
    if (far){
      problm <- rbind(problm, uniqCrs)
    }
  }
}
write.csv(problm, paste0('Data/problem_core_IDs_for_IF_',day,'.csv'))
badID <- unique(problm$coreUniq)

cors2use <- ! corsUniq %in% badID
corsUniqCln <- corsUniq[cors2use]
corInfo <- function(cr){
  crRows <- which(occ$coreUniq==cr)
  crDf <- occ[crRows,]
  ageRng <- range(crDf$age)
  c(crDf$longitude[1], crDf$latitude[1], ageRng)
}
corMat <- sapply(corsUniqCln, corInfo)
corAtts <- data.frame(t(corMat))
colnames(corAtts) <- c('long','lat','lad','fad')

# rasterize to the resolution of GCM data
rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
corPts <- SpatialPoints(corAtts[,c('long','lat')], CRS(llPrj))
corAtts$cell <- cellFromXY(rEmpt, corPts) 

corNm <- paste0('Data/core_rangethrough_data_',day,'.csv')
write.csv(corAtts, corNm)

# Absence data --------------------------------------------------------------

# Any original abundance of NA is suspect. Surprisingly, the abundance
# for some of these records is > 0; those records aren't reliable.
naAbund <- is.na(occ$orig.abundance)
# summary(occ$abundance[naAbund])
occ <- occ[!naAbund,]

uniqAbund <- unique(occ$orig.abundance)
numAbund <- as.numeric(uniqAbund)
nonNum <- uniqAbund[is.na(numAbund)]
nonNum

# Some original abundance codings are ambiguous. Data in these categories should be omitted:
bad <- c('-', '/', '*', 'rw?', '#0.0', '?1', 'p?')
badBool <- occ$orig.abundance %in% bad
occ <- occ[!badBool,]

# There are ~150 original abundance codings that imply a species presence, but abundance is 0.
# Modify the abundance so that these records are not discarded.
ok <- setdiff(nonNum, bad)
okBool <- occ$orig.abundance %in% ok
occ$abundance[okBool] <- 1

# Lastly, there is the case that abundance is recorded as 0; these are presumably absences.
zeroBool <- occ$abundance==0
occ <- occ[!zeroBool,]

# Species-level data ------------------------------------------------------

# remove T. cavernula and excelsa, which originated <800 ka
tooYng <- which(occ$earliest < .8)
occ <- occ[-tooYng,]

# excise sp traits into a separate, sp-level file to call later in analysis
traits <- c('morphDes','spinose','structure','aperture','latest','earliest',
            'mph','mphRef','eco','ecoRef','geo','geoRef',
            'log.area','Approximate.diameter') # 'Average.Area.microm2',
sppDat <- occ[,c('species',traits)]
dupes <- duplicated(sppDat$species)
sppDat <- sppDat[!dupes,]

# Depth ranges ------------------------------------------------------------
# add modern depth ranges to sp-level file

# Read in modern depth ranges: Rebotim et al 2017, table 4
dpthTbl <- read.table('Data/Rebotim_et_al_depth_ranges.txt',
                      header=TRUE, stringsAsFactors=FALSE)
abund <- dpthTbl$MaxAbund
# remove single-letter code that follows the numeric abund value
penult <- nchar(abund) - 1
dpthTbl$MaxAbund <- substr(abund, 1, penult)
dpthTbl$MaxAbund <- as.numeric(dpthTbl$MaxAbund)

# Add G trilobus data from Loncaric et al. 2006
trilobus <- c('Globigerinoides.trilobus',
              NA,208,17.4,NA,7.0,NA,'Surface',NA)
dpthTbl <- rbind(dpthTbl,trilobus)

# Add G (H) theyeri data from Schiebel et al. 2004
# found only from 40-60 m, very rare
theyeri <- c('Hirsutella.theyeri',
             NA,0.2,50,NA,15,NA,'Surface',NA)
dpthTbl <- rbind(dpthTbl,theyeri)

# Add depth ranges from Schiebel and Hemleben 2017

# B digitata is subsurface
digitata <- grep('Beella.digitata',dpthTbl$Species)
dpthTbl$DepthHabitat[digitata] <- 'Subsurface'
dpthTbl$ALD[digitata] <- 150
dpthTbl$VD[digitata] <- 75

# S dehiscens is sub-thermocline, subsurface
dehiscens <- c('Sphaeroidinella.dehiscens',
               NA,NA,150,NA,75,NA,'Subsurface',NA)

# G ungulata has same distribution as G menardii; 
# density up to 2 per cubic m
ungulata <- c('Globorotalia.ungulata',
         NA,NA,46.2,NA,25.2,NA,'Surface',NA)
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

# G tumida ranges lives in upper 80m, ALD of 50 read off fig 4b
tumida <- c('Globorotalia.tumida',
            NA,30,50,NA,15,NA,'Surface',NA)

# G conglobatus also lives in upper 80m, ALD of 50 (fig 5c)
conglob <- c('Globigerinoides.conglobatus',
             NA,3,50,NA,15,NA,'Surface',NA)

# G hexagonus is sub-thermocline, 100-200m
hexag <- c('Globorotaloides.hexagonus',
           NA,2,150,NA,75,NA,'Subsurface',NA)

# G adamsi is sub-termocline, 100-200 m
# Schiebel et al. 2004 sampled adamsi mainly 40-100m,
# but didn't sample below 100 m
adamsi <- c('Globigerinella.adamsi',
            NA,2,125,NA,75,NA,'Surface.subsurface',NA)

# G conglomerata is subsurface, 60-100m, ALD ca 80 (fig 4a)
# Schiebel et al. 2004: G conglomerata is surface, ALD=42.5, VD=30.8
glom <- c('Globoquadrina.conglomerata',
          NA,2.4,42.5,NA,30.8,NA,'Surface.subsurface',NA)

dpthTbl <- rbind(dpthTbl,tumida,conglob,hexag,adamsi,glom)

# Synonymize according to Microtax for congruence with occurrence database
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
spDatNm <- paste0('Data/foram_spp_data_',day,'.csv')
write.csv(sppDat, spDatNm, row.names = FALSE)

# cut down file size
irrel <- c('site','hole','core','section','sample.top','sample.type',
           'preservation','processing','Macro.micro')
cols2toss <- colnames(occ) %in% c(traits, irrel)
occ <- occ[,!cols2toss]

# Time binning ------------------------------------------------------------

brk <- seq(0, 800-tRes, by=tRes)
bins <- data.frame(t=brk, b=brk+tRes, mid=brk+tRes/2)

# age 'zero' = 1950, and some observations are more recent (i.e. negative age)
bins$t[1] <- -0.1

for (r in 1:nrow(occ)){
  a <- occ$age[r]
  bin <- which(bins$t < a & bins$b >= a)
  aBin <- bins$mid[bin]
  occ$bin[r] <- aBin
}

# No breeding populations at shallow depths, so any occurrences there are suspect
shall <- which(occ$water.depth < 150)
occ <- occ[-shall,]

# Rasterize to the resolution of GCM data
palCoords <- occ[,c('pal.long','pal.lat')]
pts <- SpatialPointsDataFrame(palCoords, data = occ, proj4string = CRS(llPrj))
occ$cell_number <- cellFromXY(rEmpt, pts) 
occ$centroid_lat <- occ$centroid_long <- NA
occ[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, occ$cell_number)

# remove records with imprecise age estimates
ageMods <- c('Berggren1977','Ericson1968','GTSBlow1969','Raffi2006')
unconstrnd <- union(which(occ$age.model %in% ageMods), 
                    which(occ$rng.age > tRes/1000))
occ <- occ[-unconstrnd,]

spp <- unique(occ$species)

# Shorten to unique occs --------------------------------------------------
# subset by bin, then species, then find unique species-cells combinations

# omit duplicate cell-species combinations
slcSpSmry <- function(dat, sp){
  spSlc <- dat[dat$species==sp,] 
  uniqRecord <- function(x){ 
    first <- which(spSlc$cell_number==x)[1]
    spSlc[first,] 
  }
  tempList <- lapply(unique(spSlc$cell_number), uniqRecord)
  do.call('rbind', tempList)
}

# apply slcSpSmry over species in a stage
slcSmry <- function(dat, bin){
  slc <- dat[dat$bin==bin,]
  tempList <- lapply(unique(slc$species), slcSpSmry, dat=slc)
  do.call('rbind', tempList)
}

nCore <- detectCores() - 1
pt1 <- proc.time()
registerDoParallel(nCore)
stageDfList <- foreach(bin=bins$mid) %dopar% 
  slcSmry(bin=bin, dat=occ)
stopImplicitCluster()
pt2 <- proc.time()
pt2-pt1

fin <- do.call('rbind', stageDfList)

# Cleaning ----------------------------------------------------------------

# Subset occurrences such that each sp has >5 occs per bin
saveRows <- function(sp, bin, df){
  spRows <- which(df$species==sp & df$bin==bin)
  if (length(spRows)>5){
    spRows
  }
}  
keepRowsL <- sapply(spp, function(x){
  sapply(bins$mid, function(b) saveRows(sp=x, bin=b, df=fin))
} 
)
keepRows <- unlist(keepRowsL)
fin <- fin[keepRows,]

dfNm <- paste('Data/foram-uniq-occs_latlong_',tRes, 'ka_', day, '.csv', sep='')
write.csv(fin, file=dfNm, row.names=FALSE)
