library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(phytools)
library(paleoPhylo)
library(ape)
library(ggplot2)

doTrunc <- FALSE

# Data import -------------------------------------------------------------

tRes <- 8
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

# Read in occurrence data
occ <- read.csv('Data/foram_data_800ka_IF191204.csv', stringsAsFactors=FALSE)

# easier to work on scale of ka not ma
occ$age <- occ$age*1000

# Extract core range-through data -----------------------------------------

# No breeding populations at shallow depths, so any occurrences there are suspect
shall <- which(occ$water.depth < 150)
occ <- occ[-shall,]

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
# write.csv(problm, paste0('Data/problem_core_IDs_for_IF_',day,'.csv'))
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
coreAtts <- data.frame(t(corMat))
colnames(coreAtts) <- c('long','lat','lad','fad')

# rasterize to the resolution of GCM data
rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
corPts <- SpatialPoints(coreAtts[,c('long','lat')], CRS(llPrj))
coreAtts$cell <- cellFromXY(rEmpt, corPts) 

corNm <- paste0('Data/core_rangethrough_data_',day,'.csv')
write.csv(coreAtts, corNm)

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
palCoords <- occ[,c('pal.long','pal.lat')]
pts <- SpatialPointsDataFrame(palCoords, data = occ, proj4string = CRS(llPrj))
occ$cell_number <- cellFromXY(rEmpt, pts) 
occ[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, occ$cell_number)

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
registerDoParallel(nCore)
stageDfList <- foreach(bin=bins$mid) %dopar% 
  slcSmry(bin=bin, dat=occ)
stopImplicitCluster()

df <- do.call('rbind', stageDfList)

# Phylo data compatability ------------------------------------------------

# Limit analysis to species included in tree from Aze et al. 2011
trRaw <- read.csv('Data/Aze_et_al_2011_bifurcating_tree_data.csv', stringsAsFactors=FALSE)
trRaw$Species.name <- gsub('Globigerinoides sacculifer', 'Trilobatus sacculifer', trRaw$Species.name)
trRaw$Species.name <- gsub('Globigerinoides trilobus', 'Trilobatus trilobus', trRaw$Species.name)
tr_paleo <- with(trRaw, 
                 as.paleoPhylo(Species.code, Ancestor.code, Start.date, End.date)
)
tr <- buildApe(tr_paleo)
sppCodes <- tr$tip.label
rowOrdr <- match(sppCodes, trRaw$Species.code)
tr$tip.label <- trRaw$Species.name[rowOrdr]

# 9 species are not present in the phylogeny, mostly because microporiferate
sppAll <- unique(df$species)
lostSpp <- setdiff(sppAll, tr$tip.label)

# Beella megastoma is arguably the same species as B. digitata,
# and Truncorotalia crassula is arguably senior synonym to crassaformis
# (Schiebel and Hemleben 2017). The depth ranges for both are unknown.
lostSpp <- c(lostSpp, 'Beella megastoma', 'Truncorotalia crassula')

spp <- setdiff(sppAll, lostSpp)
rows2toss <- ! df$species %in% spp
df <- df[!rows2toss,]
row.names(df) <- as.character(1:nrow(df))

# Combine enviro and spp data ---------------------------------------------

envNm <- c('ann_temp_ym_dpth'
           #'month_temp_range',
           #'month_temp_max',
           #'month_temp_min',
           #'ann_otracer14_ym_dpth',
           #'ann_mixLyrDpth_ym_uo',
           #'ann_salinity_ym_dpth',
           #'ann_W_ym_dpth'
)
envNmShort <- sapply(envNm, function(txt){
  paste(strsplit(txt,'_')[[1]][2:3], collapse='_')
}) 

# Note: code below can deal with envNm that's a vector

llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Extract at each of 4 depths: 0, 40, 78, 164 m
# The lower 3 correspond to surface, surface-subsurface, and subsurface species
# and the values are based on the mean Average Living Depth of species in the dataset.
dpths <- c(1,4,6,8)
# Compare with the Average Living Depth of species in the dataset:
sppDat <- read.csv('Data/foram_spp_data_20-03-23.csv', stringsAsFactors=FALSE)
sameSpp <- sppDat$species %in% spp
sppDat <- sppDat[sameSpp,]
zones <- unique(sppDat$DepthHabitat)
for (z in zones){
  zRows <- which(sppDat$DepthHabitat==z)
  avg <- mean(sppDat$ALD[zRows])
  avg <- round(avg)
  print(paste(z, avg))
}
# [1] "Subsurface 164"
# [1] "Surface.subsurface 93"
# [1] "Surface 49"

source('raster_brick_import_fcn.R')

modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
bins <- unique(df$bin)
bins <- sort(bins)

addEnv <- function(bin, dat, mods, binCol, cellCol, prj, envNm){
  slcBool <- dat[,binCol] == bin
  slc <- dat[slcBool,]
  slcCells <- slc[,cellCol]
  slcEnv <- getBrik(bin=bin, envNm=envNm, mods=mods)
  
  for (i in 1:length(envNm)){
    env <- envNm[i]
    envVals <- raster::extract(slcEnv[[env]], slcCells)
    # Rows = points of extraction, columns = depth layers  
    envVals <- envVals[,dpths]
    nmOld <- envNmShort[i]
    nmNew <- paste(nmOld, c('0m','surf','surfsub','sub'), sep='_')
    slc[,nmNew] <- envVals
  } 
  slc
}

pkgs <- c('sp','raster') 
registerDoParallel(nCore)
sppEnv <- foreach(bin=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  addEnv(bin=bin, envNm=envNm, dat=df, mods=modId,
         binCol='bin', cellCol='cell_number', prj=llPrj)
stopImplicitCluster()

# Env at preferred depth --------------------------------------------------

# Species without at least 7 consecutive time bins can't be used in evo model analysis.
# These also lack depth data, so remove them now before calculating in-habitat temp.
keepSpp <- character()
binL <- bins[2] - bins[1]
enuf <- rep(binL, 7)
enufTxt <- paste0(enuf, collapse='')
spp <- unique(df$species)
for (s in spp){
  spBool <- df$species==s
  spDf <- df[spBool,]
  spB <- sort(unique(spDf$bin))
  bDiff <- diff(spB)
  diffTxt <- paste0(bDiff, collapse='')
  srch <- grep(enufTxt,diffTxt)
  if (length(srch) > 0){
    keepSpp <- c(keepSpp, s)
  }
}
keepBool <- df$species %in% keepSpp
df <- df[keepBool,]

habitatCol <- paste0(envNmShort, '_hab')
sppEnv[,habitatCol] <- NA
for (s in spp){
  sRow <- which(sppDat$species==s)
  habitat <- sppDat$DepthHabitat[sRow]
  if (is.na(habitat)) next else
  if (habitat=='Surface'){
    habNm <- paste0(envNmShort, '_surf')
  }
  if (habitat=='Surface.subsurface'){
    habNm <- paste0(envNmShort, '_surfsub')
  }
  if (habitat=='Subsurface'){
    habNm <- paste0(envNmShort, '_sub')
  }
  
  sBool <- sppEnv$species==s
  sppEnv[sBool, habitatCol] <- sppEnv[sBool, habNm]
}

# Clean -------------------------------------------------------------------

# remove records where environment is unknown
envCol <- c(habitatCol, paste0(envNmShort, '_0m'))
if (length(envCol)==1){
  naRows <- is.na(sppEnv[,envCol])
} else {
  naRows <- apply(sppEnv[,envCol], 1, function(r)
    any(is.na(r))
  )
}
sppEnv <- sppEnv[!naRows,]

allEnvCols <- c(paste(envNmShort, c('surf','surfsub','sub'), sep='_'), envCol)
df <- sppEnv[,c('species','bin','cell_number',
                'centroid_long','centroid_lat',
                'coreUniq',allEnvCols)]

# inspect correlation between temperature at surface (0m) and near preferred depth
cor(sppEnv$temp_ym_0m, sppEnv$temp_ym_hab)
# [1] 0.9557191

# Enviro at unique sampled sites ------------------------------------------

# get the sampling universe (env at range-through core sites) per bin
# coreAtts <- read.csv('Data/core_rangethrough_data_20-03-24.csv', stringsAsFactors=FALSE)
sampEnv <- function(b){
  # for each bin, find out which cores range through the 8ky interval
  inBin <- which(coreAtts$fad > (b-4) & coreAtts$lad < (b+4))
  cells <- coreAtts$cell[inBin]
  cells <- unique(cells)
  
  # get the env values at the core locations
  slcEnv <- getBrik(bin=b, envNm=envNm, mods=modId)
  for (i in 1:length(envNm)){
    env <- envNm[i]
    envVals <- raster::extract(slcEnv[[env]], cells)
    # Rows = points of extraction, columns = depth layers  
    envVals <- envVals[,dpths]
    nmOld <- envNmShort[i]
    nmNew <- paste(nmOld, c('0m','surf','surfsub','sub'), sep='_')
    colnames(envVals) <- nmNew
    sampled <- cbind(b, cells, envVals)
  }
  sampled
}
registerDoParallel(nCore)
samp <- foreach(b=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  sampEnv(b=b)
stopImplicitCluster()

# remove records where environment is unknown
naSampRows <- apply(samp[,-(1:2)], 1, function(r)
  any(is.na(r))
)
samp <- samp[!naSampRows,]

# Truncate to standard temp range -----------------------------------------

# Truncate series to the last 700 ka, encompassing 7 glacial/interglacial cycles
trimBool <- df$bin <= 700
df <- df[trimBool,]
sampTrim <- samp[,'b'] <= 700
samp <- samp[sampTrim,]
bins <- bins[bins <= 700]

if (doTrunc){
  
  minmax <- function(df, b, env){
    bBool <- df[,'b']==b
    slc <- df[bBool,]
    rng <- range(slc[,env])
    c(b, rng)
  }
  sampSmryM <- sapply(bins, minmax, df=samp, env='temp_ym_0m')
  sampSmry <- data.frame(t(sampSmryM))
  colnames(sampSmry) <- c('bin','min','max')
  uppr <- min(sampSmry$max)
  lwr <- max(sampSmry$min)
  
  p <- ggplot(data=sampSmry) + theme_bw() +
    scale_x_continuous(name='Time (ka)', expand=c(0.01,0)) +
    scale_y_continuous(name = 'MAT (degrees C)') +
    geom_linerange(aes(x=-bin, ymin=min, ymax=max), colour='red') +
    geom_linerange(aes(x=-bin, ymin=lwr, ymax=uppr), colour='black') +
    geom_hline(yintercept=uppr, colour='grey', lwd=1) +
    geom_hline(yintercept=lwr, colour='grey', lwd=1)
  
  pNm <- paste0('Figs/standardised_MAT_max_min_', day, '.pdf')
  pdf(pNm, width = 6, height=4)
  print(p)
  dev.off()
  
  trunc <- data.frame()
  for (b in bins){
    bBool <- df$bin==b
    slc <- df[bBool,]
    tooBig <- which(slc$temp_ym_0m > uppr)
    tooSmol <- which(slc$temp_ym_0m < lwr)
    out <- c(tooBig, tooSmol)
    if (length(out) > 0){
      slc <- slc[-out,]
    }
    trunc <- rbind(trunc, slc)
  }
  
  # apply truncation for sampled site data too
  truncSamp <- data.frame()
  for (b in bins){
    bBool <- samp[,'b']==b
    slc <- samp[bBool,]
    tooBig <- which(slc[,'temp_ym_0m'] > uppr)
    tooSmol <- which(slc[,'temp_ym_0m'] < lwr)
    out <- c(tooBig, tooSmol)
    if (length(out) > 0){
      slc <- slc[-out,]
    }
    truncSamp <- rbind(truncSamp, slc)
  }
  
  # * Evaluate degree of truncation -----------------------------------------
  
  df$trunc <- 'in range'
  tooBig <- which(df$temp_ym_0m > uppr)
  tooSmol <- which(df$temp_ym_0m < lwr)
  df$trunc[tooBig] <- 'high'
  df$trunc[tooSmol] <- 'low'
  
  mdrnBool <- df$bin %in% bins[1:2]
  mdrn <- df[mdrnBool,]
  old <- df[!mdrnBool,]
  old$trunc <- factor(old$trunc, levels=c('high','low','in range'))
  bars <- ggplot(data=old, aes(fill=trunc, x= - bin)) + 
    scale_x_continuous(name='time (excluding last 16 ka)', expand=c(0.01,0)) +
    scale_y_continuous(expand=c(0,0)) +
    #  theme_bw() +
    geom_bar(position="stack", width = 5) +
    scale_fill_manual(name='MAT value in relation to cutoffs', 
                      values=c('plum','gold','grey20')) +
    theme(legend.position = 'top')
  
  barNm <- paste0('Figs/truncated_data_sample_size_',day,'.pdf')
  pdf(barNm, width=6, height=4)
  print(bars)
  dev.off()
  
  # inspect the proportion of observations remaining
  nrow(trunc)/nrow(df) # all data
  table(old$trunc)['in range']/nrow(old) # excluding most recent 16 ka
  
  # end case where data are truncated to standard temperature range
} else {
  trunc <- df
  truncSamp <- samp
} 

# Clean -------------------------------------------------------------------

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

# Also re-check for per-species continuity through time 
# (at least 7 successive steps of 8ka)
keepSpp <- character()
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

# check for any gaps in the time series
binsObs <- sort(unique(trunc$bin))
if (any(diff(binsObs) != binL)) warning('discontinuous time bins')

if (doTrunc){
  obsNm <- paste0('Data/foram_MAT_occs_latlong_8ka_trunc_',day,'.csv')
  sampNm <- paste0('Data/samp_MAT_occs_latlong_8ka_trunc_',day,'.csv')
} else {
  obsNm <- paste0('Data/foram_MAT_occs_latlong_8ka_',day,'.csv')
  sampNm <- paste0('Data/samp_MAT_occs_latlong_8ka_',day,'.csv')
}
write.csv(trunc, obsNm, row.names = FALSE)
write.csv(truncSamp, sampNm, row.names = FALSE)
