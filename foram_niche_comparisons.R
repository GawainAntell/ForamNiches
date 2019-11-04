# TODO: calculate niche displacement on residuals of sp ~ sampled

library(tidyr)
library(lme4)

dat <- read.csv('Data/foram_niche_data_8ky_190912.csv',stringsAsFactors=FALSE)
bins <- unique(dat$bin)

# for theyeri and adamsi, niche is quantified in only 1 time bin
# conglomeratea and dehiscens occur in only 1 set of consecutive bins
# humilis occurs in 2 sets of consecutive bins
# Globorotalia ungulata doesn't seem to have >5 occs in any bin
noGood <- c('Hirsutella theyeri','Globigerinella adamsi',
            'Globoquadrina conglomerata','Sphaeroidinella dehiscens',
            'Globorotalia ungulata','Turborotalita humilis')
toss <- which(dat$species %in% noGood)
dat <- dat[-toss,]

spp <- unique(dat$species)

# calculate displacement for species through time
coords <- dat[,c('pc1_centroid','pc2_centroid','pc3_centroid')]
displace <- function(s){
  spBool <- dat$species==s
  spCoords <- coords[spBool,]
  spSteps <- dat$bin[spBool]
  d <- dist(spCoords)
  # diagonal is all zeroes; remove top row to get dist btwn adjacent pts
  m <- as.matrix(d)[-1,]
  pairDisp <- c(diag(m),NA)
  
  # Don't calculate distance over gaps in time series
  gap <- diff(spSteps)!=8
  pairDisp[gap] <- NA
  
  cbind(species=s, bin=spSteps, d=pairDisp)
}
dispL <- lapply(spp, displace)
dispDf <- as.data.frame(do.call(rbind,dispL))
dispDf <- na.omit(dispDf)

# Compare per-species and sampled niche space:
# for niche volume and niche centroid displacement

# Split out sampled niche space into a separate column 
# to be response variable in regression
dispWide <- spread(dispDf, key=species, value=d, convert=TRUE)
cols2merge <- setdiff(spp, 'sampled')
dispLong <- gather(dispWide, key=species, value=d, cols2merge)
dispLong <- na.omit(dispLong)
memDisp <- lmer(d ~ sampled + (1|species), data=dispLong)
summary(memDisp)
confint(memDisp, parm='beta_', method='Wald')

# save lm residuals for analysis
# (not MEM residuals bc interspecific differences should be kept)
#lmDisp <- lm(d ~ sampled, data=dispLong)
#dispLong$dResid <- lmDisp$residuals

####################################################################
# Calculate niche overlap (Schoener's D) after Broennimann et al. 2012
# PCA-occ method

library(ecospat)

pcDat <- read.csv('Data/foram_uniq_occs_latlong_8ka_wEnv_190919.csv',stringsAsFactors=FALSE)
pcDat <- pcDat[,c('species','bin','cell_number','centroid_long','centroid_lat','pc1','pc2','pc3')]

# if skipping to this section without running the sections above:
bins <- unique(pcDat$bin)
spp <- unique(pcDat$species)
noGood <- c('Hirsutella theyeri','Globigerinella adamsi',
            'Globoquadrina conglomerata','Sphaeroidinella dehiscens',
            'Globorotalia ungulata','Turborotalita humilis')
toss <- which(spp %in% noGood)
spp <- spp[-toss]

# Subset occurrences such that eas sp has >5 occs per bin
saveRows <- function(sp, bin){
  spRows <- which(pcDat$species==sp & pcDat$bin==bin)
  if (length(spRows)>5){
    spRows
  }
}  
keepRowsL <- sapply(spp, function(x){
  sapply(bins, function(b)saveRows(sp=x, bin=b))
} 
)
keepRows <- unlist(keepRowsL)
pcDat <- pcDat[keepRows,]
  
# Resolution of the gridding of the climate space
R <- 100

overlapD <- function(b1, b2, sp){
  glob <- pcDat[, c('pc1','pc2')]
  glob1bool <- pcDat$bin==b1
  glob1 <- glob[glob1bool,]
  glob2bool <- pcDat$bin==b2
  glob2 <- glob[glob2bool,]
  sp1rows <- which(pcDat$species==sp & pcDat$bin==b1)
  sp1 <- glob[sp1rows,]
  sp2rows <- which(pcDat$species==sp & pcDat$bin==b2)
  sp2 <- glob[sp2rows,]
  
  # for each species at time i and i+1
  z1 <- tryCatch(
    ecospat.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0.05, th.env=0.05),
    error = function(err){ NA }
  ) 
  z2 <- tryCatch(
    ecospat.grid.clim.dyn(glob, glob2, sp2, R, th.sp=0.05, th.env=0.05),
    error = function(err){ NA }
  ) 
  # plot(z1$z) # raster of niche in PC space
  
  if (any(is.na(list(z1, z2)))){
    NA
  } else {
    ovrlp <- ecospat.niche.overlap(z1, z2, cor=FALSE)
    ovrlp$D
  }
}

bPairs <- cbind(bins[-1], bins[-length(bins)])
Dlist <- lapply(spp, function(s){
  apply(bPairs, 1, function(x){
    overlapD(b1=x[1], b2=x[2], sp=s)
  })
})
dDf <- do.call(cbind, Dlist)
row.names(dDf) <- bins[-1]
colnames(dDf) <- spp

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')
dfNm <- paste0('Data/foram_niche_overlap_D_',day,'.csv')
write.csv(dDf, dfNm)

####################################################################
# Velocity of climate change

