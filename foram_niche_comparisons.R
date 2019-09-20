# TODO: calculate niche displacement on residuals of sp ~ sampled

library(tidyr)
library(lme4)

dat <- read.csv('Data/foram_niche_data_8ky_190912.csv',stringsAsFactors=FALSE)
bins <- unique(dat$bin)

# for theyeri and adamsi, niche is quantified in only 1 time bin
# conglomeratea and dehiscens occur in only 1 set of consecutive bins
# humilis occurs in 2 sets of consecutive bins
noGood <- c('Hirsutella theyeri','Globigerinella adamsi',
            'Globoquadrina conglomerata','Sphaeroidinella dehiscens',
            'Turborotalita humilis')
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

# Resolution of the gridding of the climate space
R <- 100

glob <- pcDat[, c('pc1','pc2')]
glob1bool <- pcDat$bin==bins[1]
glob1 <- pcDat[glob1bool,c('pc1','pc2')]
glob2bool <- pcDat$bin==bins[2]
glob2 <- pcDat[glob2bool,c('pc1','pc2')]
sp1rows <- which(pcDat$species==spp[2] & pcDat$bin==bins[1])
sp1 <- pcDat[sp1rows, c('pc1','pc2')]
sp2rows <- which(pcDat$species==spp[2] & pcDat$bin==bins[2])
sp2 <- pcDat[sp2rows, c('pc1','pc2')]

# for each species at time i and i+1
z <- ecospat.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0.05, th.env=0.05)
# plot(z$z) # raster of niche in PC space
z2 <- ecospat.grid.clim.dyn(glob, glob2, sp2, R, th.sp=0.05, th.env=0.05)
plot(z$z)
plot(z2$z)

ovrlp <- ecospat.niche.overlap(z, z2, cor=FALSE)
ovrlp$D
