library(tidyr)
library(lme4)

dat <- read.csv('Data/foram_niche_data_8ky_190912.csv',stringsAsFactors=FALSE)

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
getD <- function(s){
  spBool <- dat$species==s
  spCoords <- coords[spBool,]
  spSteps <- dat$bin[spBool]
  d <- dist(spCoords)
  # diagonal is all zeroes; remove top row to get dist btwn adjacent pts
  m <- as.matrix(d)[-1,]
  pairD <- c(diag(m),NA)
  
  # Don't calculate distance over gaps in time series
  gap <- diff(spSteps)!=8
  pairD[gap] <- NA
  
  cbind(species=s, bin=spSteps, d=pairD)
}
distL <- lapply(spp, getD)
distDf <- as.data.frame(do.call(rbind,distL))
distDf <- na.omit(distDf)

# Compare per-species and sampled niche space:
# for niche volume and niche centroid displacement

# Split out sampled niche space into a separate column 
# to be response variable in regression
distWide <- spread(distDf, key=species, value=d, convert=TRUE)
cols2merge <- setdiff(spp, 'sampled')
distLong <- gather(distWide, key=species, value=d, cols2merge)
distLong <- na.omit(distLong)
memD <- lmer(d ~ sampled + (1|species), data=distLong)
summary(memD)
confint(memD, parm='beta_', method='Wald')

hullDat <- dat 
hullDat$species <- as.factor(hullDat$species)
# 'spread' gives odd behaviour if there are more than the min # columns
hullDat <- hullDat[,c('bin','species','convex_hull_2axes')] # n_sites
volWide <- spread(hullDat, key=species, value=convex_hull_2axes, convert=TRUE) 
volLong <- gather(volWide, key=species, value=vol, cols2merge)
volLong <- na.omit(volLong)
memV <- lmer(vol ~ sampled + (1|species), data=volLong)
summary(memV)
confint(memV, parm='beta_', method='Wald')
