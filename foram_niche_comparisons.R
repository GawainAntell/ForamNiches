library(tidyr)
library(lme4)
library(ecospat)
library(VoCC)
#devtools::install_github("JorGarMol/VoCC", dependencies = FALSE, build_vignettes = FALSE)

pcDat <- read.csv('Data/foram_uniq_occs_latlong_8ka_wEnv_190919.csv',stringsAsFactors=FALSE)
pcDat <- pcDat[,c('species','bin','cell_number','centroid_long','centroid_lat','pc1','pc2','pc3')]
bins <- unique(pcDat$bin)
spp <- unique(pcDat$species)
# for theyeri and adamsi, niche is quantified in only 1 time bin
# conglomeratea and dehiscens occur in only 1 set of consecutive bins
# humilis occurs in 2 sets of consecutive bins
# Globorotalia ungulata doesn't seem to have >5 occs in any bin
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

# Calculate niche overlap (Schoener's D) after Broennimann et al. 2012, PCA-occ method

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

# example plot
sp <- 'Globigerina bulloides'
b1 <- bins[1]
b2 <- bins[2]
glob <- pcDat[, c('pc1','pc2')]
glob1bool <- pcDat$bin==b1
glob1 <- glob[glob1bool,]
glob2bool <- pcDat$bin==b2
glob2 <- glob[glob2bool,]
sp1rows <- which(pcDat$species==sp & pcDat$bin==b1)
sp1 <- glob[sp1rows,]
sp2rows <- which(pcDat$species==sp & pcDat$bin==b2)
sp2 <- glob[sp2rows,]
z1 <- ecospat.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0.05, th.env=0.05)
z2 <- ecospat.grid.clim.dyn(glob, glob2, sp2, R, th.sp=0.05, th.env=0.05)

pNm <- paste0('Figs/G.bulloides_niche_heatmap', day, '.png')
png(pNm, width=1040, height=696) 
par(mfrow=c(1,2))
plot(z1$z, main=paste(bins[1],'ka'), xlab='PC1', ylab='PC2') # raster of niche in PC space
plot(z2$z, main=paste(bins[2],'ka'), xlab='PC1')
dev.off()

####################################################################
# Velocity of climate change

