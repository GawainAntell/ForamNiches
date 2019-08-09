library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(tidyr)
library(phytools)
library(paleoPhylo)
library(ape)

# save names to put packages on all cores later
pkgs <- c('sp','raster') 

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

# Read in occurrence and environmental data

occ_all <- read.csv('Data/foram_uniq_occs_latlong190808.csv', stringsAsFactors=FALSE)

# Annual files have 19 layers, but monthly max has only 1
r_all <- list.files('Data/', recursive = TRUE)
ann_pos <- intersect(grep('teini', r_all), grep('dpth.tif', r_all))
max_mo_pos <- intersect(grep('teini', r_all), grep('max', r_all))
ann_fls <- paste0('Data/', r_all[ann_pos])
max_mo_fl <- paste0('Data/', r_all[max_mo_pos])
r <- lapply(ann_fls, function(x){
  brik <- brick(x)
  brik[[1]]
})
r <- append(r, raster(max_mo_fl))
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Extract values at surface along each env axis, for each species

coords <- occ[,c('pal.long','pal.lat')]
pts <- SpatialPoints(sp_coords, proj4string = CRS(prj))

ncores <- detectCores() - 1
registerDoParallel(ncores) 
nich_m <- foreach(env=r, .packages=pkgs, .combine=rbind) %dopar% {
  vals <- extract(env, pts)
  env_trunc_name
  cbind(sp=occ$sp, axis=env_trunc_name, value=vals)
}
stopImplicitCluster()

# Re-format for PCA input: convert to short-form
nich_long <- data.frame(nich_m, stringsAsFactors=FALSE)
#colnames(nich_long) <- c('sp', 'axis', sums)
nich <- spread(nich_long, key=axis, value=value)

#############################
# PCA of niche position

# Species must have > 5 occs to reconstruct niche with adequate precision
spp <- unique(occ$species)
spp_freq <- table(occ$species)
too_scarce <- names( which(spp_freq<6) )
spp <- setdiff(spp, too_scarce)
n_spp <- length(spp)

# correlation between max mont and ann temp is VERY high (.89):
# exclude annual temp in PCA
cor(nich[,3:7])
ann_col <- which(colnames(nich)=='ann_temp_ym_dpth')
nich <- nich[,-ann_col]

row.names(nich) <- nich$sp
nich <- nich[,-1]
Pc <- prcomp(nich, center=TRUE, scale=TRUE)
  