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


#############################
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
env_nms <- vector()
for (fl in c(ann_fls,max_mo_fl)){
  env_nm_trunc <- tail(strsplit(fl, 'mean', fixed=TRUE)[[1]], 1)
  nm_lngth <- nchar(env_nm_trunc)
  env_nm_trunc <- substr(env_nm_trunc, 9, nm_lngth-4)
  env_nms <- c(env_nms, env_nm_trunc)
}
names(r) <- env_nms
llPrj <- proj4string(r[[1]])

# Limit analysis to species included in tree from Aze et al. 2011
tr_raw <- read.csv('Data/Aze_et_al_2011_bifurcating_tree_data.csv', stringsAsFactors=FALSE)
tr_raw$Species.name <- gsub('Globigerinoides sacculifer', 'Trilobatus sacculifer', tr_raw$Species.name)
tr_raw$Species.name <- gsub('Globigerinoides trilobus', 'Trilobatus trilobus', tr_raw$Species.name)
tr_paleo <- with(tr_raw, 
                 as.paleoPhylo(Species.code, Ancestor.code, Start.date, End.date)
)
tr <- buildApe(tr_paleo)
spp_codes <- tr$tip.label
row_ordr <- match(spp_codes, tr_raw$Species.code)
tr$tip.label <- tr_raw$Species.name[row_ordr]

# Species must have > 5 occs to reconstruct niche with adequate precision
spp <- unique(occ$species)
spp_freq <- table(occ$species)
too_scarce <- names( which(spp_freq<6) )

# 8 species are not present in the phylogeny, mostly because microporiferate
lost_spp <- setdiff(spp, tr$tip.label)
rows2toss <- occ$species %in% union(lost_spp, too_scarce)
occ <- occ[!rows2toss,]

# Calculate 'available' environment from all occs in every bin
# Calculated the same as for each species, so treat as if a species
avlbl <- occ
avlbl$species <- 'available'
occ <- rbind(occ, avlbl)

spp <- unique(occ$species)
n_spp <- length(spp)

#############################
# Retrieve eigenvalues for PC rotation of available environment

avlbl_coords <- avlbl[,c('pal.long','pal.lat')]
avlbl_pts <- SpatialPoints(avlbl_coords, proj4string = CRS(llPrj))
avlbl_vals <- lapply(r, raster::extract, y=avlbl_pts)
sum_vals <- function(x, fn, na.rm=TRUE){
  if (class(x)=='matrix'){
    out <- apply(x, 1, fn, na.rm=na.rm)
  } 
  if (class(x)=='numeric'){
    out <- x
  } 
  return(out)
}
avlbl_cln <- lapply(1:3, function(i){
  f <- list(min, median, max)[[i]]
  avlbl_l <- lapply(avlbl_vals, sum_vals, fn=f)
  avlbl_m <- do.call(cbind, avlbl_l)
  
  # some points lack env data; value will be -Inf
  avlbl_m[is.infinite(avlbl_m)] <- NA
  avlbl_df <- data.frame(na.omit(avlbl_m))
  avlbl_df$stat <- c('min','med','max')[i]
  avlbl_df
}
)
avlbl_by_stat <- do.call(rbind, avlbl_cln)
colnames(avlbl_by_stat) <- c(env_nms,'stat')

# mean ann temp and max monthly temp are too collinear (cor=.8)
med_rows <- avlbl_by_stat$stat=='med'
avlbl_med <- avlbl_by_stat[med_rows,-6]
cor(avlbl_med)
ann_t_col <- which(env_nms=='ann_temp_ym_dpth')
avlbl_med <- avlbl_med[,-ann_t_col]
pc_avlbl <- prcomp(avlbl_med, center=TRUE, scale=TRUE)
pc_rot <- pc_avlbl$rotation


#############################
# Extract min, med, max values along each env axis, for each species

env_extract <- function(sp, env, dat, coord_cols, sp_col, prj){
  sp_bool <- dat[,sp_col] == sp
  sp_coords <- dat[sp_bool, coord_cols]
  sp_pts <- SpatialPoints(sp_coords, proj4string = CRS(prj))
  
  # If extracted from Brick, values are put in a matrix
  # Rows = points of extraction, columns = layers
  vals <- extract(env, sp_pts)
  mn <- min(vals, na.rm = TRUE)
  med <- median(vals, na.rm = TRUE)
  mx <- max(vals, na.rm = TRUE)
  
  # Extract the name of the environmental variable
  env_nm <- env@file@name
  env_nm_trunc <- tail(strsplit(env_nm, '\\', fixed=TRUE)[[1]], 1)
  nm_lngth <- nchar(env_nm_trunc)
  env_nm_trunc <- substr(env_nm_trunc, 8, nm_lngth-4)
  
  c(sp, env_nm_trunc, mn, med, mx)
}

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

nich_part <- nich_by_stat[,-(1:2)]
nich_part_scl <- apply(nich_part, 2, scale.default)
nich_part_pc <- nich_part_scl*pc_rot
nich_part_pc <- cbind(nich_by_stat[,1:2], nich_part_pc)

