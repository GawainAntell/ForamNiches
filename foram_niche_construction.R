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
modrn <- occ_all$bin==8
occ <- occ_all[modrn,]

# Annual files have 19 layers, but monthly max has only 1
r_all <- list.files('Data/', recursive = TRUE)
ann_pos <- intersect(grep('teini', r_all), grep('dpth.tif', r_all))
max_mo_pos <- intersect(grep('teini', r_all), grep('max', r_all))
ann_fls <- paste0('Data/', r_all[ann_pos])
max_mo_fl <- paste0('Data/', r_all[max_mo_pos])
r <- lapply(ann_fls, brick)
r <- append(r, raster(max_mo_fl))
llPrj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Species must have > 5 occs to reconstruct niche with adequate precision
spp <- unique(occ$species)
spp_freq <- table(occ$species)
too_scarce <- names( which(spp_freq<6) )
spp <- setdiff(spp, too_scarce)
n_spp <- length(spp)

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
nich <- foreach(sp=spp, .packages=pkgs, .combine=rbind) %:% 
  foreach(env=r, .combine=rbind, .inorder=FALSE) %dopar% {
  env_extract(sp=sp, env=env,
    dat=occ, coord_cols = c('pal.long','pal.lat'),
    sp_col = 'species', prj = llPrj
  )
}
stopImplicitCluster()
# runtime is only seconds, on 8 core, 3.6 GHz processor

# Re-format for PCA input: convert to long-form, then new short-form
sums <- c('min', 'med', 'max')
nich <- data.frame(nich, stringsAsFactors=FALSE)
colnames(nich) <- c('sp', 'axis', sums)
nich[,sums] <- apply(nich[,sums], 2, as.numeric)
nich_long <- gather(nich, stat, value, sums, factor_key=TRUE)
nich_by_stat <- spread(nich_long, key=axis, value=value)

#############################
# (Phylo) PCA of niche position

# correlation between max mont and ann temp is VERY high (.89):
# exclude annual temp in PCA
med_rows <- nich_by_stat$stat=='med'
nich_med <- nich_by_stat[med_rows,]
cor(nich_med[,3:5])
ann_col <- which(colnames(nich_by_stat)=='ann_temp_ym_dpth')
nich_by_stat <- nich_by_stat[,-ann_col]

# format fully-birfurcating tree from Aze et al. 2011
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

# 8 species are not present in the phylogeny, mostly because microporiferate
lost_spp <- setdiff(spp, tr$tip.label)
rows2toss <- nich_by_stat$sp %in% lost_spp
nich_by_stat <- nich_by_stat[!rows2toss,]

# for (stat in sums){
stat <- 'max'

  sbset <- nich_by_stat$stat==stat
  nich_part <- nich_by_stat[sbset,]
  row_nms <- nich_part$sp
  nich_part <- nich_part[,-(1:2)]
  nich_part_scl <- apply(nich_part, 2, scale.default)
  row.names(nich_part_scl) <- row_nms
  
  pc <- prcomp(nich_part_scl, center=FALSE, scale=FALSE)
  
  # if unscaled, scale/units differ among axes, so use corr instead of cov
  # results differ slightly if scaled+cov vs. unscaled+corr
#  pc_phy <- phyl.pca(tr, nich_part_scl, method='BM', mode='cov')
#  pc_phy1 <- pc_phy$S[,'PC1']
# }

