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

modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
env_nms <- c('ann_otracer14_ym_dpth','ann_salinity_ym_dpth','ann_temp_ym_dpth',
             'ann_W_ym_dpth','max_month_temp')
get_r <- function(bin, env_nms){
  mod_row <- modId$age_1000ka == bin
  id <- modId$id[mod_row]
  
  # Annual files have 19 layers, but monthly max has only 1
  r_all <- list.files('Data/', recursive = TRUE)
  ann_pos <- intersect(grep(id, r_all), grep('dpth.tif', r_all))
  max_mo_pos <- intersect(grep(id, r_all), grep('max', r_all))
  ann_fls <- paste0('Data/', r_all[ann_pos])
  max_mo_fl <- paste0('Data/', r_all[max_mo_pos])
  r <- lapply(ann_fls, brick)
  r <- append(r, raster(max_mo_fl))
  names(r) <- env_nms
  r
}
modrn_r <- get_r(bin=8, env_nms=env_nms)
llPrj <- proj4string(modrn_r[[1]])

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

occ <- read.csv('Data/foram_uniq_occs_latlong190808.csv', stringsAsFactors=FALSE)
bins <- unique(occ$bin)

# 8 species are not present in the phylogeny, mostly because microporiferate
spp_all <- unique(occ$species)
lost_spp <- setdiff(spp_all, tr$tip.label)

# Species must have > 5 occs in a bin to reconstruct niche with adequate precision
spp <- vector()
for (b in bins){
  slc_bool <- occ$bin==b
  slc_spp <- occ$species[slc_bool]
  slc_freq <- table(slc_spp)
  slc_ok <- names( which(slc_freq > 5) )
  spp <- union(spp, slc_ok)
}

spp <- setdiff(spp, lost_spp)
rows2toss <- ! occ$species %in% spp
occ <- occ[!rows2toss,]
row.names(occ) <- as.character(1:nrow(occ))

#############################
# Extract min, med, max values for each env axis/species/bin combo

env_extract <- function(sp, env, dat, coord_cols, sp_col, prj){
  sp_bool <- dat[,sp_col] == sp
  sp_coords <- dat[sp_bool, coord_cols]
  sp_pts <- SpatialPoints(sp_coords, proj4string = CRS(prj))
  
  # If extracted from Brick, values are put in a matrix
  # Rows = points of extraction, columns = depth layers
  vals <- raster::extract(env, sp_pts)
  mn <- min(vals, na.rm = TRUE)
  med <- median(vals, na.rm = TRUE)
  mx <- max(vals, na.rm = TRUE)
  c(mn, med, mx)
}

bin_nich <- function(spp=NULL, bin, dat, bin_col, coord_cols, sp_col, prj, env_names){
  if (is.null(spp)){
    spp <- unique(dat[,sp_col])
  }
  
  slc_bool <- dat[,bin_col] == bin
  slc_dat <- dat[slc_bool,]
  
  env_l <- get_r(bin=bin, env_nms=env_nms)
  
  nich <- matrix(nrow=0, ncol=5)
  for (sp in spp){
    sp_rows <- which(slc_dat[,sp_col]==sp)
    if (length(sp_rows) < 6){
      empt <- matrix(NA, nrow=length(env_l), ncol=5)
      empt[,1] <- sp
      empt[,2] <- env_nms
      nich <- rbind(nich, empt)
    } else {
      sp_nich <- sapply(env_l, env_extract, sp=sp, dat=slc_dat, 
                        coord_cols=coord_cols, sp_col=sp_col, prj=prj
      )
      sp_nich <- cbind(sp, env_nms, t(sp_nich))
      nich <- rbind(nich,  sp_nich)
    }
  }
  
  nich <- cbind(nich, bin)
  nich
}

# Calculate 'available' environment from all occs in every bin
# Calculated the same as for each species, so treat as if a species
# Note: some species may be at same cell in a bin
get_dupes <- function(b){
  slc_bool <- occ$bin==b
  slc <- occ[slc_bool,]
  dupes <- duplicated(slc$cell_number)
  row.names(slc[dupes,])
}
dupe_rows <- unlist(sapply(bins, get_dupes))
avlbl <- occ[-as.numeric(dupe_rows),]
avlbl$species <- 'available'
occ <- rbind(occ, avlbl)

ncores <- detectCores() - 1
pt1 <- proc.time()
registerDoParallel(ncores)

nich <- foreach(bin=bins, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
bin_nich(bin=bin, env_names=env_names, dat=occ, bin_col='bin', 
         coord_cols=c('pal.long','pal.lat'), sp_col='species', prj=llPrj
         )

stopImplicitCluster()
pt2 <- proc.time()
(pt2-pt1)/60

# Re-format for PCA input: convert to long-form, then new short-form
sums <- c('min', 'med', 'max')
nich <- data.frame(nich, stringsAsFactors=FALSE)
colnames(nich) <- c('sp', 'axis', sums, 'bin')
row.names(nich) <- as.character(1:nrow(nich))
nich[,sums] <- apply(nich[,sums], 2, as.numeric)
nich_long <- gather(nich, stat, value, sums, factor_key=TRUE)

nich_by_stat <- spread(nich_long, key=axis, value=value)
# write.csv(nich_by_stat, paste0('Data/nich_unscaled_all_axes_',day,'.csv'), row.names=FALSE)


#############################
# PCA of niche position

# Retrieve eigenvalues for PC rotation of available environment.

modrn_rows <- avlbl$bin==8
avlbl_coords <- avlbl[modrn_rows, c('pal.long','pal.lat')]
avlbl_pts <- SpatialPoints(avlbl_coords, proj4string = CRS(llPrj))
avlbl_vals <- lapply(modrn_r, raster::extract, y=avlbl_pts)
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
ann_t_col <- which(colnames(nich_by_stat)=='ann_temp_ym_dpth')
nich_by_stat <- nich_by_stat[,-ann_t_col]

# Temperature can be negative but other vars shouldn't be
negtv_t <- min(nich_by_stat$max_month_temp, na.rm=TRUE)
nich_by_stat$max_month_temp <- nich_by_stat$max_month_temp - negtv_t
no_salt <- which(nich_by_stat$ann_salinity_ym_dpth < 0)
nich_by_stat$ann_salinity_ym_dpth[no_salt] <- 0
too_slow <- which(nich_by_stat$ann_W_ym_dpth < 0)
nich_by_stat$ann_W_ym_dpth[too_slow] <- 0

nich_part <- nich_by_stat[,-(1:3)]
nich_part <- na.omit(nich_part)
pc <- prcomp(nich_part, center=TRUE, scale=TRUE) 
pc_rot <- pc_avlbl$rotation
rot_nm <- paste0('Data/modern_available_env_PC_rotation_',day,'.csv')
write.csv(pc_rot, rot_nm)

nich_pc <- stats::predict(pc, nich_by_stat)
fin <- cbind(nich_by_stat[,1:3], nich_pc)
write.csv(fin, paste0('Data/niche_PC_4axes',day,'.csv'), row.names=FALSE)
