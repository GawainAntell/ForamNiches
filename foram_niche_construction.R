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

source('Code/read_foram_data.R')

modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
envNms <- c('ann_otracer14_ym_dpth','ann_temp_ym_dpth',
             'ann_mixLyrDpth_ym_uo','month_temp_range'
            # 'ann_salinity_ym_dpth','ann_W_ym_dpth','max_month_temp'
             )
getR <- function(bin, envNms){
  modRow <- modId$age_1000ka == bin
  id <- modId$id[modRow]
  
  # Annual temp and water age files have 19 layers, 
  # but temp seasonality and mix layer depth have only 1
  rAll <- list.files('Data/', recursive = TRUE)
  txt <- paste0(id,'.*tif')
  rMod <- grep(txt, rAll)
  thickPos <- intersect(rMod, grep('otracer|ann_temp', rAll))
  flatPos <- intersect(rMod, grep('mixLyrDpth|temp_range', rAll))
  
  brickFls <- paste0('Data/', rAll[thickPos])
  flatFls <- paste0('Data/', rAll[flatPos])
  r1 <- lapply(brickFls, brick)
  r2 <- lapply(flatFls, raster)
  r <- append(r1, r2)
  names(r) <- envNms
  r
}
modrnR <- getR(bin=bins[1], envNms=envNms)
llPrj <- proj4string(modrnR[[1]])

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
  
  env_l <- get_r(bin=bin, envNms=envNms)
  
  nich <- matrix(nrow=0, ncol=5)
  for (sp in spp){
    sp_rows <- which(slc_dat[,sp_col]==sp)
    if (length(sp_rows) < 6){
      empt <- matrix(NA, nrow=length(env_l), ncol=5)
      empt[,1] <- sp
      empt[,2] <- envNms
      nich <- rbind(nich, empt)
    } else {
      sp_nich <- sapply(env_l, env_extract, sp=sp, dat=slc_dat, 
                        coord_cols=coord_cols, sp_col=sp_col, prj=prj
      )
      sp_nich <- cbind(sp, envNms, t(sp_nich))
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
colnames(avlbl_by_stat) <- c(envNms,'stat')

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
