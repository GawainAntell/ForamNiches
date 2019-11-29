library(ggplot2)
library(sp)
library(raster)

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

occ <- read.csv('Data/foram_uniq_occs_latlong_8ka_wEnv_190919.csv', stringsAsFactors=FALSE)
bins <- unique(occ$bin) #seq(4, 796, by=8)

# split data into present-day (last 16ky) vs. older
modrn_rows <- occ$bin==bins[1]
modrn <- occ[modrn_rows,]
old <- occ[!modrn_rows,]

# number of species (total or >5 occs) vs. time: line plot

countN <- function(bin,dat){
  slc_rows <- dat$bin==bin
  slc <- dat[slc_rows,]
  tbl <- table(slc$species)
  abund <- names(which(tbl > 5))
  n <- c(length(tbl), length(abund))
  names(n) <- c('n_tot','n_abund')
  return(n)
}

n_ser <- sapply(bins, countN, dat=occ)
n_df <- data.frame(t(n_ser))
n_df$t <- bins

lplot <- 
  ggplot(data=n_df) +
  theme_bw() + labs(title = '') +
  scale_x_continuous(name = 'ka (over 16 ky span)', expand=c(0,0),
                     limits=c(0, 800), breaks=seq(8,792,by=64)) +
  scale_y_continuous(name = 'n species') +
  geom_line(aes(x=t, y=n_tot, colour='All spp'), lwd=0.7) +
  geom_line(aes(x=t, y=n_abund, colour='> 5 occs'), lwd=0.7) +
  geom_point(aes(x=t, y=n_tot, colour='All spp'), size=.8) +
  geom_point(aes(x=t, y=n_abund, colour='> 5 occs'), size=.8) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 8))

lNm <- paste0('tseries_spp_count_', day, '.pdf')
pdf(lNm, width=5, height=3)
lplot
dev.off()


# relative abundance vs. time: stacked bar plot

old$bin <- as.factor(old$bin)
bplot <- 
  ggplot(data=old, aes(x=bin)) + theme_dark() +
  scale_x_discrete(name = 'ka (over 16 ky span)') +
  scale_y_continuous(name = 'n occurrences', 
                     expand=c(0,0), limits=c(0,650)) +
  geom_bar(aes(fill=species), colour='black', width = .7) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(5, 'mm')) +
  guides(fill=guide_legend(ncol=2))

bNm <- paste0('stacked_abundance_', day, '.pdf')
pdf(bNm, width=11, height=8)
bplot
dev.off()


# gif map of occurrences thorugh time

modCodes <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
r_all <- list.files('Data/gcm_annual_mean/')
r_temp <- r_all[grep('temp', r_all)]

for (b in bins){
  # read in relevant raster file
  mod <- modCodes$age_1000ka==b
  id <- modCodes$id[mod]
  fl_pos <- grep(id, r_temp)
  fl <- r_temp[fl_pos]
  rNm <- paste0('Data/gcm_annual_mean/',fl)
  r <- raster(rNm)
  
  # plot points on top
  slc <- occ[occ$bin==b,]
  dupes <- duplicated(slc$cell_number)
  slc <- slc[!dupes,]
  coords <- slc[,c('centroid_long','centroid_lat')]
  
  x <- seq(-180,180,length=288)
  y <- seq(90,-90,length=144)
  rDf <- expand.grid(x,y)
  rDf$val <- values(r)
  b3dig <- sprintf("%03d", b)
  lab <- paste(b3dig, 'ka')
  
  # set blue as low temp value
  colr <- rainbow(40)[c(30:1)] 
  
  rplot <- 
    ggplot() +
    ggtitle(lab) +
    scale_x_continuous(name='', limits=c(-180,180), expand=c(0,0)) +
    scale_y_continuous(name='', limits=c(-90,90), expand=c(0,0)) +
    geom_raster(data=rDf, aes(x=Var1, y=Var2, fill=val)) +
    geom_point(data=coords, aes(x=centroid_long, y=centroid_lat), size=3) +
    scale_fill_gradientn(limits=c(0,33), colors=colr) +
    theme(plot.title = element_text(hjust=.5, size=14),
          axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(fill='Deg. C')
    # standardize colour scale across independent rasters
  
  pNm <- paste0('Figs/gif_series/map_surface_temp_and_occs_', b3dig, 'ka_', day, '.png')
  png(pNm, width=1040, height=696) 
  print(rplot)
  dev.off()
}

# cd C:\Users\sjoh4751\Dropbox\Gwen\ForamNiches\Figs\gif_series
# magick convert -delay 50 -reverse -loop 1 map_surface_temp_and_occs*.png tseries_map_surface_temp_and_occs.gif


# gif map of BVF thorugh time

modCodes <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
r_all <- list.files('Data/BruntVaisala/')
r_tif <- r_all[grep('tif', r_all)] 

for (b in bins){
  # read in relevant raster file
  mod <- modCodes$age_1000ka==b
  id <- modCodes$id[mod]
  fl_pos <- grep(id, r_tif)
  fl <- r_tif[fl_pos]
  rNm <- paste0('Data/BruntVaisala/',fl)
  r <- raster(rNm)
  
  # plot points on top
  #slc <- occ[occ$bin==b,]
  #dupes <- duplicated(slc$cell_number)
  #slc <- slc[!dupes,]
  #coords <- slc[,c('centroid_long','centroid_lat')]
  
  x <- seq(-180,180,length=288)
  y <- seq(90,-90,length=144)
  rDf <- expand.grid(x,y)
  rDf$val <- values(r)
  b3dig <- sprintf("%03d", b)
  lab <- paste(b3dig, 'ka')
  
  # set blue as low temp value
  colr <- rainbow(40)[c(30:1)] 
  
  rplot <- 
    ggplot() +
    ggtitle(lab) +
    scale_x_continuous(name='', limits=c(-180,180), expand=c(0,0)) +
    scale_y_continuous(name='', limits=c(-90,90), expand=c(0,0)) +
    geom_raster(data=rDf, aes(x=Var1, y=Var2, fill=val)) +
    #geom_point(data=coords, aes(x=centroid_long, y=centroid_lat), size=3) +
    scale_fill_gradientn(limits=c(0,17), colors=colr) +
    theme(plot.title = element_text(hjust=.5, size=14),
          axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(fill='BVF')
  # standardize colour scale across independent rasters
  
  pNm <- paste0('Figs/gif_series/map_BVF360m_', b3dig, 'ka_', day, '.png')
  png(pNm, width=1040, height=696) 
  print(rplot)
  dev.off()
}

# In terminal/command prompt:
# cd C:\Users\sjoh4751\Dropbox\Gwen\ForamNiches\Git\Figs\gif_series
# magick convert -delay 50 -reverse -loop 1 map_BVF360m*.png tseries_map_BVF360m.gif
