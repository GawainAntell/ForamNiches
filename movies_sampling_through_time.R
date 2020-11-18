library(raster)
library(ggplot2)

day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")
modCodes <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors = FALSE)
rAll <- list.files('Data/gcm_annual_mean/')
rTemp <- rAll[grep('temp', rAll)]

# standardized sampling universe (MAT at core sites) at each of 4 depths
dList <- readRDS('Data/spp-and-sampling-data_list-by-depth_2020-11-15.rds')
occ <- dList$temp_ym_0m$samp
bins <- unique(occ$bin)

# convert cell numbers to coordinates
rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
occ[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, occ$cell)

for (b in bins){
  # read in relevant raster file
  mod <- modCodes$age_1000ka == b
  id <- modCodes$id[mod]
  flPos <- grep(id, rTemp)
  fl <- rTemp[flPos]
  rNm <- paste0('Data/gcm_annual_mean/',fl)
  r <- raster(rNm)
  
  # plot points on top
  slc <- occ[occ$bin == b, ]
  coords <- slc[, c('centroid_long','centroid_lat')]
  x <- seq(-180, 180, length = 288)
  y <- seq(90, -90, length = 144)
  rDf <- expand.grid(x, y)
  rDf$val <- values(r)
  b3dig <- sprintf("%03d", b)
  lab <- paste(b3dig, 'ka')
  
  # set blue as low temp value
  colr <- rainbow(40)[c(30:1)] 
  
  rplot <- 
    ggplot() +
    ggtitle(lab) +
    scale_x_continuous(name = '', limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(name = '', limits = c(-90, 90),   expand = c(0,0)) +
    geom_raster(data = rDf, aes(x = Var1, y = Var2, fill = val), na.rm = TRUE) +
    geom_point(data = coords, aes(x = centroid_long, y = centroid_lat), size = 3) +
    scale_fill_gradientn(limits = c(0, 33), colors = colr) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(fill = 'Deg. C')
  # standardize colour scale across independent rasters
  
  pNm <- paste0('Figs/gif_series/map_surface_temp_and_occs_', b3dig, 'ka_', day, '.png')
  png(pNm, width = 1040, height = 696) 
  print(rplot)
  dev.off()
}

# cd C:\Users\sjoh4751\Dropbox\Gwen\ForamNiches\Git\Figs\gif_series
# magick convert -delay 50 -reverse -loop 1 map_surface_temp_and_occs*.png movie_S1_SS_temp_and_occs.gif
