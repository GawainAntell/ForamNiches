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
  
  # plot points on top of map
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
    scale_fill_gradientn(limits = c(-2, 34), colors = colr) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = 'left',
          legend.margin = margin(0,-15,0,0)) +
    labs(fill = 'Deg. C')
  # standardize colour scale across independent rasters
  
  pNm <- paste0('Figs/gif_series/', b3dig, 'ka_map.png')
  png(pNm, width = 6, height = 4.5, units = 'in', res = 300) 
  print(rplot)
  dev.off()
  
  # range of temperature values at global scale vs sampled sites
  # hard-coded alternative to boxplots
  boxs <- data.frame(
    x1 = c(1, 2.5),
    x2 = c(2, 3.5),
    y1 =  c(min(rDf$val,  na.rm = TRUE), min(slc$temp_ym)),
    y2 =  c(max(rDf$val,  na.rm = TRUE), max(slc$temp_ym)),
    avg = c(mean(rDf$val, na.rm = TRUE), mean(slc$temp_ym)),
    f = c('Global ', 'Sample sites')
  )
  p <- ggplot(data = boxs) +
    theme_classic() +
    geom_hline(aes(yintercept = c(-1.4, 27.1), 
                   #     color = 'Niche comparison interval'
    ),
    linetype = 'dashed',
    size = 1) +
    geom_rect(aes(xmin = x1, xmax = x2, 
                  ymin = y1, ymax = y2, 
                  fill = f),
              color = 'black') + # ,alpha = 0.75
    geom_segment(aes(x = x1, xend = x2,
                     y = avg, yend = avg),
                 size = 1) +
    scale_y_continuous('Range in mean annual temperature',
                       limits = c(-2, 34)) +
    scale_x_continuous('', limits = c(0.5, 4), expand = c(0, 0)) +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          #     axis.title.x = element_blank(),
          legend.position = 'top',
          legend.margin = margin(5,0,5,-25),
          legend.title = element_blank()) +
    scale_fill_manual(values = c('grey','white'), 
                      aesthetics = c('fill')) 
  #  scale_fill_manual(values = ('black'), aesthetics = 'color') +
  #  guides(fill = guide_legend(nrow = 2))
  
  bxNm <- paste0('Figs/gif_series/', b3dig, 'ka_boxes.png') # 'MAT-range_global-v-samp_'
  png(bxNm, width = 2, height = 4.5, units = 'in', res = 300) 
  print(p)
  dev.off()
}

# cd C:\Users\sjoh4751\Dropbox\Gwen\ForamNiches\Git\Figs\gif_series

# combine both figures at each time step:
# magick *_map.png -extent %[fx:s.w*4/3] null: *_boxes.png -gravity east -layers composite 2panel_%03d.png
# -reverse

# convert to gif or mp4 (smaller file):
# magick convert -delay 50 -reverse -loop 1 2panel_*.png movie_S1_SS_MAT.gif
# ffmpeg -r 2 -i 2panel_%03d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p output.mp4
