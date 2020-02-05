library(ggplot2)
library(sp)
library(raster)
library(cowplot)
library(grid)
library(gridExtra)

# Data prep ---------------------------------------------------------------

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

occ <- read.csv('Data/foram_MAT_occs_latlong_8ka_200129.csv', stringsAsFactors=FALSE)
bins <- unique(occ$bin) 

# split data into present-day (last 16 ka) vs. older
modrn_rows <- occ$bin %in% bins[1:2]
modrn <- occ[modrn_rows,]
old <- occ[!modrn_rows,]

# Sample size plots -------------------------------------------------------

# * per-species sample size -----------------------------------------------
# boxplot of n occs per species, at each time bin

sppN <- function(df,b){
  bBool <- df$bin==b
  bDf <- df[bBool,]
  nmTbl <- table(bDf$species)
  spp <- names(nmTbl)
  n <- as.numeric(nmTbl)
  cbind(n=n, t=b)
}
nDat <- sapply(bins[-(1:2)], sppN, df=old)
nDf <- data.frame(do.call(rbind, nDat))
nDf$t <- factor(-nDf$t)

xlab <- rev(bins[-(1:2)])
bp <- ggplot(data=nDf, aes(x=t, y=n)) +
  geom_hline(yintercept=10, colour='blue')+
  geom_boxplot(width=0.65) +
  scale_y_continuous(name='occurrences per species') +
  scale_x_discrete(name='ka', labels=xlab) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=8)
        #      panel.grid.major = element_line(size = 0.5, linetype = 'solid', 
        #                                      colour = "white")
  )

bpNm <- paste0('Figs/n_per_sp_8ky_', day, '.pdf')
pdf(bpNm, width=9, height=5)
bp
dev.off()

# * n species total -------------------------------------------------------

countN <- function(bin,dat){
  slcBool <- dat$bin==bin
  spp <- dat$species[slcBool]
  length(unique(spp))
}

# don't count the 'species' that is from sampled locations
nSer <- sapply(bins, countN, dat=occ) - 1
nDf <- data.frame(cbind(n=nSer, t=bins))

lplot <- 
  ggplot(data=nDf) +
  theme_bw() + labs(title = '') +
  scale_x_continuous(name = 'ka (8ky-long bins)', expand=c(0,0)) +
  scale_y_continuous(name = 'n species', limits=c(0,28), expand=c(0,0)) +
  geom_line(aes(x=-t, y=n), lwd=0.7) +
  geom_point(aes(x=-t, y=n), size=.8) 

lNm <- paste0('tseries_spp_count_', day, '.pdf')
pdf(lNm, width=4, height=3)
lplot
dev.off()


# Example KDE plots -------------------------------------------------------

source('GSA_custom_ecospat_fcns.R')
nbins <- length(bins)
h.method <- "nrd0" # "SJ-ste" # "ucv"
R <- 2^8

env <- 'temp_ym_0m'
spp <- unique(occ$species)
xmax <- max(occ[,env])
xmin <- min(occ[,env])
focalB <- seq(100, 700, by=200)
# b <- 100; s <- spp[2]

for (b in focalB){
  bBool <- occ$bin==b
  slc <- occ[bBool,]
  
  plotL <- list()
  for (s in spp){
    spRows <- which(slc$species==s)
    spSlc <- slc[spRows,env]
    z1 <- tryCatch(
      kdeNiche(sp=spSlc, xmax=xmax, xmin=xmin, R=R, h.method=h.method),
      error = function(err){ list() }
    ) 
    if (length(z1)==0){next}
    
    x <- z1$x
    kd <- z1$z
    x <- c(xmin, x, xmax)
    kd <- c(0, kd, 0)
    plotDat <- data.frame(x=x, kd=kd)
    nlab <- paste0('n=', length(spSlc))
    noNa <- length(x) - length(spSlc)
    plotDat$rugHack <- c(spSlc, rep(NA, noNa))
    sPlot <- 
      ggplot(data=plotDat, aes(x=x, y=kd)) +
        theme_bw() +
        ggtitle(s) +
        geom_area(fill='grey80') +
        geom_line() +
        geom_segment(x=z1$pe, xend=z1$pe, y=0, yend=1) +
        geom_segment(x=mean(spSlc), xend=mean(spSlc), y=0, yend=1, colour='blue') +
        geom_hline(yintercept=0) +
        geom_rug(aes(x=rugHack), sides='b', length=unit(0.045, "npc")) + # ,outside = TRUE
     #   coord_cartesian(clip = "off") +
        scale_x_continuous(expand=c(0,0)) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(size=9))
    sPlot <- sPlot +
      geom_text(aes(fontface=1),label=nlab, size=3,
                x=3, y=0.9*max(plotDat$kd))
      
    plotL <- append(plotL, list(sPlot))
  }
  
  # print page of plots
  pg <- plot_grid(plotlist=plotL, nrow=6)
  kdeNm <- paste0('Figs/KDE_all_spp_',h.method,'_',b,'ka_',day,'.pdf')
  pdf(kdeNm, width=8.5, height=11)
  grid.arrange(arrangeGrob(pg))
  dev.off()
}

# Map species by depth habitat --------------------------------------------

# maps of all species, last 16 ka
# label and sort species by depth habitat
# surface layer temperature raster

mod <- modCodes$age_1000ka==8
id <- modCodes$id[mod]
flPos <- grep(id, rTemp)
fl <- rTemp[flPos]
rNm <- paste0('Data/gcm_annual_mean/',fl)
# surface-level temperature
r <- brick(rNm)[[1]]
x <- seq(-180,180,length=288)
y <- seq(90,-90,length=144)
rDf <- expand.grid(x,y)
rDf$val <- values(r)

# set blue as low temp value
colr <- rainbow(40)[c(30:1)] 

plotSp <- function(s){
  sBool <- modrn$species==s
  slc <- modrn[sBool,]
  coords <- slc[,c('centroid_long','centroid_lat')]
  
  ggplot() +
    ggtitle(s) +
    scale_x_continuous(name='', limits=c(-180,180), expand=c(0,0)) +
    scale_y_continuous(name='', limits=c(-90,90), expand=c(0,0)) +
    geom_raster(data=rDf, aes(x=Var1, y=Var2, fill=val)) +
    geom_point(data=coords, aes(x=centroid_long, y=centroid_lat), size=1) +
    scale_fill_gradientn(limits=c(0,33), colors=colr) +
    theme(plot.title = element_text(hjust=.5, size=8),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = 'none') 
}

spp <- unique(modrn$species)

# combine coordinates with species attribute data (depth habitat)
atts <- read.csv('Data/foram_spp_data_200108.csv', stringsAsFactors=FALSE)
#plotDf <- merge(modrn, atts, by.x='species')
#plotDf$DepthHabitat <- factor(plotDf$DepthHabitat, 
#                              levels=c('Subsurface','Surface.subsurface','Surface'))
deepRows <- which(atts$DepthHabitat=='Subsurface')
deep <- atts$species[deepRows]
deep <- intersect(deep, spp)
midRows <- which(atts$DepthHabitat=='Surface.subsurface')
mid <- atts$species[midRows]
mid <- intersect(mid, spp)
shalRows <- which(atts$DepthHabitat=='Surface')
shal <- atts$species[shalRows]
shal <- intersect(shal, spp)

deepPlots <- lapply(deep, plotSp)
midPlots <- lapply(mid, plotSp)
shalPlots <- lapply(shal, plotSp)

deepMulti <- plot_grid(plotlist=deepPlots, ncol=5, align='v') 
midMulti <- plot_grid(plotlist=midPlots, ncol=5, align='v') 
shalMulti <- plot_grid(plotlist=shalPlots, ncol=5, align='v')

deepPanel <- grid.arrange(arrangeGrob(deepMulti), top='Subsurface species, 160 m') 
midPanel <- grid.arrange(arrangeGrob(midMulti), top='Surface-Subsurface species, 80 m') 
shalPanel <- grid.arrange(arrangeGrob(shalMulti), top='Surface species, 40 m') 

spMapsNm <- paste0('Figs/species_maps_8ka_by_depth_',day,'.pdf')
pdf(spMapsNm, width=8.5, height=11)
grid.arrange(shalPanel, midPanel, deepPanel, 
             heights=c(1,1,0.5),
             layout_matrix=rbind(1,2,3))
dev.off()

# gifs --------------------------------------------------------------------

# * occurrences thorugh time ----------------------------------------------

modCodes <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
rAll <- list.files('Data/gcm_annual_mean/')
rTemp <- rAll[grep('temp', rAll)]

for (b in bins){
  # read in relevant raster file
  mod <- modCodes$age_1000ka==b
  id <- modCodes$id[mod]
  flPos <- grep(id, rTemp)
  fl <- rTemp[flPos]
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
  
  # also export snapshots in Mollweide projection for a few time steps
  snaps <- c(0, 200, 400, 600) + bins[1]
  ylims <- c(0,30)
  mProj <- '+proj=laea'
    #"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84"
  if (b %in% snaps){
    m <- projectRaster(r, crs=mProj)
    mNm <- paste0('Figs/Mollweide_surface_temp_', b3dig, 'ka_', day, '.png')
    png(mNm, width=2080, height=2080)
    plot(m, axes=FALSE, box=FALSE, legend=FALSE, col=colr)
    dev.off()
  }
}

# cd C:\Users\sjoh4751\Dropbox\Gwen\ForamNiches\Git\Figs\gif_series
# magick convert -delay 50 -reverse -loop 1 map_surface_temp_and_occs*.png tseries_map_surface_temp_and_occs.gif

# * BVF thorugh time ------------------------------------------------------

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
