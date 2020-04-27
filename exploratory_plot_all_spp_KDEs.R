library(ggplot2)
library(cowplot)
library(raster)

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

df <- readRDS('Data/sampled_temp_ym_truncated_by_depth_20-04-05.rds')

source('GSA_custom_ecospat_fcns.R')

coldBins <- c(28, 268, 436)
hotBins <- c(124, 332, 492)
spp <- c('Globigerinoides ruber','Hirsutella scitula')

# Figure 1 ----------------------------------------------------------------

modCodes <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors = FALSE)
rAll <- list.files('Data/gcm_annual_mean/')
rTemp <- rAll[grep('temp', rAll)]

# set blue as low temp value
colr <- rainbow(40)[c(30:1)]

b1 <- 28
b2 <- 124

for (b in c(b1, b2)){
  mod <- modCodes$age_1000ka==b
  id <- modCodes$id[mod]
  flPos <- grep(id, rTemp)
  fl <- rTemp[flPos]
  rNm <- paste0('Data/gcm_annual_mean/',fl)
  r <- raster(rNm)
  prj <- '+proj=laea' 
  m <- projectRaster(r, crs=prj)
  
  mapNm <- paste0('Figs/globe_map_',b,'ka.pdf')
  pdf(mapNm, width=20, height=20)
  plot(m, axes=FALSE, box=FALSE, legend=FALSE, col=colr)
  dev.off()
  
  # plot KDE of N pachyderma
}



# Niches for all species --------------------------------------------------

plotL <- list()

sppSurf <- df$temp_ym_surf$sp
sampSurf <- df$temp_ym_surf$samp
sppSub <- df$temp_ym_sub$sp
sampSub <- df$temp_ym_sub$samp
  
  #for (d in 1:2){
  #  sppDf <- list(sppSurf, sppSub)[[d]]
  #  sampDf <- list(sampSurf, sampSub)[[d]]
for (s in spp){
  if (s=='Globigerinoides ruber'){
    sppDf <- sppSurf
    sampDf <- sampSurf
  } 
  if (s=='Hirsutella scitula'){
    sppDf <- sppSub
    sampDf <- sampSub
  }
  xmx <- max(sampDf$temp_ym)
  xmn <- min(sampDf$temp_ym)
  
for (bin in c(hotBins, coldBins)){
  
  # calculate w
  sampBool <- sampDf$bin==bin
  sampDat <- sampDf$temp_ym[sampBool] 
  # Feflecting the sample curve doesn't change the sp KDE much
  # except that the ends turn down a bit more (more convexity).
  # Since it's more complicated and throws warnings, don't do it.
  densSamp <- density(sampDat, from=xmn, to=xmx) 
#  densSamp <- density.reflected(sampSlc, lower=xmn, upper=xmx) 
  w <- approxfun(densSamp$x, densSamp$y)
  a <- min(densSamp$x)
  b <- max(densSamp$x)
  
  sBool <- sppDf$species==s & sppDf$bin==bin
  spDat <- sppDf$temp_ym[sBool]
  # in order to plot the rug later, there must be >n discrete KDE points
  z1 <- tryCatch(
    transformEst(spDat, w = w, bw='nrd0', reflect = FALSE, 
                 a = a, b = b, n = 2^10),
    error = function(err){ list() }
  ) 
  if (length(z1)==0){ next }
  
  pe <- nichStats(z1)['pe']
  x <- z1$x
  kd <- z1$y
  x <- c(xmn, x, xmx)
  kd <- c(0, kd, 0)
  plotDat <- data.frame(x = x, kd = kd)
  nlab <- paste0('n=', length(spDat))
  noNa <- length(x) - length(spDat)
  plotDat$rugHack <- c(spDat, rep(NA, noNa))
  sPlot <- 
    ggplot(data=plotDat, aes(x=x, y=kd)) +
    theme_bw() +
    ggtitle(paste(bin, 'Ka')) +
    geom_area(fill='grey80') +
    geom_line() +
    geom_segment(x = pe, xend = pe, y = 0, yend = 1) +
    geom_segment(x = mean(spDat), xend = mean(spDat), 
                 y = 0, yend = 1, colour = 'blue') +
    geom_hline(yintercept = 0) +
    geom_rug(aes(x = rugHack), sides = 'b', length=unit(0.045, "npc")) + 
    # ,outside = TRUE) +
    #   coord_cartesian(clip = "off") +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(limits=c(0,0.08)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 9))
  sPlot <- sPlot +
    geom_text(aes(fontface = 1),label = nlab, size = 3,
              x = 3, y = 0.06) # 0.9 * max(plotDat$kd))
  
  plotL <- append(plotL, list(sPlot))
  
}
}

# print page of plots
pg <- plot_grid(plotlist=plotL, ncol=3)
kdeNm <- paste0('Figs/KDE-examples-transform-RT_',day,'.pdf')
pdf(kdeNm, width=6, height=6)
print(pg)
dev.off()

