library(ggplot2)
library(cowplot)
library(sp)
library(raster)
library(testthat)
library(kerneval)

day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")
df <- readRDS('Data/spp-and-sampling-data_list-by-depth_2020-07-21.rds')

plotDens <- function(bin, s, bw, xmn, xmx, sppDf, sampDf){
  # calculate w
  sampBool <- sampDf$bin == bin
  sampDat <- sampDf$temp_ym[sampBool] 
  # Feflecting the sample curve doesn't change the sp KDE much
  # except that the ends turn down a bit more (more convexity).
  # Since it's more complicated and throws warnings, don't do it.
  densSamp <- density(sampDat, bw = bw) # from=xmn, to=xmx 
  #  densSamp <- density.reflected(sampSlc, lower=xmn, upper=xmx) 
  w <- approxfun(densSamp$x, densSamp$y)
#  a <- min(densSamp$x)
#  b <- max(densSamp$x)
  
  sBool <- sppDf$species == s & sppDf$bin == bin
  spDat <- sppDf$temp_ym[sBool]
  # in order to plot the rug later, there must be >n discrete KDE points
  z1 <- tryCatch(
    transdens(spDat, w = w, bw = bw, reflect = FALSE, 
              a = xmn, b = xmx, n = 2^10),
    error = function(err){ list() }
  ) 
  if (length(z1)==0){ next }
  
#  pa <- max(z1$y)
#  pePos <- which.max(z1$y)
  
  x <- z1$x
  kd <- z1$y
  x <- c(xmn, x, xmx)
  kd <- c(0, kd, 0)
  plotDat <- data.frame(x = x, kd = kd)
  nlab <- paste0('n=', length(spDat))
  noNa <- length(x) - length(spDat)
  plotDat$rugHack <- c(spDat, rep(NA, noNa))
  sPlot <- 
    ggplot(data = plotDat, aes(x = x, y = kd)) +
    theme_bw() +
    ggtitle(paste(bin, 'ka')) +
    geom_area(fill = 'lightblue1') +
    geom_line() +
    # geom_segment(x = pe, xend = pe, y = 0, yend = 1,
    #              colour = 'blue') +
    # geom_segment(x = mean(spDat), xend = mean(spDat), 
    #              y = 0, yend = 1, colour = 'grey20') +
    geom_hline(yintercept = 0) +
    geom_rug(aes(x = rugHack), colour = 'grey20',
             sides = 'b', length = unit(0.045, "npc")) + 
    # ,outside = TRUE) +
    #   coord_cartesian(clip = "off") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 0.11)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 9))
  sPlot <- sPlot +
    geom_text(aes(fontface = 'plain', family = 'mono'), 
              label = nlab, size = 3.5, # 3.5mm = 10pt
              x = 0.15 * max(plotDat$x), y = 0.07)
  sPlot
}

# Figure 1 ----------------------------------------------------------------

modCodes <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors = FALSE)
rAll <- list.files('Data/gcm_annual_mean/')
rTemp <- rAll[grep('temp', rAll)]

# set blue as low temp value
colr <- rainbow(40)[c(30:1)]

b1 <- 28
b2 <- 124
samp0m <- df$temp_ym_0m$samp
spp0m <- df$temp_ym_0m$sp
pachy <- 'Neogloboquadrina pachyderma'

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
  
  # use N pachyderma as example to go with line drawing
  sPlot <- plotDens(bin = b, s = pachy,
           sampDf = samp0m, sppDf = spp0m,
           xmx=max(samp0m$temp_ym),
           xmn=min(samp0m$temp_ym)
           )
  sShrt <- strsplit(pachy, ' ')[[1]][2]
  spPlotNm <- paste0('Figs/KDE-example-', sShrt, '-', b, 'ka_', day, '.pdf')
  pdf(spPlotNm, width=2.5, height=2)
  print(sPlot)
  dev.off()
}

# Glacial/IG niches -------------------------------------------------------

coldBins <- c(28, 268, 436)
hotBins <- c(4, 124, 412)

# pick 1 surface and 1 subsurface species as examples
spp <- c('Neogloboquadrina dutertrei','Hirsutella scitula')

sppSurf <- df$temp_ym_surf$sp
sampSurf <- df$temp_ym_surf$samp
sppSub <- df$temp_ym_sub$sp
sampSub <- df$temp_ym_sub$samp

plotL <- list()  

for (s in spp){
  if (s=='Neogloboquadrina dutertrei'){
    sppDf <- sppSurf
    sampDf <- sampSurf
    xmn <- -1.4
    xmx <- 27.1
  } 
  if (s=='Hirsutella scitula'){
    sppDf <- sppSub
    sampDf <- sampSub
    xmn <- -0.8
    xmx <- 23.8
  }
  
for (bin in c(hotBins, coldBins)){
  
  sbPlot <- plotDens(bin, s, 'SJ-ste', xmn, xmx, sppDf, sampDf)
  plotL <- append(plotL, list(sbPlot))
  
}
}

# combine plots onto 1 page with title for each species
top <- plot_grid(plotlist = plotL[1:6], ncol = 3)
topTtl <- ggdraw() + 
  draw_label(
    'Neogloboquadrina dutertrei',
    fontface = 'bold.italic',
    x = 0, hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
btm <- plot_grid(plotlist = plotL[7:12], ncol = 3)
btmTtl <- ggdraw() +
  draw_label(
    'Hirsutella scitula',
    fontface = 'bold.italic',
    x = 0, hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

kdeNm <- paste0('Figs/KDE-examples-transform-SJste_',day,'.pdf')
pdf(kdeNm, width=6, height=6)
plot_grid(
  topTtl, top, btmTtl, btm,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)
dev.off()

# Evaluate degree of truncation -------------------------------------------

ss <- df$temp_ym_0m$sp
ss$trunc <- 'in range'
tooBig <- which(ss$temp_ym > 27.1)
tooSmol <- which(ss$temp_ym < -1.4)
ss$trunc[tooBig] <- 'high'
ss$trunc[tooSmol] <- 'low'
table(ss$trunc)/nrow(ss)
