library(ggplot2)
library(cowplot)
library(raster)

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

df <- read.csv('Data/foram_uniq_occs_latlong_8ka_20-04-05.csv',
               stringsAsFactors = FALSE)
samp <- read.csv('Data/samp_uniq_occs_latlong_8ka_20-04-05.csv',
                 stringsAsFactors = FALSE)

source('GSA_custom_ecospat_fcns.R')
bins <- unique(df$bin)
nbins <- length(bins)
spp <- unique(df$species)
env <- 'temp_ym_0m'
xmx <- max(df[,env])
xmn <- min(df[,env])
focalB <- seq(100, 700, by=200)
# b <- 100; s <- spp[2]

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

for (bin in focalB){
  bBool <- df$bin==bin
  slc <- df[bBool,]
  
  # calculate w
  sampBool <- samp$b==bin
  sampSlc <- samp[sampBool,env]
  # Feflecting the sample curve doesn't change the sp KDE much
  # except that the ends turn down a bit more (more convexity).
  # Since it's more complicated and throws warnings, don't do it.
  densSamp <- density(sampSlc, from=xmn, to=xmx) 
#  densSamp <- density.reflected(sampSlc, lower=xmn, upper=xmx) 
  w <- approxfun(densSamp$x, densSamp$y)
  a <- min(densSamp$x)
  b <- max(densSamp$x)
  
  plotL <- list()
  for (s in spp){
    spRows <- which(slc$species==s)
    spSlc <- slc[spRows,env]
    z1 <- tryCatch(
      transformEst(spSlc, w = w, bw='SJ-ste', reflect = FALSE, a = a, b = b),
      error = function(err){ list() }
    ) 
    if (length(z1)==0){ next }
    
    pe <- nichStats(z1)['pe']
    x <- z1$x
    kd <- z1$y
    x <- c(xmn, x, xmx)
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
      geom_segment(x=pe, xend=pe, y=0, yend=1) +
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
  kdeNm <- paste0('Figs/KDE_all_spp_transform-SJ_refl_',bin,'ka_',day,'.pdf')
  pdf(kdeNm, width=8.5, height=11)
  print(pg)
  dev.off()
}
