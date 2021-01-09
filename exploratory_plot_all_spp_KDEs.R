library(ggplot2)
library(cowplot)
library(sp)
library(raster)
library(testthat)
library(kerneval)

day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")
df <- readRDS('Data/spp-and-sampling-data_list-by-depth_2020-11-15.rds')

# for the main text figure, put n-label in title not on panel
# but do plot lines to indicate optimum and mean MAT
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
  pePos <- which.max(z1$y)
  pe <- z1$x[pePos]
  
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
    theme_classic() +
    ggtitle(paste(bin, 'ka, n =', length(spDat))) +
    geom_area(fill = 'orange') +
    geom_line() +
    geom_segment(x = pe, xend = pe, y = 0, yend = 1,
                colour = 'orange3') +
    geom_segment(x = mean(spDat), xend = mean(spDat), 
                 y = 0, yend = 1, colour = 'grey20') +
    geom_segment(x = -1.4, xend = -1.4, linetype='dashed',
                 y = 0, yend = 1, colour = 'grey50') +
    geom_segment(x = 27.1, xend = 27.1, linetype='dashed',
                 y = 0, yend = 1, colour = 'grey50') +
    geom_rug(aes(x = rugHack), colour = 'grey20', size = 0.4,
             sides = 'b', length = unit(0.045, "npc")) + 
    geom_hline(yintercept = 0) +
    # ,outside = TRUE) +
    #   coord_cartesian(clip = "off") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 0.05)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 9))
  # sPlot <- sPlot +
  #   geom_text(aes(fontface = 'plain', family = 'mono'), 
  #             label = nlab, size = 3.5, # 3.5mm = 10pt
  #             x = 0.15 * max(plotDat$x), y = 0.07)
  sPlot
}

# Figure 1 ----------------------------------------------------------------

modCodes <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors = FALSE)
rAll <- list.files('Data/gcm_annual_mean/')
rTemp <- rAll[grep('temp', rAll)]

# set blue as low temp value
colr <- rainbow(40)[c(30:1)]

b1 <- 28
b2 <- 212
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
  sPlot <- plotDens(bin = b, s = pachy, bw = 'SJ-ste',
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

# calculate H distance

# species samples
sBool1 <- spp0m$species == pachy & spp0m$bin == b1
spDat1 <- spp0m$temp_ym[sBool1]
sBool2 <- spp0m$species == pachy & spp0m$bin == b2
spDat2 <- spp0m$temp_ym[sBool2]

# estimate sampling correction
sampBool1 <- samp0m$bin == b1
sampDat1 <- samp0m$temp_ym[sampBool1] 
densSamp1 <- density(sampDat1, bw = 'SJ-ste') 
w1 <- approxfun(densSamp1$x, densSamp1$y)
sampBool2 <- samp0m$bin == b2
sampDat2 <- samp0m$temp_ym[sampBool2] 
densSamp2 <- density(sampDat2, bw = 'SJ-ste') 
w2 <- approxfun(densSamp2$x, densSamp2$y)

xmn <- 1.4
xmx <- 27.1
d1 <- transdens(spDat1, w = w1, reflect = FALSE, 
                a = xmn, b = xmx, bw = 'SJ-ste')
d2 <- transdens(spDat2, w = w2, reflect = FALSE, 
                a = xmn, b = xmx, bw = 'SJ-ste')
h <- hell(d1, d2) 
round(h, 2)
# [1] 0.08

# Glacial/IG niches -------------------------------------------------------
# Supplemental Figure 1

coldBins <- c(28, 268, 436)
hotBins <- c(4, 124, 412)

# pick 1 surface and 1 subsurface species as examples
spp <- c('Neogloboquadrina dutertrei','Hirsutella scitula')

sppSurf  <- df$temp_ym_surf$sp
sampSurf <- df$temp_ym_surf$samp
sppSub   <- df$temp_ym_sub$sp
sampSub  <- df$temp_ym_sub$samp

plotL <- list()  

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
  z1 <- transdens(spDat, w = w, bw = bw, reflect = FALSE, 
                  a = xmn, b = xmx, n = 2^10)
  xlim1 <- min(sampDf$temp_ym)
  xlim2 <- max(sampDf$temp_ym)
  
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
    geom_segment(x = xmn, xend = xmn, linetype='dashed',
                 y = 0, yend = 1, colour = 'grey50') +
    geom_segment(x = xmx, xend = xmx, linetype='dashed',
                 y = 0, yend = 1, colour = 'grey50') +
    geom_rug(aes(x = rugHack), colour = 'grey20', size = 0.4,
             sides = 'b', length = unit(0.045, "npc")) + 
    geom_hline(yintercept = 0) +
    scale_x_continuous(expand = c(0, 0), limits = c(xlim1, xlim2)) +
    scale_y_continuous(limits = c(0, 0.12)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 9))
  sPlot <- sPlot +
    geom_text(aes(fontface = 'plain', family = 'mono'), 
              label = nlab, size = 3.5, # 3.5mm = 10pt
              x = 0.15 * xlim2, y = 0.07)
  sPlot
}

for (s in spp){
  if (s=='Neogloboquadrina dutertrei'){
    sppDf <- sppSurf
    sampDf <- sampSurf
    xmn <- -1.3
    xmx <- 26.9
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

# H-value reference examples ----------------------------------------------
# Supplemental Figure 2

normPdf <- function(x, shft){
  dnorm(x, (25 + shft), sd=9.7)
}

cauchyPdf <- function(x, shft){
  shft <- 1.5^(shft)*50/1.5^(50)
  dcauchy(x, location=(25+shft), scale=0.39)
}
# With the given location and scale, 99% of Cauchy distribution mass falls from 0-50 
# p <- 0.99 + pcauchy(0, location=25, scale=0.39)
# qcauchy(p, location=25, scale=0.39)

chsqPdf <- function(x, shft){
  dchisq((x-shft), df=4)
  # dchisq((x-shft)/3.76, df=4)
}
# The chi-sq distribution with 4 degrees of freedom is short-tailed, 
# so one needn't shift a paired chi-sq distribution as far away to achieve non-overlap
# qchisq(0.99, df=4)*3.76

# mixture of normal distributions, as in Chiu 1996, but this one is less peaky
mixPdf <- function(x, shft){
  d1 <- 0.75 * dnorm(x, (27+shft), 6)
  d2 <- 0.25 * dnorm(x, (11+shft), 4)
  d1 + d2
}

# Return H for a given pair of niche shapes (PDFs) and positions.
overlapr <- function(fun1, fun2, shft, a, b){
  pdf1 <- switch(fun1, 'norm'=normPdf, 'cauchy'=cauchyPdf, 'chsq'=chsqPdf, 'mix'=mixPdf)
  pdf2 <- switch(fun2, 'norm'=normPdf, 'cauchy'=cauchyPdf, 'chsq'=chsqPdf, 'mix'=mixPdf)
  intgrnd <- function(x) sqrt(pdf1(x, 0) * pdf2(x, shft))
  int <- integral(intgrnd, a, b, no_intervals = 12) 
  # for norm vs norm, shft=10, H is .024 instead of 0.35 if integrated over -Inf to Inf
  H2 <- (1 - int)
  #  I <- 1 - H^2
  #  data.frame(fun1, fun2, shft, H)
  sqrt(H2)
}

shftVals <- seq(0, 50, by=0.5)
Hseq <- vector()
for(s in shftVals){
  H <- overlapr('norm', 'cauchy', s, a = 0, b = 100)
  Hseq <- c(Hseq, H)
}
normCauch <- data.frame(shift = shftVals, H = Hseq)

# *plots -------------------------------------------------------------------

u <- 1:100
norm1 <- normPdf(u, 0)
plotEmpt <- ggplot() +
  theme_classic() +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.05), expand = c(0, 0),
                     name = 'Density') +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# for H = 0.25, shift norm dist by 6.9
H1 <- overlapr('norm','norm', 6.9, a = 0, b = 100)
ttl1 <- paste('H =', round(H1, 2))
norm2 <- normPdf(u, 6.9)
df1 <- data.frame(x = u, pdf1 = norm1, pdf2 = norm2)
plot1 <- plotEmpt +
  ggtitle(ttl1) +
  geom_line(data = df1, aes(x = x, y = pdf1), size = 1) +
  geom_line(data = df1, aes(x = x, y = pdf2), size = 1) 
plot1

# for H = 0.5, shift norm dist by 14.7
norm3 <- normPdf(u, 14.7)
H2 <- overlapr('norm','norm', 14.7, a = 0, b = 100)
ttl2 <- paste('H =', sprintf('%0.2f', round(H2, 2)))
df2 <- data.frame(x = u, pdf1 = norm1, pdf2 = norm3)
plot2 <- plotEmpt +
  ggtitle(ttl2) +
  geom_line(data = df2, aes(x = x, y = pdf1), size = 1) +
  geom_line(data = df2, aes(x = x, y = pdf2), size = 1) 
plot2

# norm vs. mixed distribution, H = 0.25
mix <- mixPdf(u, 7.6)
H3 <- overlapr('norm','mix', 7.6, a = 0, b = 100)
ttl3 <- paste('H =', round(H3, 2))
df3 <- data.frame(x = u, pdf1 = norm1, pdf2 = mix)
plot3 <- plotEmpt +
  ggtitle(ttl3) +
  geom_line(data = df3, aes(x = x, y = pdf1), size = 1) +
  geom_line(data = df3, aes(x = x, y = pdf2), size = 1) 
plot3

# norm vs. chi-sq, H = 0.57
chsq <- chsqPdf(u, 19.5)
H4 <- overlapr('norm','chsq', 19.5, a = 0, b = 100)
ttl4 <- paste('H =', round(H4, 2))
df4 <- data.frame(x = u, pdf1 = norm1, pdf2 = chsq)
plot4 <- plotEmpt +
  ggtitle(ttl4) +
  geom_line(data = df4, aes(x = x, y = pdf1), size = 1) +
  geom_line(data = df4, aes(x = x, y = pdf2), size = 1) +
  scale_y_continuous('Density', limits = c(0, 0.2), expand = c(0, 0))
plot4

day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")
panelNm <- paste0('H-reference-values_4-panels_', day, '.pdf')
pdf(panelNm, width = 5, height = 5)
plot_grid(plot1, plot2, plot3, plot4, 
          ncol = 2, scale = 0.9,
          labels = 'AUTO', label_size = 14,
          label_x = 0.05, vjust = 2.5
)
dev.off()
