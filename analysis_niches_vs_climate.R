library(ggplot2)
library(fUnitRoots) 
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(cowplot)

ss <- TRUE

# Data prep ---------------------------------------------------------------
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

# foram data
spAttr <- read.csv('Data/foram_spp_data_20-04-05.csv', stringsAsFactors=FALSE)
if (ss){
  df <- read.csv('Data/foram_niche_sumry_metrics_0m_20-04-05.csv', stringsAsFactors=FALSE)
} else {
  df <- read.csv('Data/foram_niche_sumry_metrics_20-04-05.csv', stringsAsFactors=FALSE)
}
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]
bins <- unique(df$bin)
spp <- unique(df$sp)
nspp <- length(spp)
binL <- bins[1] - bins[2]
nCore <- detectCores() - 1

# standardized sampling universe (MAT at range-through core sites) at each of 4 depths
truncEnv <- readRDS('Data/sampled_temp_ym_truncated_by_depth_20-04-05.rds')
envNm <- 'temp_ym'

# mean MAT over the globe, at a standard grid of lat-long points
glob <- read.csv('Data/global_surface_MAT_at_grid_pts_4ka.csv')
cols <- paste0('X',bins)
globMean <- colMeans(glob[,cols])

# Scatterplot -------------------------------------------------------------

# Mean species H overlap vs delta global MAT in each bin.
# H does not indicate direction, only magnitude of niche overlap,
# so should compare it with absolute differences in global MAT.

sumH <- function(b){
  bBool <- df$bin==b
  slc <- df[bBool,]
  avg <- apply(slc[,c('h','pe','m')], 2, mean, na.rm=TRUE)
  binN <- nrow(slc)
  c(bin=b, avg['h'], avg['pe'], avg['m'], nSpp=binN)
}
bybin <- sapply(bins, sumH)
bybin <- data.frame(t(bybin))
bybin$glob <- globMean

# all H values are NA at most recent time step
Hseq <- bybin[-nrow(bybin),]

# delta MAT time series is stationary but H series is NOT
delta <- diff(globMean)
absDelta <- abs(delta)
Hseq$absDelta <- absDelta
adfTest(absDelta) 
acf(absDelta) 
adfTest(Hseq$h)
plot(acf(Hseq$h), ci.type="ma")
# account for non-stationarity and autocorrelation
arH <- arima(Hseq$h, order=c(1,0,0))
resid <- as.numeric(arH$residuals)
adfTest(resid)
acf(resid)
arH$coef
lmH <- lm(resid ~ absDelta) 
acf(lmH$residuals)
cor.test(resid, absDelta, method='pear')

xmx <- max(Hseq$absDelta) * 1.1
scatrCor <- cor(resid, absDelta, method='pear')
scatrLab <- paste('r =', round(scatrCor, 2))
deltaPlot <- 
  ggplot(data=Hseq, aes(x=absDelta, y=h)) +
  theme_bw() +
  scale_y_continuous('Mean Hellinger\'s H', limits=c(0, 0.52), expand=c(0,0)) +
  scale_x_continuous('Magnitude of global MAT change (C)',
                     limits=c(0,xmx), expand=c(0,0)) +
  geom_point() 
#  theme(axis.title.y = element_blank())

deltaPlot <- deltaPlot +
  geom_text(label=scatrLab, size=3, x=xmx*0.8, y=0.45)

# Global time series ------------------------------------------------------

# contrast the niche overlap between extreme situations:
# glaciation peak and terminus, for 8 cycles
# (compared to peak vs. peak and terminus vs. terminus)

# find the local max and min global MAT timing in each 100ky interval
ints <- data.frame(yng=c(seq(0,400,by=100), 480, 560), # 690
                   old=c(seq(100,400,by=100), 480, 560, 690) # 800
)
#ints$maxAge <- ints$minAge <- NA
for (r in 1:nrow(ints)){
  bounds <- ints[r,]
  inInt <- which(bins > bounds$yng & bins <= bounds$old)
  minPos <- which.min(globMean[inInt])
  minAge <- bins[inInt[minPos]]
  maxPos <- which.max(globMean[inInt])
  maxAge <- bins[inInt[maxPos]]
  range(globMean[inInt])
  ints[r,c('minAge','maxAge','minT','maxT')] <- 
    c(minAge, maxAge, range(globMean[inInt]))
}

# plot time series of global MAT
globDat <- data.frame(bins, globMean)
globTseries <- ggplot() +
  theme_bw() +
  scale_y_continuous('Global MAT (C)') +
  scale_x_continuous('Time (Ka)', expand=c(0,0),
                     limits=c(-700,0), breaks=seq(-700,0,by=100),
                     labels=paste(seq(700,0,by=-100))) +
  geom_line(data=globDat, aes(x=-bins, y=globMean)) +
  geom_point(data=globDat, aes(x=-bins, y=globMean)) +
  geom_point(data=ints, aes(x=-minAge, y=minT), colour='deepskyblue', size=2) +
  geom_point(data=ints, aes(x=-maxAge, y=maxT), colour='firebrick2', size=2)

# Extreme comparisons -----------------------------------------------------

# prepare a framework of every pairwise comparison to compute
# (every warm vs. cold, warm vs. warm, and cold vs. cold interval)
wc <- expand.grid(X1=ints$minAge, X2=ints$maxAge)
wc$type <- 'cold-warm'
cc <- combn(ints$minAge, 2)
cc <- data.frame(t(cc))
cc$type <- 'cold-cold'
ww <- combn(ints$maxAge, 2)
ww <- data.frame(t(ww))
ww$type <- 'warm-warm'
intPairs <- rbind(wc, cc, ww) 
colnames(intPairs) <- c('t1','t2','type')
comps <- c('cold-cold','warm-warm','cold-warm')
intPairs$type <- factor(intPairs$type, levels = comps)
colr <- c('cold-cold'="deepskyblue",
          'warm-warm'="firebrick2",
          'cold-warm'="purple3")

# inspect the mean delta MAT for each comparison type
globDiff <- function(pair){
  nm1 <- paste0('X', pair['t1'])
  nm2 <- paste0('X', pair['t2'])
  abs(globMean[nm1] - globMean[nm2])
}
intPairs$deltaMAT <- apply(intPairs[,c('t1','t2')], 1, globDiff)
for (typ in comps){
  compBool <- intPairs$type==typ
  compDelt <- intPairs$deltaMAT[compBool]
  mDelt <- round(mean(compDelt),2)
  print(paste(typ, mDelt))
}

source('GSA_custom_ecospat_fcns.R')
source('species_kde_buildr.R')

pairL <- list()
for (i in 1:nrow(intPairs)){
  b1 <- intPairs$t1[i]
  b2 <- intPairs$t2[i]
  entry <- c(b1, b2)
  pairL <- append(pairL, list(entry))
}

# warning - this could take 1 - 2 hours
# pkg <- c('pracma','GoFKernel')
# pt1 <- proc.time()
# registerDoParallel(nCore)
# if (ss){
#  kdeSum <- foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar% 
#   kde(truncEnv[[1]], bPair, envNm)
# } else {
# kdeSum <- foreach(dat=truncEnv[2:4], .combine=rbind, .inorder=FALSE, .packages=pkg) %:% 
#  foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar% 
#  kde(dat, bPair, envNm)
# }
# stopImplicitCluster()
# pt2 <- proc.time()
# pt2-pt1
# nas <- is.na(kdeSum$bin)
# kdeSum <- kdeSum[!nas,]
# if (ss){
#  sumNm <- paste0('Data/foram_niche_xtreme_comparisons_0m_',day,'.csv')
# } else {
#  sumNm <- paste0('Data/foram_niche_xtreme_comparisons_hab_',day,'.csv')
# }
# write.csv(kdeSum, sumNm, row.names = FALSE)

# if the script has already been run once, read in the intermediate products instead
if (ss){
 kdeSum <- read.csv('Data/foram_niche_xtreme_comparisons_0m_20-04-07.csv')
} else {
 kdeSum <- read.csv('Data/foram_niche_xtreme_comparisons_hab_20-04-07.csv')
}

# take the mean h for each bin combination
sumH <- function(bPair, dat){
  intBool <- which(dat$bin==bPair[1] & dat$bin2==bPair[2])
  int <- dat[intBool,]
  avgH <- mean(int$h, na.rm=TRUE)
  data.frame(t1=bPair[1], t2=bPair[2], avgH)
}
Hlist <- lapply(pairL, sumH, dat=kdeSum)
Hdf <- do.call(rbind, Hlist)

intPairs <- merge(intPairs, Hdf)
#mxH <- max(intPairs$avgH) * 1.1
ovpBoxs <- ggplot(data=intPairs) +
  theme_bw() +
  scale_y_continuous(limits=c(0,0.52), expand=c(0,0)) +
  geom_boxplot(aes(x=type, y=avgH, fill=type)) +
  scale_fill_manual(values=colr) +
  theme(axis.title=element_blank(),
        axis.text.y = element_blank(),
        legend.position='none')

# Multipanel plots --------------------------------------------------------

if (ss){
  # main text figure
  
  # align bottom 2 plot panels vertically, then top and bottom panel horizontally
  alignLft <- align_plots(globTseries, deltaPlot, align = 'v', axis = 'l')
  alignLwr <- align_plots(alignLft[[2]], ovpBoxs, align='h')
  lwrRow <- plot_grid(
    alignLwr[[1]], alignLwr[[2]],
    ncol = 2,
    rel_widths = c(1.13,1),
    labels=c('B','C'),
    label_size = 12,
    label_x = c(0.16,0.05),
    vjust=2.3
  )
  mlti <- plot_grid(alignLft[[1]], lwrRow, 
                    ncol=1, 
                    labels=c('A',''), 
                    label_size = 12,
                    label_x = 0.085, 
                    vjust=2.3)
  
  panelsNm <- paste0('Figs/H-vs-climate_panels_main_',day,'.pdf')
  pdf(panelsNm, width=7, height=4)
  print(mlti)
  dev.off()
  
} else {
  # make supplemental figure: panels B and C only
  alignLwr <- align_plots(deltaPlot, ovpBoxs, align='h')
  lwrRow <- plot_grid(
    alignLwr[[1]], alignLwr[[2]],
    ncol = 2,
    rel_widths = c(1.13,1),
    labels='AUTO',
    label_size = 12,
    label_x = c(0.17,0.05),
    vjust=2.3
  )
  panelsNm <- paste0('Figs/H-vs-climate_panels_SI_',day,'.pdf')
  pdf(panelsNm, width=6, height=3)
  print(lwrRow)
  dev.off()
}
