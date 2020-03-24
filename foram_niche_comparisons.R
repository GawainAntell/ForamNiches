library(PerformanceAnalytics)
library(paleoTS)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

doTrunc <- TRUE

# Data prep ---------------------------------------------------------------
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y-%m-%d')

if (doTrunc){
  df <- read.csv('Data/foram_niche_sumry_metrics_trunc_20-03-24.csv', stringsAsFactors=FALSE)
  samp <- read.csv('Data/samp_MAT_occs_latlong_8ka_trunc_20-03-24.csv',
                   stringsAsFactors = FALSE)
  v <- 'temp_ym_0m'
  # See note on niche construction script: if truncating temp by surface values,
  # it doesn't make sense to use in-habitat temperature downstream.
} else {
  df <- read.csv('Data/foram_niche_sumry_metrics_20-03-24.csv', stringsAsFactors=FALSE)
  samp <- read.csv('Data/samp_MAT_occs_latlong_8ka_20-03-24.csv',
                   stringsAsFactors = FALSE)
  v <- 'temp_ym_0m' # 'temp_ym_hab'
}
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]
vShrt <- paste(strsplit(v,'_')[[1]], collapse='')
bins <- unique(df$bin)
spp <- unique(df$sp)
binL <- bins[1] - bins[2]
nspp <- length(spp)

# Deal with incomplete sp ts ----------------------------------------------

spL <-  vector('list', length = nspp)
for (i in 1:nspp){
  s <- spp[i]
  spBool <- df$sp==s
  sp <- df[spBool,]
  j <- 1
  sp$sp[1] <- paste0(s, j)
  penult <- nrow(sp) - 1
  for (r in 1:penult){
    age <- sp$bin[r]
    nextAge <- sp$bin[r+1]
    if (age !=  nextAge+binL){
      j <- j + 1
    }
    sp$sp[r+1] <- paste0(s, j)
  }
  spL[[i]] <- sp
}
splitSpp <- do.call(rbind, spL)

# Check for per-species continuity - ditch any lineage chuncks of < 7 time bins
# Could reduce this to 6 bins but would have to modify paleoTS code
# TODO is this necessary anymore? should be taken care of in revised data prep script
longSpp <- character()
enuf <- rep(binL, 7)
enufTxt <- paste0(enuf, collapse='')
epithets <- unique(splitSpp$sp)
for (s in epithets){
  spBool <- splitSpp$sp==s
  spDf <- splitSpp[spBool,]
  spB <- sort(unique(spDf$bin))
  bDiff <- diff(spB)
  diffTxt <- paste0(bDiff, collapse='')
  srch <- grep(enufTxt,diffTxt)
  if (length(srch) > 0){
    longSpp <- c(longSpp, s)
  }
}
keepBool <- splitSpp$sp %in% longSpp
nich <- splitSpp[keepBool,]

# Sampled vs global MAT corr ----------------------------------------------

# calculate mean MAT over globe, at set grid of points
glob <- read.csv('Data/global_surface_MAT_at_grid_pts_4ka.csv')
cols <- paste0('X',bins)
globMean <- colMeans(glob[,cols])

# combine with sampling n and mean MAT
sampSumry <- sapply(bins, function(b){
  slcBool <- samp$b==b
  slc <- samp[slcBool,]
  temp <- mean(slc$temp_ym_0m)
  n <- nrow(slc)
  c(temp,n)
})
sampSumry <- data.frame(t(sampSumry))
colnames(sampSumry) <- c('mat','n')

# autocorrelation in residuals of main relationship of interest
l <- lm(sampSumry$mat~globMean)
r <- resid(l)
acf(r)

corVars <- cbind(n=sampSumry$n, sampMean=sampSumry$mat, globMean)

# There are WAY more sampling points in the last 2 time steps -
# these are very influential points, so exclude from cor tests.
# But could also include these points - they aren't really outliers.
excl <- nrow(corVars) - 0:1
corVars <- corVars[-excl,]

if (doTrunc){
  corNm <- paste0('Figs/correlations_sampled_vs_global_trunc_',day,'.pdf')
} else {
  corNm <- paste0('Figs/correlations_sampled_vs_global_',day,'.pdf')
}
pdf(corNm, width=5, height=5)
chart.Correlation(corVars, histogram=TRUE, pch=19, method='pearson')
dev.off()

# Evo models --------------------------------------------------------------
# TODO: use sampled environment from same depth as species, tsSamp object

# alternatively, load saved versions and skip running this section
# (if running 9 models, this is slow, so it's worthwhile to save and load)
# mods4 <- readRDS("Data/Data/4_evo_model_fits_trunc_tempym0m200214.rds")
# mods4 <- readRDS("Data/Data/4_evo_model_fits_tempymhab200214.rds")

evoFit <- function(s, sampStats, nmods='four'){
  mNm <- sampStats[1]
  sdNm <- sampStats[2]
  spBool <- nich$sp==s
  sp <- nich[spBool,]
  # also subset the sampling time series to the same bins, for comparison
  sBins <- sp$bin
  sampBool <- nich$bin %in% sBins & nich$sp=='sampled1'
  samp <- nich[sampBool,]
  
  # ages must start at 0
  samp$scaledT <- sp$scaledT <- 1:nrow(sp) -1
  
  # save metadata about the time series
  l <- nrow(sp)
  strt <- sp$bin[1]
  
  ts <- as.paleoTS(mm = sp[,mNm], vv = sp[,sdNm]^2, 
                   nn = sp$n, tt = sp$scaledT, 
                   oldest = 'first', reset.time = FALSE)
  tsSamp <- as.paleoTS(mm = samp$m_temp_ym_0m, vv = samp$sd_temp_ym_0m^2, 
                       nn = samp$n, tt = samp$scaledT, 
                       oldest = 'first', reset.time = FALSE)
  
  # some sample variances are not equal, so use pool=FALSE for all sequences
  if (l < 14 | nmods=='four'){
    modsSp <- fit4models(ts, method='AD', silent=TRUE, pool=FALSE) 
    modsSamp <- fit4models(tsSamp, method='AD', silent=TRUE, pool=FALSE) 
  } else {
    modsSp <- fit9models(ts, method='AD', silent=TRUE, pool=FALSE)
    modsSamp <- fit9models(tsSamp, method='AD', silent=TRUE, pool=FALSE)
  }
  # From the package documentation:
  # 'Method = "Joint" is a full likelihood approach, considering each time-series as a joint sample
  # from a multivariate normal distribution. Method = "AD" is a REML approach that uses the differences
  # between successive samples. They perform similarly, but the Joint approach does better under
  # some circumstances (Hunt, 2008).'
  # From Hunt 2008:
  # 'joint parameterization  is  better  able  to  correctly  identify  directional  trends. 
  # This advantage increases with sequence length, and is most pronounced when sampling error is high'
  
  # HOWEVER, although joint might be better for most series, for Pulleniatina obliquiloculata3 
  # a fatal error occurs (matrix not full rank?) unless pool and hess arguments are modified
  # Error in optim(p0, fn = logL.joint.URW, control = cl, method = meth, lower = c(NA,  : 
  # non-finite value supplied by optim
  
  wts <- modsSp$modelFits$Akaike.wt
  maxMod <- which.max(wts)
  modNm <- row.names(modsSp$modelFits)[maxMod]
  w <- max(wts)
  params <- modsSp$parameters[modNm][[1]]
  
  # check what model would be predicted by sampling alone for the given time bins
  sampMx <- which.max(modsSamp$modelFits$Akaike.wt)
  sampEvo <- row.names(modsSamp$modelFits)[sampMx]
  
  out <- data.frame(sp=s, var=mNm, bestMod=modNm, weight=w, l=l, start=strt, samplingMod=sampEvo)
  out$params <- list(params)
  out
}

statsNm <- paste(c('m','sd'), v, sep='_')

mods4l <- lapply(longSpp, evoFit, sampStats = statsNm, nmods = 'four') # 'nine'
mods4 <- do.call('rbind', mods4l)

if (doTrunc){
  mods4nm <- paste0('Data/4_evo_model_fits_trunc_',vShrt, day, '.rds')
} else {
  mods4nm <- paste0('Data/4_evo_model_fits_',vShrt, day, '.rds')
}
# saveRDS(mods4, mods4nm)

table(mods4$bestMod)

# Spp vs sampling evo mode ------------------------------------------------

modsSpBool <- mods4$sp != 'sampled1'
modsSp <- mods4[modsSpBool,]

evoModes <- c('StrictStasis','Stasis','URW','GRW')
# if using 9 models, add 'Punc-1','Stasis-URW','Stasis-GRW','URW-Stasis','GRW-Stasis'
xy <- expand.grid(spMode=evoModes, sampMode=evoModes, stringsAsFactors=FALSE)
xy$n <- NA
for (i in 1:nrow(xy)){
  x <- xy$sampMode[i]
  y <- xy$spMode[i]
  same <- which(modsSp$samplingMod==x & modsSp$bestMod==y)
  n <- length(same)
  xy$n[i] <- n
}
xy$spMode <- factor(xy$spMode, levels=evoModes)
xy$sampMode <- factor(xy$sampMode, levels=evoModes)
empty <- xy$n==0
xy <- xy[!empty,]

bubbl <- 
  ggplot(data=xy, aes(x=sampMode, y=spMode, size=n)) +
  theme_bw() +
  geom_text(aes(label=n), nudge_x=0.05, nudge_y=0.05) +
  scale_x_discrete(name = 'Sampling model', drop=FALSE) +
  scale_y_discrete(name = 'Species model', drop=FALSE) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        legend.position = 'none') +
  scale_size(range=c(4,8))

if (doTrunc){
  bubblNm <- paste0('Figs/evo_mode_bubble_matrix_trunc_', vShrt, day, '.pdf')
} else {
  bubblNm <- paste0('Figs/evo_mode_bubble_matrix_', vShrt, day, '.pdf')
}
pdf(bubblNm, width=4, height=4)
print(bubbl)
dev.off()


# Intra-sp niche overlap --------------------------------------------------

consec <- ! is.na(df$h)
intraH <- df[consec,]
intraH$shortNm <- sapply(intraH$sp, function(txt){
  splt <- strsplit(txt, ' ')
  splt[[1]][2]
} )

intra <- 
  ggplot(data=intraH, aes(x=shortNm, y=h)) +
  scale_y_continuous('Hellinger\'s H', limits=c(0,1), expand=c(0,0)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())

if (doTrunc){
  ovrlpNm <- paste0('Figs/overlap_H_boxplots_by_species_trunc_',vShrt,day,'.pdf')
} else {
  ovrlpNm <- paste0('Figs/overlap_H_boxplots_by_species_',vShrt,day,'.pdf')
}
pdf(ovrlpNm, width=6, height=4)
intra
dev.off()

# * Intra overlap vs evo var ----------------------------------------------

getV <- function(x){
  xmod <- x$bestMod
  if (xmod=='StrictStasis'){
    v <- 0
  }
  if (xmod=='Stasis'){
    plist <- x['params']
    v <- plist[[1]]['omega']
  }
  if (xmod %in% c('URW','GRW')){
    plist <- x['params']
    v <- plist[[1]]['vstep']
  }
  return(v)
}
mods4$stepvar <- apply(mods4, 1, getV)

# subset the model-summary dataframe to the focal environmental variable,
v4rows <- grep(v, mods4$var)
mods4v <- mods4[v4rows,]

# summarise H for the same sequences as evo models are fit
# (the other plots consider H even from intervals outside sequences of 7+ bins)
getH <- function(s){
  sBool <- nich$sp==s
  sDf <- nich[sBool,]
  median(sDf[,yNm], na.rm=TRUE)
}
longSpp <- sort(longSpp)
seqH <- sapply(longSpp, getH)
hDf <- data.frame(sp=longSpp, h=seqH)

mods4v <- merge(mods4v, hDf, 'sp')
vh <- 
  ggplot(mods4v, aes(x=h, y=stepvar)) +
  scale_y_continuous(name='estimated evo model variance') +
  scale_x_continuous(name='median H among time intervals') +
  geom_point(aes(col=bestMod)) +# , size=3
  theme_bw()

if (doTrunc){
  vhNm <- paste0('Figs/overlap_H_vs_evo_mod_variance_trunc_',vShrt,day,'.pdf')
} else {
  vhNm <- paste0('Figs/overlap_H_vs_evo_mod_variance_',vShrt,day,'.pdf')
}
pdf(vhNm, width=5, height=4)
print(vh)
dev.off()

# * Overlap by ecomorph ---------------------------------------------------

# summarise niche data for each species as the median among intervals, to plot
spSmry <- function(s, dat){
  sBool <- dat$sp==s
  sDf <- dat[sBool,]
  numCol <- ! colnames(dat) %in% c('bin','sp','n','shortNm','class')
  out <- apply(sDf[,numCol], 2, median, na.rm=TRUE)
  out <- data.frame(t(out))
  out$nInt <- nrow(sDf)
  out$species <- s
  out
}
medsL <- lapply(spp, spSmry, dat=intraH)
medsDf <- do.call(rbind, medsL)

# combine with species attribute data
atts <- read.csv('Data/foram_spp_data_200108.csv', stringsAsFactors=FALSE)
plotDf <- merge(medsDf, atts, by.x='species')
plotDf$spinose <- as.factor(plotDf$spinose)
plotDf$DepthHabitat <- factor(plotDf$DepthHabitat, 
                              levels=c('Subsurface','Surface.subsurface','Surface'))
plotDf$eco <- factor(plotDf$eco, levels=paste(5:1))
ecoLabs <- c('Mixed layer, symb',
             'Mixed layer, no symb',
             'Thermocline',
             'Sub-thermocline',
             'High latitude')
# From Aze et al. 2011, ecotype codes are:
# 1 = open ocean, mixed layer, trop/subtrop, w symbionts
# 2 = open ocean, mixed layer, trop/subtrop, w/o symbionts
# 3 = open ocean, thermocline
# 4 = open ocean, sub-thermocline
# 5 = high lat
# 6 = upwelling/high productivity

ymx <- max(plotDf$h)*1.05
empt <- ggplot(data=plotDf, aes(y=h)) +
  theme_bw() +
  scale_y_continuous(limits=c(0,ymx), expand=c(0,0)) +
  theme(axis.title.x = element_blank())

# add n-labels aligned at right of plot area
lablr <- function(v){
  lab <- paste('n =', length(v))
  data.frame(y=0.025, label=lab)
}

p1 <- 
  empt +
  geom_boxplot(aes(x=DepthHabitat), fill='grey') +
  scale_x_discrete(name='Preferred depth') +
  theme(axis.text.x = element_blank()) +
  stat_summary(
    aes(x=DepthHabitat),
    fun.data = lablr, 
    geom = "text"
  ) + 
  coord_flip()

p2 <- empt +
  geom_boxplot(aes(x=eco), fill='grey') +
  # reverse the order of labels because of the axis flip
  scale_x_discrete(name='Ecotype', labels=rev(ecoLabs)) +
  theme(axis.text.x = element_blank()) +
  stat_summary(
    aes(x=eco),
    fun.data = lablr, 
    geom = "text"
  ) + 
  coord_flip()

p3 <- empt +
    geom_boxplot(aes(x=spinose), fill='grey') +
    scale_x_discrete(name='Morphotype') +
    stat_summary(
      aes(x=spinose),
      fun.data = lablr, 
      geom = "text"
    ) + 
    coord_flip()

smallMult <- plot_grid(p1, p2, p3,
                       rel_heights=c(0.7, 1, 0.5),
                       ncol=1, align='v') 

if (doTrunc){
  multNm <- paste0('Figs/overlap_H_by_ecomorph_trunc_',vShrt,day,'.pdf')
} else {
  multNm <- paste0('Figs/overlap_H_by_ecomorph_',vShrt,day,'.pdf')
}
pdf(multNm, width=4, height=7)
grid.arrange(arrangeGrob(smallMult), bottom='Mean Hellinger\'s H among species') 
dev.off()

# * Delta global MAT vs overlap -------------------------------------------
# note that H overlap does not indicate direction, only magnitude, of niche change
# so should compare it with absolute differences in global MAT
sumH <- function(b){
  bBool <- intraH$bin==b
  slc <- intraH[bBool,]
  avgH <- mean(slc[,'h'])
  binN <- nrow(slc)
  c(bin=b, avgH=avgH, nSpp=binN)
}
Hseq <- sapply(bins, sumH)
Hseq <- data.frame(t(Hseq))
# all values NA at most recent time step
Hseq <- Hseq[-nrow(Hseq),]

delta <- diff(globMean)
Hseq$absDelta <- abs(delta)

ymx <- max(Hseq$avgH) * 1.1
xmx <- max(Hseq$absDelta) * 1.1
scatrCor <- cor(Hseq$absDelta, Hseq$avgH, method='spear')
scatrLab <- paste('rho =', round(scatrCor, 3))
deltaPlot <- 
  ggplot(data=Hseq, aes(x=absDelta, y=avgH)) +
  theme_bw() +
  scale_x_continuous('Absolute change in global MAT (C)',
                     limits=c(0,xmx), expand=c(0,0)) +
  scale_y_continuous('Mean H among boundary-crossers', 
                     limits=c(0,ymx), expand=c(0,0)) +
  geom_point()
deltaPlot <- deltaPlot +
  geom_text(aes(fontface=1), label=scatrLab,# size=3,
            x=xmx*0.8, y=ymx*0.8)
if (doTrunc){
  deltaNm <- paste0('Figs/delta-MAT-vs-mean-H_trunc_',vShrt,'_to',cutoff,'ka_',day,'.pdf')
} else {
  deltaNm <- paste0('Figs/delta-MAT-vs-mean-H_',vShrt,'_',day,'.pdf')
}
pdf(deltaNm, width=6, height=4)
print(deltaPlot)
dev.off()

# * G-IG overlap comparisons ----------------------------------------------

# contrast the niche overlap between extreme situations:
# glaciation peak and terminus, for 8 cycles
# (compared to peak vs. peak and terminus vs. terminus)

# find the local max and min global MAT timing in each 100ky interval
ints <- data.frame(yng=c(seq(0,400,by=100), 480, 560, 690),
                   old=c(seq(100,400,by=100), 480, 560, 690, 800)
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
  scale_y_continuous('global MAT (C)') +
  scale_x_continuous('time (ka), 8ky intervals', expand=c(0,0),
                     limits=c(-800,0), breaks=seq(-800,0,by=100),
                     labels=paste(seq(800,0,by=-100))) +
  geom_line(data=globDat, aes(x=-bins, y=globMean)) +
  geom_point(data=globDat, aes(x=-bins, y=globMean)) +
  geom_point(data=ints, aes(x=-minAge, y=minT), colour='deepskyblue', size=2) +
  geom_point(data=ints, aes(x=-maxAge, y=maxT), colour='firebrick2', size=2)

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
intPairs$type <- factor(intPairs$type, levels = c('cold-cold','warm-warm','cold-warm'))
colr <- c('cold-cold'="deepskyblue",
          'warm-warm'="firebrick2",
          'cold-warm'="purple3")

# calculate the delta MAT for given bin pairs
globDiff <- function(pair){
  nm1 <- paste0('X', pair['t1'])
  nm2 <- paste0('X', pair['t2'])
  abs(globMean[nm1] - globMean[nm2])
}
intPairs$deltaMAT <- apply(intPairs[,c('t1','t2')], 1, globDiff)
mxDelta <- max(intPairs$deltaMAT) * 1.1
deltaBoxs <- 
  ggplot(data=intPairs) +
  theme_bw() +
  scale_y_continuous('abs diff in global MAT',
                     limits=c(0,mxDelta), expand=c(0,0)
                     ) +
  geom_boxplot(aes(x=type, y=deltaMAT, fill=type)) +
  theme(axis.title.x=element_blank(),
        legend.position='none') +
  scale_fill_manual(values=colr)

# * * KDE overlap for given bin pairs -------------------------------------
source('GSA_custom_ecospat_fcns.R')

nicher <- function(dat, b1, b2, s, env, xmn, xmx, w1, w2, ...){
  sp1rows <- which(dat$species==s & dat$bin==b1)
  sp1 <- dat[sp1rows,env]
  sp2rows <- which(dat$species==s & dat$bin==b2)
  sp2 <- dat[sp2rows,env]
  d1 <- JonesEst(sp1, w = w1, reflect = TRUE, a = xmn, b = xmx, ...) # bw='brt', 
  d2 <- JonesEst(sp2, w = w2, reflect = TRUE, a = xmn, b = xmx, ...) # bw='brt', 
  hell(d1, d2) 
}

# make an outer nested function, to find the survivors of the given time bin pair
# run nicher inside. Take the mean H among survivors; repeat for all time bin pairs.
xmn <- min(samp[,v])
xmx <- max(samp[,v])
bPairs <- as.matrix(intPairs[,c('t1','t2')])

intPairs$avgH <- apply(bPairs, 1, function(x){
  # estimate bias function for each time bin
  # based on sampling distribution, with boundary reflection
  sampRows1 <- which(samp[,'b']==x[1])
  sampRows2 <- which(samp[,'b']==x[2])
  samp1 <- samp[sampRows1,v]
  samp2 <- samp[sampRows2,v]
  # TODO consider making this NOT reflected
  densSamp1 <- density.reflected(samp1, lower=xmn, upper=xmx) 
  densSamp2 <- density.reflected(samp2, lower=xmn, upper=xmx) 
  w1 <- approxfun(densSamp1$x, densSamp1$y)
  w2 <- approxfun(densSamp2$x, densSamp2$y)
  
  # Adjust for any species observations that fall just slightly
  # outside the range of sample density estimation
  # (because of the discretisation of the kernel estimation).
  # Especially a problem for 4ka, upper extreme.
  estMx <- min(max(densSamp1$x), max(densSamp2$x))
  estMn <- max(min(densSamp1$x), min(densSamp2$x))
  pairBool <- samp[,'b'] %in% x
  adj <- samp[pairBool,]
  tooLow <- adj[,v] < estMn
  adj[tooLow,v] <- estMn
  tooHot <- adj[,v] > estMx
  adj[tooHot,v] <- estMx
  
  # find species that are sufficiently sampled in both bins
  spp1 <- df$sp[df$bin==x[1]]
  spp2 <- df$sp[df$bin==x[2]]
  surv <- intersect(spp1, spp2)
  
  # TODO find out why this breaks
  sppH <- sapply(surv, function(s){
    nicher(dat=adj, b1=x[1], b2=x[2], s=s, env=v, 
           w1=w1, w2=w2, xmn=xmn, xmx=xmx)
  })
  mean(sppH)
})

# boxplot of mean H overlap per pairwise comparison category
mxH <- max(intPairs$avgH) * 1.1
ovpBoxs <- ggplot(data=intPairs) +
  theme_bw() +
  scale_y_continuous('mean per-sp, per-interval H',
                     limits=c(0,mxH), expand=c(0,0)
                     ) +
  geom_boxplot(aes(x=type, y=avgH, fill=type)) +
  scale_fill_manual(values=colr) +
  theme(axis.title.x=element_blank(),
        legend.position='none')

# anova(lm(avgH ~ type, data=intPairs))

# * * Combine plot panels -------------------------------------------------

# align time series with MAT plot on left, since y-axes differ in digits
pAlignd <- align_plots(deltaBoxs, globTseries, align = 'v', axis = 'l')

mult_row <- plot_grid(
  pAlignd[[1]], ovpBoxs,
  ncol = 2,
  rel_widths = c(1,1),
  labels=c('B','C')
)
panels <- plot_grid(
  pAlignd[[2]], mult_row,
  labels = c('A', ''),
  ncol=1
  #  rel_heights = c(1, 0.5)
)

if (doTrunc){
  panelsNm <- paste0('Figs/EES-extremes-panels_trunc_',vShrt,'_',day,'.pdf')
} else {
  panelsNm <- paste0('Figs/EES-extremes-panels_',vShrt,'_',day,'.pdf')
}
pdf(panelsNm, width=8, height=8)
panels
dev.off()

# Global MAT vs. species optima -------------------------------------------

optCor <- function(s, var, dat){
  spRows <- grep(s, dat$sp)
  sDat <- dat[spRows,]
  bNms <- paste0('X',sDat$bin)
  gDat <- globMean[bNms]
  cor(gDat, sDat[,var], method='spear')
}

corsM <- sapply(spp, optCor, dat=df, var='m_temp_ym_0m')
corsPe <- sapply(spp, optCor, dat=df, var='pe')
trt <- c(rep('mean',nspp), rep('pref env', nspp))
cors <- data.frame(sp=spp, cor=c(corsM, corsPe), metric=trt)
boxes <- 
  ggplot(data=cors) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_boxplot(aes(x=trt, y=cor)) +
  scale_y_continuous(name='rho corr., sp vs. global MAT', 
                     limits=c(-1,1), expand=c(0,0)) + # c(-0.7,0.7)
  theme(axis.title.x=element_blank())

if (doTrunc){
  boxesNm <- paste0('Figs/global-MAT-corr-w-sp_trunc_',vShrt,'_ka_',day,'.pdf')
} else {
  boxesNm <- paste0('Figs/global-MAT-corr-w-sp_',vShrt,'_ka_',day,'.pdf')
}
pdf(boxesNm, width=4, height=4)
print(boxes)
dev.off()
  
#  t.test(corsM, corsPe, paired = TRUE)

# Time series -------------------------------------------------------------

# Missing sp-bin combinations have already been removed from the df object.
# Add them back in so that the time series are plotted with gaps.

combos <- expand.grid(spp, bins)
colnames(combos) <- c('sp','bin')
df <- merge(combos, df, all.x=TRUE, by=c('sp','bin'))

mCurr <- paste('m', v, sep='_')
sdCurr <- paste('sd', v, sep='_')
df$m <- df[,mCurr]
df$se <- df[,sdCurr]/sqrt(df$n)

spSort <- sort(spp)
df$sp <- factor(df$sp, levels=spSort)

tsPlot <- 
  ggplot(data=df)+
  theme_bw() +
  scale_x_continuous(name='Time (8ka intervals)', expand=c(0.01,0)) +
  geom_line(aes(x=-bin, y=m), size=0.5) + 
  geom_point(aes(x=-bin, y=m), size=0.7) + 
  geom_linerange(aes(x=-bin, ymin=m-se, ymax=m+se), size=0.4)+ 
  facet_wrap(~sp) +
  theme(strip.text.x = element_text(size = 6)) 

if (doTrunc){
  tsNm <- paste0('Figs/time_series_species_mean.se_trunc_', vShrt, day, '.pdf')
} else {
  tsNm <- paste0('Figs/time_series_species_mean.se_', vShrt, day, '.pdf')
}
pdf(tsNm, width=8, height=8)
print(tsPlot)
dev.off()
