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
  df <- read.csv("Data/foram_niche_sumry_metrics_trunc_200214.csv", stringsAsFactors=FALSE)
  allPts <- read.csv('Data/foram_MAT_occs_latlong_8ka_trunc_200213.csv')
  v <- 'temp_ym_0m'
  # See note on niche construction script: if truncating temp by surface values,
  # it doesn't make sense to use in-habitat temperature downstream.
} else {
  df <- read.csv("Data/foram_niche_sumry_metrics_200214.csv", stringsAsFactors=FALSE)
  allPts <- read.csv('Data/foram_MAT_occs_latlong_8ka_200213.csv')
  v <- 'temp_ym_hab'
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

# calculate mean absolute latitude at sampling points
allSampBool <- allPts$species=='sampled'
allSamp <- allPts[allSampBool,]
absLat <- sapply(bins, function(b){
  slcBool <- allSamp$bin==b
  slc <- allSamp[slcBool,]
  abslat <- abs(slc$centroid_lat)
  mean(abslat)
})

# calculate mean MAT over globe, at set grid of points
glob <- read.csv('Data/global_surface_MAT_at_grid_pts_4ka.csv')
cols <- paste0('X',bins)
globMean <- colMeans(glob[,cols])

# combine with sampling n and mean MAT
sampBool <- nich$sp=='sampled1'
samp <- nich[sampBool,]
sampMean <- samp$m_temp_ym_0m
logSampN <- log(samp$n)

# no autocorrelation in residuals of main relationship of interest
l <- lm(sampMean~globMean)
r <- resid(l)
acf(r)

corVars <- cbind(logSampN, absLat, sampMean, globMean)

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

# Inter-spp niche overlap -------------------------------------------------

pairH <- read.csv('Data/foram_species_pairs_KDE_H_200214.csv', stringsAsFactors=FALSE)
# Watch out - not normally distributed because of bounds at 0 and 1.

pairH$bin <- factor(pairH$bin, levels = bins)

# ggplot calls the variables directly - cannot give names as characters
# so duplicate the relevant enviro variable to a column with the consistent name 'y'
yNm <- paste('h', v, sep='_')
pairH$y <- pairH[,yNm]

inter <- 
  ggplot(data=pairH, aes(x=bin, y=y)) + 
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Intra-sp niche overlap --------------------------------------------------

keep <- df$sp != 'sampled'
intraH <- df[keep,]
intraH$y <- intraH[,yNm]
consec <- ! is.na(intraH$y)
intraH <- intraH[consec,]
intraH$shortNm <- sapply(intraH$sp, function(txt){
  splt <- strsplit(txt, ' ')
  splt[[1]][2]
} )

intra <- 
  ggplot(data=intraH, aes(x=shortNm, y=y)) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

y.grob <- textGrob('Hellinger\'s H distance', gp=gpar(fontface='bold', fontsize=15), rot=90)
doubl <- plot_grid(inter, intra, nrow=1, align='h', rel_widths = c(1,0.35))

if (doTrunc){
  ovrlpNm <- paste0('Figs/overlap_H_boxplots_inter_vs_intraspecific_trunc_',vShrt,day,'.pdf')
} else {
  ovrlpNm <- paste0('Figs/overlap_H_boxplots_inter_vs_intraspecific_',vShrt,day,'.pdf')
}
pdf(ovrlpNm, width=9, height=5)
grid.arrange(arrangeGrob(doubl, left = y.grob))
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

empt <- ggplot(data=plotDf, aes(y=y)) +
  theme_bw() +
  scale_y_continuous(limits=c(0,0.3), expand=c(0,0)) +
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
# note that H overlap does not indicate directly, only magnitude, of niche change
# so should compare it with absolute differences in global MAT
for (cutoff in c(156, 800)){
  sumH <- function(b){
    bBool <- intraH$bin==b
    slc <- intraH[bBool,]
    avgH <- mean(slc[,yNm])
    binN <- nrow(slc)
    c(bin=b, avgH=avgH, nSpp=binN)
  }
  Hseq <- sapply(bins, sumH)
  Hseq <- data.frame(t(Hseq))
  # all values NA at most recent time step
  Hseq <- Hseq[-nrow(Hseq),]
  
  delta <- diff(globMean)
  Hseq$absDelta <- abs(delta)
  
  # go back in time only as far as cutoff
  seqKeep <- Hseq$bin <= cutoff
  Hseq <- Hseq[seqKeep,]
  
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
}

# Global MAT vs. species optima -------------------------------------------

for (cutoff in c(156, 800)){
  dfKeep <- df$bin <= cutoff
  dfCut <- df[dfKeep,]
  realSpp <- setdiff(spp, 'sampled')
  nReal <- length(realSpp)
  
  optCor <- function(s, var, dat){
    spRows <- grep(s, dat$sp)
    sDat <- dat[spRows,]
    bNms <- paste0('X',sDat$bin)
    gDat <- globMean[bNms]
    mCol <- paste(var,v,sep='_')
    cor(gDat, sDat[,mCol], method='spear')
  }
  
  corsM <- sapply(realSpp, optCor, dat=dfCut, var='m')
  corsPe <- sapply(realSpp, optCor, dat=dfCut, var='pe')
  trt <- c(rep('mean',length(realSpp)), rep('pref env', length(realSpp)))
  cors <- data.frame(sp=realSpp, cor=c(corsM, corsPe), metric=trt)
  boxes <- 
    ggplot(data=cors) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_boxplot(aes(x=trt, y=cor)) +
    scale_y_continuous(name='rho corr., sp vs. global MAT', 
                       limits=c(-1,1), expand=c(0,0)) + # c(-0.7,0.7)
    theme(axis.title.x=element_blank())
  
  if (doTrunc){
    boxesNm <- paste0('Figs/global-MAT-corr-w-sp_trunc_',vShrt,'_to',cutoff,'ka_',day,'.pdf')
  } else {
    boxesNm <- paste0('Figs/global-MAT-corr-w-sp_',vShrt,'_to',cutoff,'ka_',day,'.pdf')
  }
  pdf(boxesNm, width=4, height=4)
  print(boxes)
  dev.off()
  
#  t.test(corsM, corsPe, paired = TRUE)
}

# * Correlation strength vs n ---------------------------------------------

# At lower n, KDE correction for bias has higher MISE, so induced correlation
# with global temperatures could be higher.
sumN <- function(s){
  sBool <- df$sp==s
  nseq <- df$n[sBool]
  c(avg=mean(nseq), med=median(nseq))
}
spN <- sapply(realSpp, sumN)

peBool <- cors$metric=='pref env'
peCors <- cors[peBool,]
peCors$avgN <- spN['avg',]
peCors$medN <- spN['med',]

medPlot <- 
  ggplot(data=peCors) +
  theme_bw() +
  scale_x_continuous('median N') +
  scale_y_continuous('correlation with global MAT') +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_point(aes(x=medN, y=cor)) 
avgPlot <- 
  ggplot(data=peCors) +
  theme_bw() +
  scale_x_continuous('mean N') +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_point(aes(x=avgN, y=cor)) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())
  
if (doTrunc){
  pairNm <- paste0('Figs/corr_w_global_vs_sample_size_trunc_', vShrt, day, '.pdf')
} else {
  pairNm <- paste0('Figs/corr_w_global_vs_sample_size_', vShrt, day, '.pdf')
}
pdf(pairNm, width=6, height=4)
plot_grid(medPlot, avgPlot, nrow=1, rel_widths = c(1, 0.85))
dev.off()

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

spSort <- c('sampled',sort(realSpp))
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
