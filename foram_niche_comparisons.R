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
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

#envVars <- c('temp_ym_hab','temp_ym_0m')

if (doTrunc){
  df <- read.csv("Data/foram_niche_sumry_metrics_trunc_200214.csv", stringsAsFactors=FALSE)
  v <- 'temp_ym_0m'
  # See note on niche construction script: if truncating temp by surface values,
  # it doesn't make sense to use in-habitat temperature downstream.
} else {
  df <- read.csv("Data/foram_niche_sumry_metrics_200214.csv", stringsAsFactors=FALSE)
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

vRows <- grep(v, mods$var)
modsV <- mods[vRows,]

evoModes <- c('StrictStasis','Stasis','URW','GRW','Punc-1',
              'Stasis-URW','Stasis-GRW','URW-Stasis','GRW-Stasis')
xy <- expand.grid(spMode=evoModes, sampMode=evoModes, stringsAsFactors=FALSE)
xy$n <- NA
for (i in 1:nrow(xy)){
  x <- xy$sampMode[i]
  y <- xy$spMode[i]
  same <- which(modsV$samplingMod==x & modsV$bestMod==y)
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
  geom_point() +
  scale_x_discrete(name = 'Sampling model', drop=FALSE) +
  scale_y_discrete(name = 'Species model', drop=FALSE) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_size(range=c(2,10), breaks = c(1,3,5,7,9))

bubblNm <- paste0('Figs/evo_mode_bubble_matrix_', vShrt, day, '.pdf')
pdf(bubblNm, width=5, height=4.4)
print(bubbl)
dev.off()

# Inter- vs intra- spp overlap --------------------------------------------

# * Inter-spp niche overlap -----------------------------------------------

pairH <- read.csv('Data/foram_species_pairs_KDE_H_200213.csv', stringsAsFactors=FALSE)
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

# * Intra-sp niche overlap ------------------------------------------------

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

ovrlpNm <- paste0('Figs/overlap_H_boxplots_inter_vs_intraspecific_',vShrt,day,'.pdf')
pdf(ovrlpNm, width=9, height=5)
grid.arrange(arrangeGrob(doubl, left = y.grob))
dev.off()

# * * Intra overlap vs evo mode -------------------------------------------

realSpp <- setdiff(spp, 'sampled')

# arrange boxplot in order of inceasing niche lability/larger H
meanH <- function(s,dat){
  sBool <- dat$sp==s
  sDat <- dat[sBool,]
  hNm <- paste('h',v,sep='_')
  mean(sDat[,hNm])
}
hVect <- sapply(realSpp, meanH, dat=intraH)
spOrdr <- names(sort(hVect))
intraH$sp <- factor(intraH$sp, levels=spOrdr)

# add colour-codes for evo model type:
# strict and loose stasis vs. non-stasis vs. mix
# note: no species were fit with punct-1, so this isn't classified
stBool <- modsV$bestMod %in% c('StrictStasis','Stasis')
stSp <- modsV$sp[stBool]
labBool <- modsV$bestMod %in% c('GRW','URW')
labSp <- modsV$sp[labBool]
mixBool <- modsV$bestMod %in% 
  c('GRW-Stasis','URW-Stasis','Stasis-GRW','Stasis-URW')
mixSp <- modsV$sp[mixBool]

for (s in realSpp){
  sBool <- intraH$sp==s
  
  isSt <- length(grep(s, stSp))
  isLab <- length(grep(s, labSp))
  isMix <- length(grep(s, mixSp))
  # check for cases where a species has multiple segments
  # and different model types are fit to them
  tests <- c(isSt!=0, isLab!=0, isMix!=0)
  if (sum(tests) > 1){
    intraH$class[sBool] <- 'NA'
  } else {
    if (isSt > 0){
      intraH$class[sBool] <- 'Static '
    }
    if (isLab > 0){
      intraH$class[sBool] <- 'Labile '
    }
    if (isMix > 0){
      intraH$class[sBool] <- 'Complex '
    }
  }
  if ((sum(tests) == 0)) stop(paste('no classification for', s))
}

colr <- c('#3333FF','#FFFF33','#00CC00','#FFFFFF')
names(colr) <- c('Static ','Labile ','Complex ','NA')

evoHplot <- 
  ggplot(data=intraH, aes(x=sp, y=y)) +
  theme_bw() +
  scale_y_continuous(name='Hellinger\'s H', limits=c(0,1), expand=c(0,0)) + 
  geom_boxplot(aes(fill=class)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(name='Evo model class', values=colr, 
                    breaks=c('Static ','Labile ','Complex ','NA')) 

evoHnm <- paste0('Figs/overlap_H_boxplots_intraspecific_vs_evo_mode_',vShrt,day,'.pdf')  
pdf(evoHnm, width=7, height=4)
print(evoHplot)
dev.off()

# * * Intra overlap vs evo var --------------------------------------------

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

vhNm <- paste0('Figs/overlap_H_vs_evo_mod_variance_',vShrt,day,'.pdf')
pdf(vhNm, width=5, height=4)
vh
dev.off()

# * * Overlap by ecomorph -------------------------------------------------

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

multNm <- paste0('Figs/overlap_H_by_ecomorph_',vShrt,day,'.pdf')
pdf(multNm, width=4, height=7)
grid.arrange(arrangeGrob(smallMult), bottom='Mean Hellinger\'s H among species') 
dev.off()

# Sampled vs. species optima ----------------------------------------------

nReal <- length(realSpp)

optCor <- function(s, var, dat){
  spRows <- grep(s, dat$sp)
  sDat <- dat[spRows,]
  bNms <- paste0('X',sDat$bin)
  gDat <- globMean[bNms]
  mCol <- paste(var,v,sep='_')
  cor(gDat, sDat[,mCol], method='spear')
}

corsM <- sapply(realSpp, optCor, dat=df, var='m')
corsPe <- sapply(realSpp, optCor, dat=df, var='pe')
trt <- c(rep('mean',length(realSpp)), rep('pref env', length(realSpp)))
cors <- data.frame(cor=c(corsM, corsPe), metric=trt)
yNm <- paste('rho corr., sp',metric,' vs. global MAT')
boxes <- 
  ggplot(data=cors) +
  geom_boxplot(aes(x=trt, y=cor)) +
  scale_y_continuous(name=yNm, limits=c(-1,1), expand=c(0,0)) + # c(-0.7,0.7)
  theme(axis.title.x=element_blank())
boxesNm <- paste0('Figs/global_MAT_corr_w_sp_', metric, '_', vShrt, day, '.pdf')
pdf(boxesNm, width=4, height=4)
print(boxes)
dev.off()

# Time series -------------------------------------------------------------

# Missing sp-bin combinations have already been removed from the df object.
# Add them back in so that the time series are plotted with gaps.

combos <- expand.grid(spp, bins)
colnames(combos) <- c('sp','bin')
df <- merge(combos, df, all.x=TRUE, by=c('sp','bin'))

mCurr <- paste('m', v, sep='_')
sdCurr <- paste('sd', v, sep='_')
nCurr <- paste('n', v, sep='_')
df$m <- df[,mCurr]
df$se <- df[,sdCurr]/sqrt(df[,nCurr])

spSort <- c('sampled',sort(realSpp))
df$sp <- factor(df$sp, levels=spSort)

# this chunk draws on objects from 'spp vs sampling evo mode' section
for (s in spp){
  sBool <- df$sp==s
  
  isSt <- length(grep(s, stSp))
  isLab <- length(grep(s, labSp))
  isMix <- length(grep(s, mixSp))
  # check for cases where a species has multiple segments
  # and different model types are fit to them
  tests <- c(isSt!=0, isLab!=0, isMix!=0)
  if (sum(tests) > 1){
    df$class[sBool] <- 'Complex '
  } else {
    if (isSt > 0){
      df$class[sBool] <- 'Static '
    }
    if (isLab > 0){
      df$class[sBool] <- 'Labile '
    }
    if (isMix > 0){
      df$class[sBool] <- 'Complex '
    }
  }
  if ((sum(tests) == 0)) stop(paste('no classification for', s))
}

colr <- c('#3333FF','#FFFF33','#00CC00')
names(colr) <- c('Static ','Labile ','Complex ')

tsPlot <- ggplot(data=df)+
  theme_bw() +
  scale_x_continuous(name='Time (8ka intervals)', expand=c(0.01,0)) +
  geom_line(aes(x=-bin, y=m, col=class), size=0.5) +
  geom_point(aes(x=-bin, y=m, col=class), size=0.7) +
  geom_linerange(aes(x=-bin, ymin=m-se, ymax=m+se, col=class), size=0.4)+
  facet_wrap(~sp) +
  theme(strip.text.x = element_text(size = 6)) +
  theme(legend.position = 'bottom') +
  scale_fill_manual(name='Evo model class', values=colr, 
                    breaks=c('Static ','Labile ','Complex ')) 

tsNm <- paste0('Figs/time_series_species_mean.se_', vShrt, day, '.pdf')
pdf(tsNm, width=8, height=9)
tsPlot
dev.off()
