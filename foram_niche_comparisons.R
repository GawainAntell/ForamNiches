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

# Data prep ---------------------------------------------------------------
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

df <- read.csv("Data/foram_niche_sumry_metrics_raw_values_200116.csv", stringsAsFactors=FALSE)
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]

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
allPts <- read.csv('Data/foram_MAT_occs_latlong_8ka_200116.csv')
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
sampMean <- samp$m
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

corNm <- paste0('Figs/correlations_sampled_vs_global_',day,'.pdf')
pdf(corNm, width=5, height=5)
chart.Correlation(corVars, histogram=TRUE, pch=19, method='pearson')
dev.off()

# Evo models --------------------------------------------------------------

evoFit <- function(s){
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
  
  ts <- as.paleoTS(mm = sp$m, vv = sp$sd^2, nn = sp$n, tt = sp$scaledT, 
                   oldest = 'first', reset.time = FALSE)
  tsSamp <- as.paleoTS(mm = samp$m, vv = samp$sd^2, nn = samp$n, tt = samp$scaledT, 
                       oldest = 'first', reset.time = FALSE)
  
  if (l < 14){
    mods <- fit4models(ts, method='Joint', silent=TRUE)
    modsSamp <- fit4models(tsSamp, method='Joint', silent=TRUE)
  } else {
    mods <- fit9models(ts, method='Joint', silent=TRUE)
    modsSamp <- fit9models(tsSamp, method='Joint', silent=TRUE)
  }
  # From the package documentation:
  # 'Method = "Joint" is a full likelihood approach, considering each time-series as a joint sample
  # from a multivariate normal distribution. Method = "AD" is a REML approach that uses the differences
  # between successive samples. They perform similarly, but the Joint approach does better under
  # some circumstances (Hunt, 2008).'
  # From Hunt 2008:
  # 'joint parameterization  is  better  able  to  correctly  identify  directional  trends. 
  # This advantage increases with sequence length, and is most pronounced when sampling error is high'
  
  wts <- mods$modelFits$Akaike.wt
#  mods$modelFits[order(wts),]
  maxMod <- which.max(wts)
  modNm <- row.names(mods$modelFits)[maxMod]
  w <- max(wts)
  params <- mods$parameters[modNm][[1]]
  
  # check what model would be predicted by sampling alone for the given time bins
  sampMx <- which.max(modsSamp$modelFits$Akaike.wt)
  sampEvo <- row.names(modsSamp$modelFits)[sampMx]
  
  out <- data.frame(sp=s, bestMod=modNm, weight=w, l=l, start=strt, samplingMod=sampEvo)
  out$params <- list(params)
  out
}

# save package names to put  on all cores
pkgs <- c('paleoTS') 

pt1 <- proc.time()
ncores <- detectCores() - 1
registerDoParallel(ncores)
  mods <- foreach(s=longSpp, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar% evoFit(s)
stopImplicitCluster()
pt2 <- proc.time()
pt2-pt1
# only 2 min runtime

table(mods$bestMod)

# * Sampling model --------------------------------------------------------

# ages must start at 0
samp$scaledT <- 1:nrow(samp) -1
ts <- as.paleoTS(mm = samp$m, vv = samp$sd^2, nn = samp$n, tt = samp$scaledT, 
                 oldest = 'first', reset.time = FALSE)
sampMods <- fit9models(ts, method='Joint', silent=TRUE)
sampMods

# * Plot sampled MAT t-series ---------------------------------------------

samp$se <- samp$sd/sqrt(samp$n)
puncParams <- sampMods$parameters$`Punc-1`
jump <- puncParams['shift1']
y1 <- puncParams['theta1']
y2 <- puncParams['theta2']
linep <- ggplot(data=samp) +
  theme_bw() +
  scale_y_continuous(name='mean annual temp at surface') +
  scale_x_continuous(expand=c(0.01,0), name='ka')+
  geom_linerange(aes(x=-bin, ymax=m+se, ymin=m-se)) +
  geom_line(aes(x=-bin, y=m)) +
  geom_segment(x=-bins[1], xend=-bins[jump], y=y1, yend=y1, colour='blue') +
  geom_segment(x=0, xend=-bins[jump], y=y2, yend=y2, colour='blue') 

lineNm <- paste0('Figs/sampled_MAT_tseries_wPunc_',day,'.pdf')
pdf(lineNm, width=6, height=4)
linep
dev.off()

# * Inspection of mods by evo type ----------------------------------------
puncBool <- mods$bestMod=='Punc-1'
punc <- mods[puncBool,]
punc$params

grwBool <- mods$bestMod=='GRW'
grw <- mods[grwBool,]
grw$params

stasisBool <- mods$bestMod=='StrictStasis'
stasis <- mods[stasisBool,]

spModBool <- mods$sp != 'sampled1'
spMods <- mods[spModBool,]
spMods$bestMod <- factor(spMods$bestMod, levels=c('StrictStasis','GRW','Punc-1','Stasis-GRW'))
ymx <- max(table(spMods$bestMod))*1.1
bp <- ggplot(data=spMods, aes(bestMod)) + 
  theme_minimal() + #theme_bw() +
  scale_y_continuous(name='n species sequences', expand=c(0,0), limits = c(0,ymx), labels=NULL) +
  scale_x_discrete(name='Evolutionary model', expand=c(0,0)) +
  geom_bar() +
  geom_text(aes(label=..count..),stat="count",position=position_stack(.7),colour='white')
bpNm <- paste0('Figs/evo_mode_barplot_',day,'.pdf')
pdf(bpNm, width=4, height=3)
bp
dev.off()

# Inter- vs intra- spp overlap --------------------------------------------

# * Inter-spp niche overlap -----------------------------------------------

pairD <- read.csv('Data/foram_species_pairs_KDE_D_200116.csv', stringsAsFactors=FALSE)
# Watch out - not normally distributed because of bounds at 0 and 1.
# But some values = 1 do occur, so can't do logit transformation.

pairD$bin <- factor(pairD$bin, levels = bins)
inter <- 
  ggplot(data=pairD, aes(x=bin, y=d)) +
  scale_y_continuous(limits=c(0,1.05), expand=c(0,0)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# * Intra-sp niche overlap ------------------------------------------------

kdeFull <- read.csv("Data/foram_niche_sumry_metrics_KDE_200116.csv", stringsAsFactors=FALSE)
kdeKeep <- kdeFull$sp != 'sampled'
kde <- kdeFull[kdeKeep,]
consec <- ! is.na(kde$d)
kde <- kde[consec,]
kde$shortNm <- sapply(kde$sp, function(txt){
  splt <- strsplit(txt, ' ')
  splt[[1]][2]
} )

intra <- 
  ggplot(data=kde, aes(x=shortNm, y=d)) +
  scale_y_continuous(limits=c(0,1.05), expand=c(0,0)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

y.grob <- textGrob('Schoener\'s D overlap', gp=gpar(fontface='bold', fontsize=15), rot=90)
doubl <- plot_grid(inter, intra, nrow=1, align='h', rel_widths = c(1,0.35))

ovrlpNm <- paste0('Figs/overlap_D_boxplots_inter_vs_intraspecific',day,'.pdf')
pdf(ovrlpNm, width=9, height=5)
grid.arrange(arrangeGrob(doubl, left = y.grob))
dev.off()

# Sampled vs. species optima ----------------------------------------------
realSpp <- setdiff(spp, 'sampled')
optCor <- function(s, dat){
  spBool <- dat$sp==s
  sDat <- dat[spBool,]
  bNms <- paste0('X',sDat$bin)
  gDat <- globMean[bNms]
  cor(gDat, sDat$pe, method='spear')
}
cors <- sapply(realSpp, optCor, dat=kdeFull)
summary(cors)
sum(cors>0)/length(cors)
boxplot(cors)
