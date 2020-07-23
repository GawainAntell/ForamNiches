library(paleoTS)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(tidyr)

# set whether to run analyses on surface-level or in-habitat niches
ss <- TRUE

# Data prep ---------------------------------------------------------------
day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")

# foram data
spAttr <- read.csv('Data/foram-spp-data_2020-07-21.csv', stringsAsFactors=FALSE)
if (ss){
  df <- read.csv('Data/niche-sumry-metrics_nrd0_SS_2020-07-22.csv', stringsAsFactors=FALSE)
} else {
  df <- read.csv('Data/niche-sumry-metrics_nrd0_hab_2020-07-22.csv', stringsAsFactors=FALSE)
}
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]
bins <- unique(df$bin)
spp <- unique(df$sp)
nspp <- length(spp)
binL <- bins[1] - bins[2]
envNm <- 'temp_ym'

# standardized sampling universe (MAT at core sites) at each of 4 depths
dList <- readRDS('Data/spp-and-sampling-data_list-by-depth_2020-07-21.rds')
# get mean, sd, and n across sites in order to fit evo models to sampling seqs
bsum <- function(d){
  bMat <- sapply(bins, function(b){
    bBool <- d$samp$bin==b
    slc <- d$samp[bBool,]
    m <- mean(slc$temp_ym)
    sdv <- sd(slc$temp_ym)
    nsite <- nrow(slc)
    cbind(b, m, sdv, nsite)
  })
  bDf <- data.frame(t(bMat))
  colnames(bDf) <- c('bin','m','sd','nsite')
  bDf
}
samp <- lapply(dList, bsum)

# Deal with incomplete sp ts ----------------------------------------------

# rename species so non-consecutive sequences have a numerical suffix
spL <-  vector('list', length = nspp)
for (i in 1:nspp){
  spNm <- spp[i]
  spBool <- df$sp==spNm
  sp <- df[spBool,]
  j <- 1
  sp$sp[1] <- paste0(spNm, j)
  penult <- nrow(sp) - 1
  for (r in 1:penult){
    age <- sp$bin[r]
    nextAge <- sp$bin[r+1]
    if (age != nextAge + binL){
      j <- j + 1
    }
    sp$sp[r+1] <- paste0(spNm, j)
  }
  spL[[i]] <- sp
}
splitSpp <- do.call(rbind, spL)

# Check for per-species continuity - ditch any lineage chuncks of < 7 time bins
longSpp <- character()
enuf <- rep(binL, 6)
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
tsDf <- splitSpp[keepBool,]

# Evo models --------------------------------------------------------------

# Convert an environmental variable name to habitat name
trnslt <- function(txt){ 
  switch(txt, 
         Surface='temp_ym_0m',
         Surface.subsurface='temp_ym_surfsub',
         Subsurface='temp_ym_sub'
  )
}

evoFit <- function(s, dat, sampDat, spDat, nmods='four', method='AD'){
  spBool <- dat$sp==s
  sp <- dat[spBool,]
  # also subset the sampling time series to the same bins, for comparison
  sBins <- sp$bin
  
  # lookup the species' habitat and use corresponding depth for sampling data
  sShort <- substr(s, 1, (nchar(s)-2))
  attrRow <- grep(paste0(sShort,'*'), spDat$species)
  hab <- spDat$DepthHabitat[attrRow]
  lvl <- trnslt(hab)
  sampDpth <- sampDat[[lvl]]
  sampBool <- sampDpth$bin %in% sBins
  sampSame <- sampDpth[sampBool,]
  
  # correlate species response and global MAT.
  # too much autocorrelation in residuals of mean vs. MAT so use KDE optimum
  # peLm <- lm(sp$pe ~ sampSame$m)
  # acf(peLm$residuals)
  sCor <- cor(sp$pe, sampSame$m, method='pear')
  
  # ages must start at 0
  sampSame$scaledT <- sp$scaledT <- 1:nrow(sp) -1
  
  # save metadata about the time series
  l <- nrow(sp)
  strt <- sp$bin[1]
  
  ts <- as.paleoTS(mm = sp$m, vv = sp$sd^2, 
                   nn = sp$n, tt = sp$scaledT, 
                   oldest = 'first', reset.time = FALSE)
  tsSamp <- as.paleoTS(mm = sampSame$m, 
                       vv = sampSame$sd^2, 
                       nn = sampSame$nsite, 
                       tt = sampSame$scaledT, 
                       oldest = 'first', reset.time = FALSE)
  
  # some sample variances are not equal, so use pool=FALSE for all sequences
  if (l > 13 & nmods=='nine'){
    modsSp <- fit9models(ts, method=method, silent=TRUE, pool=FALSE)
    modsSamp <- fit9models(tsSamp, method=method, silent=TRUE, pool=FALSE)
  } else {
    modsSp <- fit4models(ts, method=method, silent=TRUE, pool=FALSE) 
    modsSamp <- fit4models(tsSamp, method=method, silent=TRUE, pool=FALSE) 
  }
  # From the package documentation:
  # 'Method = "Joint" is a full likelihood approach, considering each time-series as a joint sample
  # from a multivariate normal distribution. Method = "AD" is a REML approach that uses the differences
  # between successive samples. They perform similarly, but the Joint approach does better under
  # some circumstances (Hunt, 2008).'
  # From Hunt 2008:
  # 'joint parameterization  is  better  able  to  correctly  identify  directional  trends. 
  # This advantage increases with sequence length, and is most pronounced when sampling error is high'
  
  # HOWEVER, although joint might be better for most series, for 'Globigerinella calida6' 
  # a fatal error occurs (matrix not full rank) unless pool and hess arguments are modified
  # Error in optim(p0, fn = logL.joint.URW, control = cl, method = meth, lower = c(NA,  : 
  # non-finite value supplied by optim
  
  modes <- c('GRW','URW','Stasis','StrictStasis')
  wts <- sapply(modes, function(m) modsSp$modelFits[m,'Akaike.wt'])
  modNm <- names(which.max(wts))
  
  # check what model would be predicted by sampling alone for the given time bins
  sampMx <- which.max(modsSamp$modelFits$Akaike.wt)
  sampEvo <- row.names(modsSamp$modelFits)[sampMx]
  
  out <- data.frame(sp = s, l = l, start = strt, r = sCor, 
                    bestMod = modNm, t(wts), samplingMod = sampEvo)
  out
}

mods4l <- lapply(longSpp, evoFit, dat = tsDf, 
                 sampDat = samp, spDat = spAttr, nmods = 'four') 
# longSpp2 <- setdiff(longSpp, 'Globigerinella calida6')
# mods4l <- lapply(longSpp2, evoFit, dat = tsDf, method='Joint',
#                  sampDat = samp, spDat = spAttr, nmods = 'four') 
mods4 <- do.call('rbind', mods4l)

# when G calida6 excluded, using method='Joint', SURFACE
# GRW       Stasis StrictStasis          URW 
# 6            8            7           10
# IN-HABITAT, joint
# GRW       Stasis StrictStasis          URW 
# 2            9            8           12
table(mods4$bestMod)
mean(mods4$r); sd(mods4$r)

# Spp vs sampling evo mode ------------------------------------------------

evoModes <- c('StrictStasis','Stasis','URW','GRW')
# if using 9 models, add:
# 'Punc-1','Stasis-URW','Stasis-GRW','URW-Stasis','GRW-Stasis'

xy <- expand.grid(spMode=evoModes, sampMode=evoModes, stringsAsFactors=FALSE)
xy$n <- NA
for (i in 1:nrow(xy)){
  x <- xy$sampMode[i]
  y <- xy$spMode[i]
  same <- which(mods4$samplingMod==x & mods4$bestMod==y)
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

if (ss){
  bubblNm <- paste0('Figs/evo-mode-bubble-matrix_SS_', day, '.pdf')
} else {
  bubblNm <- paste0('Figs/evo-mode-bubble-matrix_hab_', day, '.pdf')
}
pdf(bubblNm, width=4, height=4)
print(bubbl)
dev.off()

# Spp time series ---------------------------------------------------------

# Missing sp-bin combinations have already been removed from the df object.
# Add them back in so that the time series are plotted with gaps.
combos <- expand.grid(spp, bins)
colnames(combos) <- c('sp', 'bin')
df <- merge(combos, df, all.x = TRUE, by = c('sp', 'bin'))
spSort <- sort(spp)
df$sp <- factor(df$sp, levels = spSort)

xtick <- paste(seq(600, 0, by = -200))
xbreak <- seq(-600, 0, by = 200)
tsPlot <- 
  ggplot(data = df)+
  theme_bw() +
  scale_y_continuous(name = 'Mean annual temperature at occupied sites') +
  scale_x_continuous(name = 'Age (ka)', expand = c(0.01, 0),
                     breaks = xbreak, labels = xtick) +
  geom_ribbon(aes(x = -bin, ymin = m-sd, ymax = m+sd), 
              fill = 'grey50', alpha = 0.5, size = 0.4) + 
  geom_line(aes(x = -bin, y = m), size = 0.5) + 
  facet_wrap(~sp) +
  theme(strip.text.x = element_text(size = 6)) 

if (ss){
  tsNm <- paste0('Figs/species-time-series_SS_', day, '.pdf')
} else {
  tsNm <- paste0('Figs/species-time-series_hab_', day, '.pdf')
}
pdf(tsNm, width=7.5, height=8)
print(tsPlot)
dev.off()

# Sampling time series ----------------------------------------------------

plotL <- list()
nlvl <- length(samp)
lvlNms <- c(' Sea level',' Surface habitat',
            ' Surface-subsurface habitat',' Subsurface habitat')
for (i in 1:nlvl){
  sampPlot <- 
    ggplot(data=samp[[i]])+
    theme_bw() +
    ggtitle(lvlNms[i]) +
    scale_y_continuous(limits = c(2.5, 27.5), expand = c(0,0)) +
    scale_x_continuous(expand = c(0.01,0), breaks = xbreak, labels = xtick) +
    geom_ribbon(aes(x = -bin, ymin = m-sd, ymax = m+sd), 
                fill = 'grey50', alpha = 0.5, size = 0.4) + 
    geom_line(aes(x = -bin, y = m), size = 0.5) +
    theme(axis.title = element_blank())
  plotL <- append(plotL, list(sampPlot))
}

# export as compound plot
pg <- plot_grid(plotlist = plotL, nrow = 2, labels = 'AUTO', 
                label_size = 14, hjust = -0.7) 
yGrob <- textGrob('Temperature at sample sites', 
                  gp = gpar(fontface="bold"), rot = 90)
xGrob <- textGrob('Age (ka)', 
                   gp = gpar(fontface="bold"))
sampNm <- paste0('Figs/sampling-time-series_by-depth_', day, '.pdf')
pdf(sampNm, width = 6, height = 6)
grid.arrange(arrangeGrob(pg, left = yGrob, bottom = xGrob))
dev.off()

# Relative model support --------------------------------------------------
# stacked bar chart of support for each evo model, for each sequence

# first 3 letters of sp epithet are unique except for ruber and rubescens
# but BEWARE if crassula and crassaformis are later both included
mods4$sp <- gsub('Globoturborotalita rubescens6','Globoturborotalita rbs6',
                 mods4$sp)
mods4$seq <- ''
for (i in 1:nrow(mods4)){
  binom <- mods4$sp[i]
  nmVect <- strsplit(binom, ' ')[[1]]
  abrv <- paste(substr(nmVect[1], 1, 2), 
                 substr(nmVect[2], 1, 3), sep='.')
  # if the species already has sequences, add a distinct suffix
  nrep <- length(grep(abrv, mods4$seq))
  abrvUniq <- paste0(abrv, nrep+1)
  mods4$seq[i] <- abrvUniq
}

abrvOrdr <- order(mods4$seq)
mods4 <- mods4[abrvOrdr,]
modLong <- pivot_longer(mods4, cols = all_of(evoModes), 
                        names_to = 'model', values_to = 'weight')
modLong$model <- factor(modLong$model, levels = rev(evoModes))
modLong$seq <- factor(modLong$seq, levels = rev(unique(modLong$seq)))

colr <- c('StrictStasis'='#283593', 'Stasis'='#5dade2', 
          'URW'='#FFC300', 'GRW'='#FF5733')

bars <- 
  ggplot() +
  theme_minimal() +
  # use position=fill to indicate data as %, so rounding errors don't result in sum > 1  
  geom_bar(data = modLong, aes(x = seq, y = weight, fill = model), 
           position = 'fill', stat = 'identity', width = 0.75) +
  scale_x_discrete(name = 'Species sequence') +
  scale_y_continuous(name = 'Model support (AIC weight)        ', 
                     expand = c(0, 0), limits = c(0, 1.2),
                     breaks = seq(0, 1, by = 0.25)) +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank()) +
  scale_colour_manual(name = element_blank(), values = colr, 
                      aesthetics = 'fill', limits = evoModes,
                      labels = c('Strict stasis','Stasis','Random walk','Directional walk')) +
  guides(fill = guide_legend(nrow = 2))
barsFlip <- 
  bars +
  geom_text(aes(x = seq, y = 1.1, label = l), 
            data = mods4, hjust = 1, size = 3.2) +
  coord_flip()

if (ss){
  barNm <- paste0('Figs/evo-model-support-barplot_SS_',day,'.pdf')
} else {
  barNm <- paste0('Figs/evo-model-support-barplot_hab_',day,'.pdf')
}
# the column width is 3.4252 in, but the exported plot has white space to trim
pdf(barNm, width=3.4252+0.3, height=6)
barsFlip
dev.off()
