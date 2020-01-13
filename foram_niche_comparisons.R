library(paleoTS)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
#library(VoCC)
#devtools::install_github("JorGarMol/VoCC", dependencies = FALSE, build_vignettes = FALSE)

# Data prep ---------------------------------------------------------------
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

df <- read.csv("Data/foram_niche_sumry_metrics_raw_values_200113.csv", stringsAsFactors=FALSE)
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]

# if using KDE data, rename niche variables for consistency with those from raw values:
nc <- ncol(df)
if (nc > 5){
  colnames(df)[(nc-1) : nc] <- c('m','sd')
}

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

# Sampling model ----------------------------------------------------------

sampBool <- nich$sp=='sampled1'
samp <- nich[sampBool,]
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

# Inspection of mods by evo type ------------------------------------------
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

#source('ecospat.grid.clim.dyn.GSA.fcn.R')

df <- read.csv('Data/foram_MAT_occs_latlong_8ka_200113.csv',stringsAsFactors=FALSE)
bins <- unique(df$bin)
nbins <- length(bins)
spp <- unique(df$species)
env <- 'mat'
h.method <- "nrd0" # "SJ-ste" # "ucv"
R <- 1000

interSppD <- function(b,df,R,h.method){
  globBool <- df$species=='sampled'
  glob <- df[globBool,env]
  
  glob1rows <- which(df$species=='sampled' & df$bin==b)
  glob1 <- df[glob1rows,env]
  
  bSppRows <- which(df$species!='sampled' & df$bin==b)
  bSpp <- unique(df$species[bSppRows])
  
  # Construct KDE of all species
  kdeL <- lapply(bSpp, function(s){
    spRows <- which(df$species==s & df$bin==b)
    sp1 <- df[spRows,env]
    z <- GSA.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0, th.env=0, h.method=h.method)
  }
  )
  names(kdeL) <- bSpp
  
  # Compute D for all pairs
  d <- numeric()
  for (s1 in bSpp){
    for (s2 in bSpp){
      if (s1==s2) {next} else{
        ovrlpL <- GSA.ecospat.niche.overlap(kdeL[[s1]], kdeL[[s2]], cor=FALSE)
        d <- c(d, ovrlpL$D)
      }
    }
  }
  
  # Save a summary overlap value among species 
  # Watch out - not normally distributed because of bounds at 0 and 1.
  # But some values = 1 do occur, so can't do logit transformation.
  # Safest approach is to work with the median isntead of mean.
  median(d)
}

# Find the median niche overlap among species pairs in each time bin
tMeds <- sapply(bins, interSppD, df=df, R=R, h.method=h.method)

# * Intra-sp niche overlap ------------------------------------------------

kde <- read.csv("Data/foram_niche_sumry_metrics_KDE_200108.csv", stringsAsFactors=FALSE)
kdeSpp <- unique(kde$sp)
kdeSpp <- setdiff(kdeSpp, 'sampled')
spMeds <- sapply(kdeSpp, function(s){
  sBool <- kde$sp==s
  slc <- kde[sBool,]
  median(slc$d, na.rm=T)
})

summary(tMeds)
summary(spMeds)
t.test(tMeds, spMeds)

typeVect <- c(rep('intra',length(spMeds)), rep('inter',length(tMeds)) )
dDf <- data.frame(d=c(spMeds, tMeds), type=typeVect)
dDf$type <- factor(dDf$type)
table(dDf$type)
ymin <- min(dDf$d)*.97
dBox <- ggplot(data=dDf, aes(x=type, y=d)) +
  theme_minimal() +
  scale_y_continuous(name='Schoener\'s D overlap', limits=c(ymin, 1), expand=c(0,0)) +
  scale_x_discrete(name='', labels=c('Within species','Between species')) +
  geom_boxplot(fill='grey')
boxNm <- paste0('Figs/overlap_inter_vs_intra_medianD_',day,'.pdf')
pdf(boxNm, width=4, height=4)
  dBox
dev.off()

# Sampled vs. species optima ----------------------------------------------

optCor <- function(s, dat){
  sampBool <- dat$sp=='sampled'
  samp <- dat[sampBool,]
  # acf(samp$pe) # surprisingly small AC
  spBool <- dat$sp==s
  sDat <- dat[spBool,]
  sampSame <- samp$bin %in% sDat$bin
  samp <- samp[sampSame,]
  cor(samp$pe, sDat$pe, method='spear')
}
cors <- sapply(kdeSpp, optCor, dat=kde)
summary(cors)
sum(cors>0)/length(cors)

