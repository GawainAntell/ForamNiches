#library(lme4)
library(paleoTS)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
#library(VoCC)
#devtools::install_github("JorGarMol/VoCC", dependencies = FALSE, build_vignettes = FALSE)

# Data prep ---------------------------------------------------------------

df <- #read.csv("Data/foram_niche_sumry_metrics_KDE_200106.csv", stringsAsFactors=FALSE)
       read.csv("Data/foram_niche_sumry_metrics_raw_values_200106.csv", stringsAsFactors=FALSE)
ordr <- order(df$bin, decreasing = TRUE)
df <- df[ordr,]

# if using KDE data, rename niche variables for consistency with those from raw values:
nc <- ncol(df)
if (nc > 5){
  colnames(df)[(nc-1) : nc] <- c('m','v')
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
  # ages must start at 0
  sp$scaledT <- 1:nrow(sp) -1
  
  # save metadata about the time series
  l <- nrow(sp)
  strt <- sp$bin[1]
  
  ts <- as.paleoTS(mm = sp$m, vv = sp$v, nn = sp$n, tt = sp$scaledT, 
                   oldest = 'first', reset.time = FALSE)
  
  if (l < 14){
    mods <- fit4models(ts, method='Joint', silent=TRUE)
  } else {
    mods <- fit9models(ts, method='Joint', silent=TRUE)
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
  
  out <- data.frame(sp=s, bestMod=modNm, weight=w, l=l, start=strt)
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

# Punctuated mods ---------------------------------------------------------
puncBool <- mods$bestMod=='Punc-1'
punc <- mods[puncBool,]
punc$params

grwBool <- mods$bestMod=='GRW'
grw <- mods[grwBool,]
grw$params
