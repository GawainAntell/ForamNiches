#library(tidyr)
library(lme4)
library(adehabitatMA)
library(adehabitatHR)
library(ecospat)
library(paleoTS)
#library(VoCC)
#devtools::install_github("JorGarMol/VoCC", dependencies = FALSE, build_vignettes = FALSE)

# Data prep ---------------------------------------------------------------

source('ecospat.grid.clim.dyn.GSA.fcn.R')

df <- read.csv('Data/foram_MAT_occs_latlong_8ka_200106.csv',stringsAsFactors=FALSE)
bins <- unique(df$bin)
nbins <- length(bins)
spp <- unique(df$species)

env <- 'mat'
h.method <- c("nrd0","ucv","SJ-ste")[1]
# Resolution of the gridding of the climate space. Ecospat default is 100.
R <- 1000

# Calculate niche overlap (Schoener's D), peak abundance, preferred enviro, & tolerance
nicher <- function(b1, b2, sp, env, h.method){

  globBool <- df$species=='sampled'
  glob <- df[globBool,env]

  glob1rows <- which(df$species=='sampled' & df$bin==b1)
  glob1 <- df[glob1rows,env]

  glob2rows <- which(df$species=='sampled' & df$bin==b2)
  glob2 <- df[glob2rows,env]

  sp1rows <- which(df$species==sp & df$bin==b1)
  sp1 <- df[sp1rows,env]

  sp2rows <- which(df$species==sp & df$bin==b2)
  sp2 <- df[sp2rows,env]
  
  # for each species at time i and i+1
  z1 <- tryCatch(
    GSA.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0, th.env=0, h.method=h.method),
    error = function(err){ list() }
  ) 
  z2 <- tryCatch(
    GSA.grid.clim.dyn(glob, glob2, sp2, R, th.sp=0, th.env=0, h.method=h.method),
    error = function(err){ list() }
  ) 

  # the species may be absent in one or both bins, in which z1|2 is an empty vector
  if (length(z1)==0){
    data.frame(bin=NA, sp=NA, n=NA, d=NA, pa=NA, pe=NA, tol=NA)
  } else {
    n <- length(sp1rows)
    if (length(z2)==0){
      data.frame(bin=b1, sp=sp, n=n, d=NA, pa=z1$pa, pe=z1$pe, tol=z1$t)
    } else{
      ovrlp <- GSA.ecospat.niche.overlap(z1, z2, cor=FALSE)
      data.frame(bin=b1, sp=sp, n=n, d=ovrlp$D, pa=z1$pa, pe=z1$pe, tol=z1$t)
    }
  }
}

# Calculate niche metrics -------------------------------------------------

# the older bin is column 1, the younger is column 2
bPairs <- cbind(bins[-1], bins[-length(bins)])
nichL <- lapply(spp, function(s){
  l <- apply(bPairs, 1, function(x){
    nicher(b1=x[1], b2=x[2], sp=s, env=env, h.method=h.method)
  })
  do.call(rbind, l)
})
nich <- do.call(rbind, nichL)
# remove NA rows
nas <- is.na(nich$bin)
nich <- nich[!nas,]

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')
dfNm <- paste0('Data/foram_niche_sumry_metrics_',day,'.csv')
# write.csv(nich, dfNm, row.names=FALSE)

# Evo models --------------------------------------------------------------

# nich <- read.csv("Data/foram_niche_sumry_metrics_200106.csv", stringsAsFactors=FALSE)

#nich$bin <- - nich$bin
ordr <- order(nich$bin)
nich <- nich[ordr,]

sp <- 'sampled' #for (sp in spp){}

spBool <- nich$sp==sp
sp <- nich[spBool,]
# ages must start at 0
sp$scaledT <- 1:nrow(sp) -1
ts <- as.paleoTS(mm = sp$pe, vv = sp$tol^2, nn = sp$n, tt = sp$scaledT, 
                 oldest = 'first', reset.time = FALSE)
mods <- fit9models(ts, method='Joint', silent=TRUE)

mods$modelFits
# From the package documentation:
# 'Method = "Joint" is a full likelihood approach, considering each time-series as a joint sample
# from a multivariate normal distribution. Method = "AD" is a REML approach that uses the differences
# between successive samples. They perform similarly, but the Joint approach does better under
# some circumstances (Hunt, 2008).'
# From Hunt 2008:
# 'joint parameter-ization  is  better  able  to  correctly  identify  directional  trends. 
# This advantage increases with sequence length, and is most pronounced when sampling error is high'

