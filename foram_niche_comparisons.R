#library(tidyr)
library(lme4)
library(adehabitatMA)
library(adehabitatHR)
library(ecospat)
library(paleoTS)
#library(VoCC)
#devtools::install_github("JorGarMol/VoCC", dependencies = FALSE, build_vignettes = FALSE)


# Evo models --------------------------------------------------------------

nich <- read.csv("Data/foram_niche_sumry_metrics_200106.csv", stringsAsFactors=FALSE)

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

