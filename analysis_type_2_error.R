library(sp)
library(raster)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
library(cowplot)

# Data prep ---------------------------------------------------------------

source('species_kde_buildr.R')
day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")
nSim <- 10
bw <- 'SJ-ste'
envNm <- 'temp_ym'
nCore <- detectCores() - 1

# standardized sampling universe (MAT at core sites) at each of 4 depths
dList <- readRDS('Data/spp-and-sampling-data_list-by-depth_2020-07-21.rds')
df <- dList$temp_ym_0m$sp
bins <- unique(df$bin)
spp <- unique(df$sp)

# Define study intervals --------------------------------------------------

# glacial/interglacial ages from table S3
# excluding the oldest/youngest, for sample size problems
minAge <- c(28, 156, 268, 356, 436, 532) # 668
maxAge <- c(124, 212, 332, 412, 492, 588) # 4
ints <- data.frame(minAge = minAge, maxAge = maxAge)
xtremes <- c(minAge, maxAge)

# list every pairwise comparison for which to compute H distance during KDE
# (every warm vs. cold, warm vs. warm, and cold vs. cold interval)
wc <- expand.grid(X1 = ints$minAge, 
                  X2 = ints$maxAge)
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

# Define 12ka reference data ----------------------------------------------

# get occurrences at 12 ka
ka12bool <- df$bin == 12
ka12 <- df[ka12bool, ]

# extract cells for each species as a list
getCells <- function(s){
  sBool <- ka12$species == s
  ka12$cell[sBool]
}
cellList <- sapply(spp, getCells)

# extract temperature at all 12 ka cells for all glacial/interglacial bins

source('raster_brick_import_fcn.R')
modId <- read.csv('Data/gcm_model_codes.csv', stringsAsFactors=FALSE)
cells12ka <- unique(ka12$cell)

# Modified function from "foram occ data prep" script 
# to accept a list of raster cells (those from 12 ka)
# instead of using the actual/observed cells in each bin.
# Surface temperatures only.
addEnv <- function(bin, dat, mods, binCol, cells, env, dpths){ 
  slcEnv <- getBrik(bin = bin, envNm = env, mods = mods)
  envVals <- raster::extract(slcEnv[[env]], cells)
  # Rows = points of extraction, columns = depth layers  
  envVals <- envVals[,dpths]
  
  # Infer environment if it's missing and some of the adjacent 9 cells have values
  naVals <- sapply(envVals, function(x) any(is.na(x)) )
  if (sum(naVals) > 0){
    naCoords <- xyFromCell(slcEnv[[env]], cells[naVals])
    # distance to corner cell is 196 km for 1.25-degree resolution (~111 km/degree)
    fillr <- raster::extract(slcEnv[[env]], naCoords, 
                             buffer = 200*1000, fun = mean)
    envVals[naVals] <- fillr[,dpths]
  }
  envVals
}

pkgs <- c('sp','raster') 
registerDoParallel(nCore)
cellEnv <- foreach(bin=xtremes, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar%
  addEnv(bin = bin, dat = df, mods = modId, binCol = 'bin', cells = cells12ka, 
         env = envNm, dpths = 1) 
stopImplicitCluster()
# output is temperature as a matrix of bins (rows) by cells (columns)
row.names(cellEnv) <- xtremes
colnames(cellEnv) <- paste0('c', cells12ka)

naCount <- sum(is.na(cellEnv))
if (naCount > 0) stop('Deal with NA paleo-temperatures at 12 ka sites')

# Simulate past occurrences -----------------------------------------------

# pick the observed number of occs for a focal species in a focal bin,
# sampled from its 12 ka occurrences without replacement
spSimulatr <- function(s, dat, envVect, samplePot){
  slcN <- sum(dat[,'species'] == s)
  sCells <- samplePot[[s]]
  if (slcN > length(sCells)){
    stop(paste(s, 'too rare'))
  }
  simOccs <- sample(sCells, slcN, replace = FALSE)
  valNms <- paste0('c', simOccs)
  simEnv <- envVect[valNms]
  data.frame('species' = s, 'temp_ym' = simEnv, 'cell' = simOccs)
}

# replicate over species in a single bin
binSimulatr <- function(b, dat, samplePot, envMat){
  bBool <- dat[,'bin'] == b
  slc <- dat[bBool,]
  bSpp <- unique(slc[, 'species'])
  bEnv <- envMat[paste(b),]
  outL <- lapply(bSpp, spSimulatr, dat = slc, 
                 envVect = bEnv, samplePot = samplePot)
  outDf <- do.call('rbind', outL)
  outDf$bin <- b
  outDf[,c('species','bin','temp_ym')]
}

# replicate over bins in a single simulation iteration
simulatr <- function(bins, dat, samplePot, envMat){
  binL <- lapply(xtremes, binSimulatr, dat = df, 
                 samplePot = samplePot, envMat = envMat)
  df <- do.call('rbind', binL)
  list(df)
}

simList <- replicate(nSim,
  simulatr(bins = xtremes, dat = df, samplePot = cellList, envMat = cellEnv)
)

# KDE ---------------------------------------------------------------------

# determine the standard axis limits for the sea surface
samp <- dList[[1]]$samp
sampSmry <- sapply(bins, minmax, df = samp, env = envNm)
xmx <- min(sampSmry[2,])
xmn <- max(sampSmry[1,])

# calculate sampling density in the single reference bin (12 ka)
sampRows12 <- which(samp$bin == 12)
samp12 <- samp[sampRows12, envNm]
densSamp12 <- density(samp12, bw = bw)
w12 <- approxfun(densSamp12$x, densSamp12$y)

# Modify the kde function from "species kde buildr" script so that it 
# takes the given sampling bias function (constructed for the 12 ka bin).
# That way there's no need to supply 'samp' in the data list object,
# or spend time estimatating sampling density in each bin.
kde4sim <- function(dat, bPair, envNm, bw = 'nrd0', xmn, xmx, w){
  b1 <- bPair[1]
  b2 <- bPair[2]
  zoneSp <- unique(dat$species)
  sList <- lapply(zoneSp, function(s){
    nicher(dat = list(sp = dat), bw = bw, s = s, envNm = envNm,
           b1 = b1, b2 = b2, w1 = w, w2 = w, xmn = xmn, xmx = xmx)
  })
  do.call(rbind, sList)
}

pairL <- list()
for (i in 1:nrow(intPairs)){
  b1 <- intPairs$t1[i]
  b2 <- intPairs$t2[i]
  entry <- c(b1, b2)
  pairL <- append(pairL, list(entry))
}

# warning - this could take hours
pkg <- c('pracma','GoFKernel','kerneval')
pt1 <- proc.time()
registerDoParallel(nCore)
kdeSum <- foreach(dat=simList, .inorder=FALSE, .packages=pkg) %:%
  foreach(bPair=pairL, .combine=rbind, .inorder=FALSE, .packages=pkg) %dopar%
    kde4sim(dat, bPair, envNm, bw = bw, xmn = xmn, xmx = xmx, w = w12)
stopImplicitCluster()
pt2 <- proc.time()
(pt2-pt1)/60

# each simulation run will be a different dataframe in the list
# note that there will be all-NA rows in each dataframe to remove later
outNm <- paste0('Data/niche-xtremes-simulations_',bw,'_SS_',day,'.rds')
saveRDS(kdeSum, outNm)

# TIME-SAVNG OPTION:
# if the script has already been run once, the niche summaries were exported
# so read them in here instead of running the code chunk above

# kdeSum <- readRDS('Data/niche-xtremes-simulations_SJ-ste_SS_2020-08-04.rds')
# kdeSum <- readRDS('Data/niche-xtremes-simulations_SJ-ste_SS_2020-08-05.rds')

# Results -----------------------------------------------------------------

# take the mean H for each bin combination (rows)
# repeated for each iteration (columns)
sumH <- function(bPair, dat){
  intBool <- which(dat$bin == bPair[1] & dat$bin2 == bPair[2])
  int <- dat[intBool,]
  mean(int$h, na.rm = TRUE)
}
Hlist <- lapply(kdeSum, function(dat){
  tempList <- lapply(pairL, sumH, dat = dat)
  Hvect <- unlist(tempList)
})
Hmat <- do.call(cbind, Hlist)

# count the proportion of trials where signal of niche change is found
ccRows <- which(intPairs$type == 'cold-cold')
wwRows <- which(intPairs$type == 'warm-warm')
cwRows <- which(intPairs$type == 'cold-warm')
isDetected <- apply(Hmat, 2, function(dat){
  ccAvg <- mean(dat[ccRows])
  wwAvg <- mean(dat[wwRows])
  cwAvg <- mean(dat[cwRows])
  cwAvg > ccAvg & cwAvg > wwAvg
})
sum(isDetected)/nSim

# multi-panel plot

colr <- c('cold-cold' = 'deepskyblue',
          'warm-warm' = 'firebrick2',
          'cold-warm' = 'purple3')
plotdat <- cbind(intPairs, Hmat)
colnames(plotdat)[-(1:3)] <- paste0('rep',1:nSim)

boxer <- function(iter, noX = FALSE, noY = FALSE){
  colNm <- paste0('rep', iter)
  dat <- plotdat[,c('type', colNm)]
  colnames(dat) <- c('type', 'y')
  p <- ggplot(data = dat) +
    theme_bw() +
    scale_y_continuous('Mean niche H distance',
      limits = c(0, 0.5), expand = c(0, 0)) +
    geom_boxplot(aes(x = type, y = y, fill = type)) +
    scale_fill_manual(values = colr) +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 6)
          )
  if (noX){
    p <- p +
      theme(axis.text.x = element_blank())
  }
  if (noY){
    p <- p +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank())
  }
  p
}

# fine-tune the plotting parameters
# don't duplicate axis labels in interior plots
topL <- (boxer(1, noX = TRUE))
topC <- boxer(2, noX = TRUE, noY = TRUE)
topR <- boxer(3, noX = TRUE, noY = TRUE)
lwrL <- boxer(4)
lwrC <- boxer(5, noY = TRUE)
lwrR <- boxer(6, noY = TRUE)
boxNm <- paste0('Figs/H-vs-climate-extreme_simulated_',day,'.pdf')
pdf(boxNm, width = 6, height =4)
plot_grid(
  topL, topC, topR,
  lwrL, lwrC, lwrR,
  ncol = 3, rel_widths = c(1.2, 1, 1),
  labels = 'AUTO', label_size = 12,
  hjust = c(-5, -2, -2, -5, -2, -2.4), 
  vjust = c(2, 2, 2, 1.85, 1.85, 1.85)
)
dev.off()

