#library(tidyr)
library(lme4)
library(adehabitatMA)
library(adehabitatHR)
library(ecospat)
#library(VoCC)
#devtools::install_github("JorGarMol/VoCC", dependencies = FALSE, build_vignettes = FALSE)
source('ecospat.grid.clim.dyn.GSA.fcn.R')

df <- read.csv('Data/foram_MAT_occs_latlong_8ka_200103.csv',stringsAsFactors=FALSE)
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
    data.frame(bin=NA, sp=NA, d=NA, pa=NA, pe=NA, tol=NA)
  } else {
    if (length(z2)==0){
      data.frame(bin=b1, sp=sp, d=NA, pa=z1$pa, pe=z1$pe, tol=z1$t)
    } else{
      ovrlp <- GSA.ecospat.niche.overlap(z1, z2, cor=FALSE)
      data.frame(bin=b1, sp=sp, d=ovrlp$D, pa=z1$pa, pe=z1$pe, tol=z1$t)
    }
  }
}

# the older bin is column 1, the younger is column 2
bPairs <- cbind(bins[-1], bins[-length(bins)])
nichL <- lapply(spp, function(s){
  l <- apply(bPairs, 1, function(x){
    nicher(b1=x[1], b2=x[2], sp=s, env=env, h.method=h.method)
  })
  do.call(rbind, l)
})
nichDf <- do.call(rbind, nichL)

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')
dfNm <- paste0('Data/foram_niche_sumry_metrics_',day,'.csv')
write.csv(nichDf, dfNm, row.names=FALSE)

####################################################################
# Data visualisation: least-squares cross-validation vs. subjective h

# Use 3 time intervals as examples
focalPos <- round(nbins*c(0.25,0.5,0.75))
focalPts <- bins[focalPos]
binL <- bins[2]-bins[1]

for (pt in focalPts){
  b1 <- pt + binL
  b2 <- pt
  slcBool <- df$bin %in% c(b1, b2)
  slc <- df[slcBool,]
  glob <- slc[, c('pc1','pc2')]
  glob1bool <- slc$bin==b1
  glob1 <- glob[glob1bool,]
  glob2bool <- slc$bin==b2
  glob2 <- glob[glob2bool,]
  
  # Compare the rarest vs. commonest species
  spList1 <- unique(slc$species[glob1bool])
  spList2 <- unique(slc$species[glob2bool])
  spBoth <- intersect(spList1, spList2)
  spBoth <- setdiff(spBoth, 'sampled')
  slcBoth <- slc$species[slc$species %in% spBoth]
  nmTbl <- sort(table(slcBoth))
  rare <- names(which.min(nmTbl))
  abund <- names(which.max(nmTbl))

# Compare subjective vs. least-squares cross-validation bandwidth
for (sp in c(rare, abund)){
  sp1rows <- which(slc$species==sp & slc$bin==b1)
  sp1 <- glob[sp1rows,]
  sp2rows <- which(slc$species==sp & slc$bin==b2)
  sp2 <- glob[sp2rows,]
  
  glob <- as.matrix(glob)
  glob1 <- as.matrix(glob1)
  glob2 <- as.matrix(glob2)
  mask <- ascgen(SpatialPoints(cbind((0:(R))/R, (0:(R)/R))), 
                 nrcol = R - 2, count = FALSE)
  
  xmin <- min(glob[, 1])
  xmax <- max(glob[, 1])
  ymin <- min(glob[, 2])
  ymax <- max(glob[, 2])

  glob1r <- data.frame(cbind((glob1[, 1] - xmin)/abs(xmax - xmin), 
                             (glob1[, 2] - ymin)/abs(ymax - ymin))
  )
  glob2r <- data.frame(cbind((glob2[, 1] - xmin)/abs(xmax - xmin), 
                             (glob2[, 2] - ymin)/abs(ymax - ymin))
  )
  spr1 <- data.frame(cbind((sp1[, 1] - xmin)/abs(xmax - xmin), 
                           (sp1[, 2] - ymin)/abs(ymax - ymin)))
  spr2 <- data.frame(cbind((sp2[, 1] - xmin)/abs(xmax - xmin), 
                           (sp2[, 2] - ymin)/abs(ymax - ymin)))

  # Follow Silverman's recommendation of using
  # least-squares cross-validation within a range of the subjective h.
  # See getvolumeHD from kernelUD for argument details.
  sp1dens <- kernelUD(SpatialPoints(spr1), h = "LSCV", 
                      hlim=c(0.5, 2), grid = mask, kern = "bivnorm")
  sp2dens <- kernelUD(SpatialPoints(spr2), h = "LSCV", 
                      hlim=c(0.5, 2), grid = mask, kern = "bivnorm")
  #glob1dens <- kernelUD(SpatialPoints(glob1r), h = "LSCV", 
  #                       hlim=c(0.5, 2), grid = mask, kern = "bivnorm")
  # save h values to plot
  getH <- function(dens, spr){
    lscv <- dens@h$h
    subjSigma <-  0.5*(var(spr$X1)+var(spr$X2))
    subj <- sqrt(subjSigma) * nrow(spr)^(-1/6)
    out <- c(lscv, subj)
    names(out) <- c('lscv','subjective')
    out
  }
  sp1h <- getH(sp1dens, spr1)
  sp2h <- getH(sp2dens, spr2)
  
  # Subjective bandwidth
  z1 <- ecospat.grid.clim.dyn(glob, glob1, sp1, R, th.sp=0.05, th.env=0.05)
  z2 <- ecospat.grid.clim.dyn(glob, glob2, sp2, R, th.sp=0.05, th.env=0.05)
  ovrlp <- ecospat.niche.overlap(z1, z2, cor=TRUE)
  d <- ovrlp$D 
  maxY <- max(z1$y)*0.8
  
  # Least-squares cross-validation bandwidth
  z1lscv <- suppressWarnings(
    ecospat.grid.clim.dyn.GSA(glob, glob1, sp1, R, th.sp=0.05, th.env=0.05)
  )
  z2lscv <- suppressWarnings(
    ecospat.grid.clim.dyn.GSA(glob, glob2, sp2, R, th.sp=0.05, th.env=0.05)
  )
  ovrlpLscv <- ecospat.niche.overlap(z1lscv, z2lscv, cor=TRUE)
  dLscv <- ovrlpLscv$D 
  maxYlscv <- max(z1lscv$y)*0.8
  
  # plot PC1 vs PC2 for the subjective and lscv bandwidths
  # list overlap on the plot, for the next time step
  if (sp==rare){suffix <- 'rare'} else {suffix <- 'common'}
  sp <- gsub(' ', '.', sp)
  pNm <- paste0('Figs/niche_heatmap_', b1, 'ka_', suffix, '_', sp, '_', day, '.pdf')
  pdf(pNm, width=4, height=8) 
    par(mfrow=c(3,1), mar=c(3,3,3,3))
    plotLSCV(sp1dens)
    abline(v=sp1h['subjective'], col='red')
    abline(v=sp1h['lscv'], col='blue')
    legend('bottomleft', lty=1, 
           col=c('red','blue'),
           legend = c('Default','LSCV'))
    plot(z1$z, main='Default h', xlab='PC1', ylab='PC2')
    lbl <- paste('overlap =', round(d,3))
    text(0, maxY, lbl)
    plot(z1lscv$z, main='LSCV h', xlab='PC1')
    lblLscv <- paste('overlap =', round(dLscv,3))
    text(0, maxYlscv, lblLscv)
  dev.off()
  
} # end loop through species
} # end loop through focal points

