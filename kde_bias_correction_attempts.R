library(pracma)
library(GoFKernel)
library(doParallel)
library(tidyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

# Source functions for:
# length-bias corrected kernels
# Hellinger's H
# Holland's niche summarisation
# revised version of GoFKernel boundary reflection function
source('GSA_custom_ecospat_fcns.R')

# Simulate for simple bias ------------------------------------------------

# in Barmi & Simonoff 2000: n = 50 or 200
# for w=x, true Chi-square distribution is k=2 or 12
# for w=1/x, k=3 or 16
nVals <- c(50, 200)
kVals <- c(2, 12)
kVals2 <- c(3, 16)
#hMeth <- c('nrd0','SJ-ste')
#nbreak <- 2^8 # number of kernel estimation points
# using 2^9 gives oly 1x10^-6 better mise for one example

# original study used 500 replicates but calculation is slow
nreps <- 250
pkgs <- c('pracma','GoFKernel')
ncores <- detectCores() - 1

# recreate table 1 from Barmi & Simonoff 2000

pt1 <- proc.time()
registerDoParallel(ncores)
tab1a <- foreach(n=nVals, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:%
  foreach(k=kVals, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar% {
    miser <- function(n,k){
      w <- function(x){ x }
      x <- rchisq(n=n, df=(k+2))
      je <- JonesEst(x,w)
      te <- transformEst(x,w)
      mseJ <- function(x) (je$f(x) - dchisq(x, df=k))^2
      miseJ <- integral(mseJ, je$lower, je$upper)
      mseT <- function(x) (te$f(x) - dchisq(x, df=k))^2
      miseT <- integral(mseT, te$lower, te$upper)
      c(miseJ, miseT)
    }
    mises <- replicate(n=nreps, miser(n=n, k=k))
    mise <- rowMeans(mises)
    data.frame(n, k, t(mise))
  }
# TODO:
# Add tweak to uniroot function within inverse, for badly behaved empirical data:
# SJ-ste, n=200, k=16, w=1/x
# (or skip over error: 'f() values at end points not of opposite sign').
#tab1b <- foreach(n=nVals, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:%
  #  foreach(hMeth=hMeths, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:%
#  foreach(k=kVals2, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar% {
#    miser <- function(n,k){
#      w <- function(x){ 1/x }
#      x <- rchisq(n=n, df=(k-2))
#      je <- JonesEst(x,w)
#      te <- transformEst(x,w)
#      mseJ <- function(x) (je$f(x) - dchisq(x, df=k))^2
#      miseJ <- integral(mseJ, je$lower, je$upper)
#      mseT <- function(x) (te$f(x) - dchisq(x, df=k))^2
#      miseT <- integral(mseT, te$lower, te$upper)
#      c(miseJ, miseT)
#    }
#    mises <- replicate(n=nreps, miser(n=n, k=k))
#    mise <- rowMeans(mises)
#    data.frame(n, k, t(mise))
#  }
stopImplicitCluster()

pt2 <- proc.time()
pt2 - pt1 

colnames(tab1a) <- # colnames(tab1b) <- 
  c('n','k','MISEJones','MISEtrans')
tab1a$fmla <- 'w=x'
#tab1b$fmla <- 'w=1/x'
# write.csv(tab1a, 'Data/BarmiSimonoff_sim_250reps.csv', row.names=FALSE)

# Simulate polynomial bias ------------------------------------------------

# True distribution is chi-sq, df=5.
# Bias function is 4-th order polynomial,
# similar to foram sampling bias.
# Positive for values from 0 to 10.
w <- function(x){
  x <- .3*(x-5.5) 
  -0.8*x^4 + 0.7*x^3 + 1.5*x^2 - 2*x + 2
}
g <- function(x){
  g1 <- function(x) w(x)*dchisq(x, df=5)
  u <- integral(g1, 0, 10)
  w(x)*dchisq(x, df=5)/u
}

for (n in c(50, 200)){
  
set.seed(10)
x <- random.function(n=n, g, lower=0, upper=10)

je <- JonesEst(x,w)
mseJ <- function(x) (je$f(x) - dchisq(x, df=5))^2
miseJ <- integral(mseJ, je$lower, je$upper)

te <- transformEst(x,w)
mseT <- function(x) (te$f(x) - dchisq(x, df=5))^2
miseT <- integral(mseT, te$lower, te$upper)

# without accounting for bias:
kde <- density(x, kernel='gaussian', bw='nrd0', cut=0)
f <- approxfun(kde$x, kde$y)
mseBad <- function(x) (f(x) - dchisq(x, df=5))^2
miseBad <- integral(mseBad, min(kde$x), max(kde$x))

# miseJ; miseT; miseBad

xmn <- min(x)
xmx <- max(x)
u <- seq(xmn,xmx,by=.1)
kdeEsts <- data.frame(x=u)
kdeEsts$true <- dchisq(u, df=5)
kdeEsts$g <- g(u)
kdeEsts$je <- je$f(u)
kdeEsts$te <- te$f(u)
kdeEsts$biased <- f(u)
estsLong <- pivot_longer(kdeEsts, cols=c('true','g','biased','je','te'), 
                         names_to='method')
estsLong$method <- factor(estsLong$method, levels=c('true','g','biased','je','te'))
numNa <- nrow(estsLong) - length(x)
estsLong$rugHack <- c(x, rep(NA, numNa))
ymx <- max(estsLong$value, na.rm=TRUE) * 1.05
colr <- c('black','grey50','blue','orange','red')
names(colr) <- c('true','g','biased','je','te')

polyBiasP <- 
  ggplot(data=estsLong, aes(x=x)) +
  theme_bw() +
  scale_y_continuous(name='Density', limits=c(0,ymx), expand=c(0,0)) +
  scale_x_continuous(limits=c(xmn,xmx), expand=c(0,0)) +
  geom_line(aes(y=value, colour=method)) +
  geom_rug(aes(x=rugHack), sides='b') + 
  scale_colour_manual(values=colr)

txt <- c(paste('biased MISE =', round(miseBad,4)),
         paste('je MISE =', round(miseJ,4)),
         paste('te MISE =', round(miseT,4))
)
polyBiasP <- polyBiasP + 
  annotate('text', x=rep(xmx-1.5,3), y=ymx-c(0.025, 0.05, 0.075), label=txt)

polyNm <- paste0('Figs/KDE_example_chisq5_',n,'obs.pdf')
pdf(polyNm, width=6, height=4)
  print(polyBiasP)
dev.off()

}

# Boundary reflection -----------------------------------------------------

for (n in c(50, 200)){
  
# simulate a broad distribution that extends below 0
set.seed(20)
x <- rchisq(n, df=8) - 5
# hist(x, breaks=10)

lower <- 0
upper <- 10
outside <- which(x < lower | x > upper)
x <- x[-outside]

# Assume artificial boundaries at x = 0 and x = 10
# Do not specify cut; leave at default (3). This means
# the function will allow gradual tapering of density beyond 
# the reflected limits. The returned values will already
# be truncated to the given interval.
# Note that the given interval may be larger than the range of X,
# in which case the returned density will be over a larger range
# than if cut=0 were specified.
kdeRefl <- density.reflected(x, lower=lower, upper=upper)
fRefl <- approxfun(kdeRefl$x, kdeRefl$y)
kdeRaw <- density(x, from=lower, to=upper)
fRaw <- approxfun(kdeRaw$x, kdeRaw$y)

xmn <- min(x)
xmx <- max(x)
u <- seq(xmn,xmx,by=.1)
kdeEsts <- data.frame(x=u)
kdeEsts$true <- dchisq(u+5, df=8)
kdeEsts$refl <- fRefl(u)
kdeEsts$unrefl <- fRaw(u)

meths <- c('true','refl','unrefl')
estsLong <- pivot_longer(kdeEsts, cols=meths, 
                         names_to='method')
estsLong$method <- factor(estsLong$method, levels=meths)
numNa <- nrow(estsLong) - length(x)
estsLong$rugHack <- c(x, rep(NA, numNa))
ymx <- max(estsLong$value, na.rm=TRUE) * 1.05
colr <- c('black','orange','blue')
names(colr) <- meths

reflectP <- 
  ggplot(data=estsLong, aes(x=x)) +
  theme_bw() +
  scale_y_continuous(name='Density', limits=c(0,ymx), expand=c(0,0)) +
  scale_x_continuous(limits=c(xmn,xmx), expand=c(0,0)) +
  geom_line(aes(y=value, colour=method)) +
  geom_rug(aes(x=rugHack), sides='b') + 
  scale_colour_manual(values=colr)

#txt <- c(paste('biased MISE =', round(miseBad,4)),
#         paste('je MISE =', round(miseJ,4)),
#         paste('te MISE =', round(miseT,4))
#)
#reflectP <- reflectP + 
#  annotate('text', x=rep(xmx-1.5,3), y=ymx-c(0.025, 0.05, 0.075), label=txt)

reflNm <- paste0('Figs/KDE_example_chisq8shifted_trunc_',length(x),'obs.pdf')
pdf(reflNm, width=6, height=4)
print(reflectP)
dev.off()

}

# Empirical data ----------------------------------------------------------

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

occ <- read.csv('Data/foram_MAT_occs_latlong_8ka_trunc_200213.csv', stringsAsFactors=FALSE)
focalB <- seq(100, 700, by=200)
env <- 'temp_ym_0m'
spp <- unique(occ$species)
xmax <- max(occ[,env])
xmin <- min(occ[,env])

for (b in focalB){
  for (meth in c('weight','transform')){
    
    # estimate bias function (w) from sampling distribution, with boundary reflection
    sampRows <- which(occ$bin==b & occ$species=='sampled')
    samp <- occ[sampRows,env]
    densW <- density.reflected(samp, # kernel='gaussian', bw='nrd0', 
                               from=xmin, to=xmax, n = 2^9
    ) 
    w <- approxfun(densW$x, densW$y)
    
    plotL <- list()
    for (s in spp){
      # construct KDE for the sp, with boundary reflection and Jones' correction
      spRows <- which(occ$bin==b & occ$species==s)
      sp <- occ[spRows,env]
      
      if (meth=='weight'){
        kdeSp <- tryCatch(
          JonesEst(sp, w, reflect = TRUE, a = xmin, b = xmax),
          error = function(err){ list() }
        ) 
      }
      if (meth=='transform'){
        kdeSp <- tryCatch(
          transformEst(sp, w, reflect = TRUE, a = xmin, b = xmax),
          error = function(err){ list() }
        ) 
      }
      
      if (length(kdeSp)==0){ next }
      pePos <- which.max(kdeSp$y)
      kdeSp$pe <- kdeSp$x[pePos]
      
      # save a plot of the KDE
      plotDat <- data.frame(x=kdeSp$x, kd=kdeSp$y)
      nlab <- paste0('n=', length(sp))
      numNa <- nrow(plotDat) - length(sp)
      plotDat$rugHack <- c(sp, rep(NA, numNa))
      sPlot <- 
        ggplot(data=plotDat, aes(x=x, y=kd)) +
        theme_bw() +
        ggtitle(s) +
        geom_area(fill='cornsilk') +
        geom_line() +
        geom_segment(x=kdeSp$pe, xend=kdeSp$pe, y=0, yend=1, colour='red') +
        geom_segment(x=mean(sp), xend=mean(sp), y=0, yend=1, colour='blue') +
        geom_hline(yintercept=0) +
        geom_rug(aes(x=rugHack), sides='b', length=unit(0.045, "npc")) + 
        scale_x_continuous(expand=c(0,0)) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(size=9))
        sPlot <- sPlot +
          geom_text(label=nlab, size=3, # fontface=1,
                    x=8, y=0.9*max(plotDat$kd))
        plotL <- append(plotL, list(sPlot))
    } # loop through species
    
    # print page of plots
    pg <- plot_grid(plotlist=plotL, nrow=6)
    kdeNm <- paste0('Figs/KDE_all_spp_',meth,'_corrected_',b,'ka_',day,'.pdf')
    pdf(kdeNm, width=8.5, height=11)
      grid.arrange(arrangeGrob(pg))
    dev.off()
  } # loop through kernel correction methods
} # loop through bins
