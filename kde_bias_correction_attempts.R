library(pracma)
library(GoFKernel)
library(doParallel)

# Simulate data -----------------------------------------------------------

# calculate MISE for chi-sq distribution with given n and df
# based on data drawn with known bias, w
miser <- function(x, n, k, hMeth, nbreak, w){
  
  muHat <- n * sum( w(x)^-1 )^-1
  
  test <- c(is.infinite(w(0)), is.na(w(0)))
  if (any(test)){
    lwr <- exp(-30)
  } else {
    lwr <- 0
  }
  Y <- sapply(x, function(u){
    integral(w, xmin = lwr, xmax = u) 
  })
  
  # kernel density estimation on transformed data
  # note that the authors used local quadratic likelihood instead
  kdeTrans <- density(Y, kernel='gaussian', bw=hMeth, cut=0, n=nbreak)
  
  # scale to integrate to one
  kdeTrans$y <- muHat * kdeTrans$y
  
  # Back-transform from the y=W(x) argument to x.
  # Build custom cdf function, since 
  # inverse() does not play well with ecdf (step fcn)
  upr <- max(kdeTrans$x)
  cdf <- function(f, lower) {
    function(z) integral(f, lower, z)
  }
  wcdf <- cdf(f = w, lower = lwr)
  inv <- inverse(wcdf, lower = lwr, upper = upr)
  
  kdeTrans$xTrans <- kdeTrans$x
  # this is the bottleneck step:
  kdeTrans$x <- sapply(kdeTrans$x, inv)
  # plot(kdeTrans)
  
  # calculate MISE (minimum integrated squared error)
  
  # if w(x) is undefined at 0, then min and max of X can off
  # (albeit only slightly) from values at which kernel can estimate 
  fHatT <- approxfun(kdeTrans$x, kdeTrans$y)
  seTrans <- function(u) (fHatT(u) - dchisq(u, df=k))^2
  a <- min(kdeTrans$x)
  b <- max(kdeTrans$x)
  miseTrans <- integral(seTrans, a, b)
  
  # weighted kernel density estimation after Jones 1991
  wts <- 1/w(x)
  wts <- wts/sum(wts)
  kdeJ <- density(x, kernel='gaussian', bw=hMeth, cut=0, n=nbreak,
                  weights=wts)
  fhatJ <- approxfun(kdeJ$x, kdeJ$y)
  seJ <- function(x) (fhatJ(x) - dchisq(x, df=k))^2
  miseJ <- integral(seJ, a, b)
  
  c(miseJ, miseTrans)
}

# in Barmi & Simonoff 2000: n = 50 or 200
# for w=x, true Chi-square distribution is k=2 or 12
# for w=1/x, k=3 or 16
nVals <- c(50, 200)
kVals <- c(2, 12)
kVals2 <- c(3, 16)
hMeths <- c('nrd0','SJ-ste')
nbreak <- 2^8 # number of kernel estimation points
# using 2^9 gives oly 1x10^-6 better mise for one example

# original study used 500 replicates but calculation is slow
nreps <- 100
pkgs <- c('pracma','GoFKernel')
ncores <- detectCores() - 1

# recreate table 1 from Barmi & Simonoff 2000

pt1 <- proc.time()
registerDoParallel(ncores)
tab1a <- foreach(n=nVals, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:%
  foreach(hMeth=hMeths, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:%
  foreach(k=kVals, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar% {
    w <- function(x){ x }
    x <- rchisq(n=n, df=(k+2))
    mises <- replicate(nreps, miser(x,n,k,hMeth,nbreak,w))
    mise <- colMeans(t(mises))
    data.frame(hMeth, n, k, t(mise))
  }
tab1b <- foreach(n=nVals, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:%
  foreach(hMeth=hMeths, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:%
  foreach(k=kVals2, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar% {
    w <- function(x){ 1/x }
    x <- rchisq(n=n, df=(k-2))
    mises <- replicate(nreps, miser(x,n,k,hMeth,nbreak,w))
    mise <- colMeans(t(mises))
    data.frame(hMeth, n, k, t(mise))
  }
stopImplicitCluster()

pt2 <- proc.time()
pt2 - pt1 

colnames(tab1a) <- colnames(tab1b) <- 
  c('hMeth','n','k','MISEJones','MISEtrans')
tab1a$fmla <- 'w=x'
tab1b$fmla <- 'w=1/x'
tab1 <- rbind(tab1a,tab1b)
# write.csv(tab1b, 'Data/BarmiSimonoff_sim_replications.csv', row.names=FALSE)

# TODO
# use a true distribution other than chi-sq
# use a non-monotonic bias function
# introduce boundary bias

# Empirical data ----------------------------------------------------------

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

occ <- read.csv('Data/foram_MAT_occs_latlong_8ka_200129.csv', stringsAsFactors=FALSE)
bins <- unique(occ$bin) 

source('GSA_custom_ecospat_fcns.R')
nbins <- length(bins)
h.method <- "nrd0" # "SJ-ste" # "ucv"
R <- 2^8

env <- 'temp_ym_0m'
spp <- unique(occ$species)
xmax <- max(occ[,env])
xmin <- min(occ[,env])

b <- 100; s <-'sampled'#  spp[2]

spRows <- which(occ$bin==b & occ$species==s)
sp <- 
  samp <- occ[spRows,env]

# Estimate kernel density -------------------------------------------------

#kdeNiche <- 
#  function (sp, samp=NULL, xmax, xmin, R, h.method="nrd0", weight=FALSE) {
sp <- as.matrix(sp)
x <- seq(from = xmin, to = xmax, length.out = R)

if (weight==TRUE){
  samp <- as.matrix(samp)
  kdeSamp <- density(samp, kernel='gaussian', bw=h.method, 
                     from=xmin, to=xmax, n=R)
  # integrate in segments of x; use these as weights for species KDE
  sampFun <- approxfun(x, kdeSamp$y)
  bounds <- data.frame(a=x[-length(x)], b=x[-1])
  pmf <- apply(bounds, 1, function(ab){
    integral(sampFun, ab[1], ab[2])
  })
  bounds$pmf <- pmf
  
  # find weight for locations of observations
  wt <- sapply(sp, function(u){
    r <- which(bounds$a <= u & bounds$b > u)
    1/bounds$pmf[r]
  } )
  wt <- wt/sum(wt)
  kdeSp <- density(sp, kernel='gaussian', bw=h.method, 
                   from=xmin, to=xmax, n=R,
                   weights=wt)
} else {
  kdeSp <- density(sp, kernel='gaussian', bw=h.method, 
                   from=xmin, to=xmax, n=R)
}

plot(kdeSamp)
plot(kdeSp)

# format function output
l <- list()
l$density <- kdeSp$y

# Add Holland 1995 niche summary metrics to output:
# peak abundance, preferred environment, & tolerance
l$pa <- max(kdeSp$y)
pe.pos <- which.max(kdeSp$y)
l$pe <- x[pe.pos]

# To approximate the KDE as a function, one also needs x-values:
l$x <- x

#    return(l)
#  }
