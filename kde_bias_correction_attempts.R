library(pracma)
library(GoFKernel)
library(doParallel)

transformEst <- function(x, w, hMeth='nrd0', nbreak=2^8){
  n <- length(x)
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
  kde <- density(Y, kernel='gaussian', bw=hMeth, cut=0, n=nbreak)
  
  # scale to integrate to one
  kde$y <- muHat * kde$y
  
  # Back-transform from the y=W(x) argument to x.
  # Build custom cdf function, since 
  # inverse() does not play well with ecdf (step fcn)
  upr <- max(kde$x)
  cdf <- function(f, lower) {
    function(z) integral(f, lower, z)
  }
  wcdf <- cdf(f = w, lower = lwr)
  inv <- inverse(wcdf, lower = lwr, upper = upr)
  
  kde$xTrans <- kde$x
  kde$x <- sapply(kde$x, inv)
  
  # if w(x) is undefined at 0, then min and max of X can off
  # (albeit only slightly) from values at which kernel can estimate 
  f <- approxfun(kde$x, kde$y)
  a <- min(kde$x)
  b <- max(kde$x)
  
  append(list(f=f, lower=a, upper=b), kde)
}

# weighted kernel density estimation after Jones 1991
JonesEst <- function(x, w, hMeth='nrd0', nbreak=2^8){
  wts <- 1/w(x)
  wts <- wts/sum(wts)
  kde <- density(x, kernel='gaussian', bw=hMeth, cut=0, n=nbreak,
                 weights=wts)
  f <- approxfun(kde$x, kde$y)
  a <- min(kde$x)
  b <- max(kde$x)
  append(list(f=f, lower=a, upper=b), kde)
}

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

# TODO
# Use a true distribution other than chi-sq
# Boundary reflection
# Use Borrajo's rule of thumb and 2 bootstrap bandwidth estimators



# Simulate polynomial bias ------------------------------------------------

#True distribution is chi-sq, df=5.
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
n <- 50
set.seed(10)
x <- random.function(n=n, g, lower=0, upper=10)

je <- JonesEst(x,w)
mseJ <- function(x) (je$f(x) - dchisq(x, df=5))^2
miseJ <- integral(mseJ, je$lower, je$upper)

te <- transformEst(x,w)
mseT <- function(x) (te$f(x) - dchisq(x, df=5))^2
miseT <- integral(mseT, te$lower, te$upper)

# without accounting for bias:
kde <- density(x, kernel='gaussian', bw='nrd0', cut=0, n=2^8)
f <- approxfun(kde$x, kde$y)
mseBad <- function(x) (f(x) - dchisq(x, df=5))^2
miseBad <- integral(mseBad, min(kde$x), max(kde$x))

miseJ; miseT; miseBad

xtrue <- rchisq(50, df=5)
toss <- xtrue> 10
xtrue <- xtrue[!toss]
hist(xtrue, xlim=c(0,10), breaks=seq(0,10,by=1), freq=FALSE)
hist(x, xlim=c(0,10), freq=FALSE)
u <- seq(0,10,by=.1)
plot(u,w(u))

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
