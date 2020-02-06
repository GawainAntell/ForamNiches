library(pracma)
library(GoFKernel)

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
  wcdf <- function(u){
    integral(w, xmin = lwr, xmax = u) 
  }
  wInv <- inverse(wcdf, lower=lwr)
  
  kdeTrans$xTrans <- kdeTrans$x
  # this is the bottleneck step:
  kdeTrans$x <- sapply(kdeTrans$x, wInv)
  # plot(kdeTrans)
  
  # calculate MISE (minimum integrated squared error)
  
  # if w(x) is undefined at 0, then min(X) can be larger
  # (albeit only slightly) than values at which kernel can estimate 
  a <- min(kdeTrans$x) 
  b <- max(x)
  fHatT <- approxfun(kdeTrans$x, kdeTrans$y)
  seTrans <- function(u) (fHatT(u) - dchisq(u, df=k))^2
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
kVals <- c(3, 16)
  # c(2, 12)
hMeths <- c('nrd0','SJ-ste')
w <- function(x){ 
  1/x 
#   x # classical length bias
  }
nbreak <- 2^8 # number of kernel estimation points
# using 2^9 gives oly 1x10^-6 better mise for one example

# for TESTING ONLY
# n <- 50; k <- 3; hMeth <- 'nrd0'

# recreate table 1 from Barmi & Simonoff 2000
tab1 <- data.frame(matrix(ncol=5, nrow=0))
for (n in nVals){
  for (k in kVals){
    for (hMeth in hMeths){
      # when w(x)=x, then g is a chi-sq density of k+2
      # when w(x)=1/x, then chi-sq has k-2
      x <- rchisq(n=n, df=(k-2))
      # rchisq(n=n, df=(k-2))
      
      # original study used 500 replicates but calculation is slow
      mises <- replicate(100, miser(x,n,k,hMeth,nbreak,w))
      mise <- colMeans(t(mises))
      newRow <- data.frame(hMeth, n, k, t(mise))
      tab1 <- rbind(tab1, newRow)
    }
  }
}
colnames(tab1) <- c('hMeth','n','k','MISEJones','MISEtrans')



