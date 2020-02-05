library(pracma)

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

# in Barmi & Simonoff 2000:
# true Chi-square distribution is k=2 or 12; n = 50 or 200
nVals <- c(50, 200)
kVals <- c(2, 12)
hMeths <- c('nrd0','SJ-ste')

# calculate MISE for chi-sq distribution with given n and df
simBS <- function(n, k, hMeth){
  w <- function(x){ x }
  
  # when w(x)=x, then g is a chi-sq density of k+2
  X <- rchisq(n=n, df=k+2)
  
  uhat <- n * sum( w(X)^-1 )^-1
  Y <- sapply(X, function(x){
    integrate(w, lower = 0, upper = x)$value
  })
  
  # kernel density estimation on transformed data
  # note that the authors used local quadratic likelihood instead
  kdeTrans <- density(Y, kernel='gaussian', bw=hMeth, cut=0)
  
  # scale to integrate to one
  kdeTrans$y <- uhat * kdeTrans$y
  
  # back-transform from the y=W(x) argument to x
  kdeTrans$xTrans <- kdeTrans$x
  kdeTrans$x <- sqrt(2 * kdeTrans$x)
  # plot(kde)
  
  # calculate MISE (minimum integrated squared error)
  a <- min(X)
  b <- max(X)
  fhat <- approxfun(kdeTrans$x, kdeTrans$y)
  se <- function(x) (fhat(x) - dchisq(x, df=k))^2
  mise <- integral(se, a, b)
  
  # kernel density estimation after Jones 1991
  
  # calculate  weights
  wts <- 1/w(X)
  wts <- wts/sum(wts)
  kdeJ <- density(X, kernel='gaussian', bw=hMeth, cut=0,
                  weights=wts)
  fhatJ <- approxfun(kdeJ$x, kdeJ$y)
  seJ <- function(x) (fhatJ(x) - dchisq(x, df=k))^2
  miseJ <- integral(seJ, a, b)
  
  c(miseJ, mise)
}

# recreate table 1 from Barmi & Simonoff 2000
tab1 <- data.frame(matrix(ncol=5, nrow=0))
for (n in nVals){
  for (k in kVals){
    for (hMeth in hMeths){
      mises <- replicate(500, simBS(n,k,hMeth))
      mise <- colMeans(t(mises))
      newRow <- data.frame(hMeth, n, k, t(mise))
      tab1 <- rbind(tab1, newRow)
    }
  }
}
colnames(tab1) <- c('hMeth','n','k','MISEJones','MISEtrans')



