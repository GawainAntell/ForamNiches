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

n <- 200
  #  upr <- 10
    
# define true species density (f) and sampling bias (w) functions

# true Chi-square distribution is k=2 in Barmi & Simonoff 2000
  #f <- rchisq(n=n, df=2)
  #hist(f)

  # Alternatively, more complicated gamma function
  # f <- function(x, a, s){
  #   1/(s^a * gamma(a)) * x^(a-1) * exp(1)^-(x/s)
  # } 
  # fx <- f(seq(0,1,by=0.01), 2, 0.25)
  # plot(seq(0,1,by=0.01), fx)

# classical length-bias, Cox 1969
w <- function(x){ x }

  #gUnscld <- function(x){
  #  dchisq(x, df=2) * w(x)
    # f(x,a=2,s=0.25) * w(x) # gamma distribution
  #}
  #u <- integrate(gUnscld, lower=0, upper=upr)$value
  #g <- function(x){ 
  #  gUnscld(x) / u
  #}
# simulate data under g, where g=f*w/u
  #gIn <- runif(n, min=0, max=upr)
  #x <- g(gIn)

# when w(x)=x, then g is a chi-sq density of k+2
X <- rchisq(n=n, df=4)

# Barmi and Simonoff transform --------------------------------------------

uhat <- n * sum( w(X)^-1 )^-1
Y <- sapply(X, function(x){
  integrate(w, lower = 0, upper = x)$value
})

# kernel density estimation on transformed data
kde <- density(Y, kernel='gaussian', bw='nrd0', cut=0) #    from=xmin, to=xmax, n=R

# scale to integrate to one
kde$y <- uhat * kde$y

# back-transform from the y=W(x) argument to x
kde$xTrans <- kde$x
kde$x <- sqrt(2 * kde$x)
plot(kde)

# calculate minimum integrated squared error
# when the smoothing parameter is chosen
# in each replication to minimize ISE.

