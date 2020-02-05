library(pracma)

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

