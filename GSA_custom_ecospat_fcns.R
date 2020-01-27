
GSA.grid.clim.dyn <- 
  function (glob, glob1, sp, R, th.sp = 0, th.env = 0, geomask = NULL, h.method="nrd0") {
    glob <- as.matrix(glob)
    glob1 <- as.matrix(glob1)
    sp <- as.matrix(sp)
    l <- list()
    xmax <- max(glob[, 1])
    xmin <- min(glob[, 1])
    x <- seq(from = xmin, to = xmax, length.out = R)
    # Arguments 'from', 'to', 'R', and 'cut' don't affect bandwidth for nrd0, SJ-ste, or ucv
    # but they do affect the roughness of the density plot - smoother (better) with defaults
    # for SJ-ste and ucv especially. Customised density behaviour is ok for nrd0.
    # Argument 'cut' is manually set to 0 in the ecospat code. Default is 3. The cut specifies
    # how quickly the density drops to 0 outside of the data range; since the data are
    # truncated to a constant temperature range across time bins, it makes sense to allow
    # the density to extend a bit past the extreme values. This is an argument to NOT
    # set cut = 0. HOWEVER, in order to compare niche overlap, the density functions
    # must be integrable over the same interval. This is readily achieved with cut=0.
    # Note that 'from' and 'to' arguments are a redundant way to specify the cut (?).
    # Setting from=xmin, to=xmax is the same as cut=0.
    z <- density(sp[, 1], 
                       kernel = "gaussian", 
                       bw=h.method,
                       from = xmin, 
                       to = xmax, 
                       n = R
                  #     cut = 0
                       )
    glob1.dens <- density(glob1[, 1], 
                          kernel = "gaussian", 
                          bw=h.method,
                          from = xmin, 
                          to = xmax, 
                          n = R
                  #        cut = 0
                          )
    # Original ecospat code re-scales the distribution by the max density value.
    # The reason why is unclear, and this step does very slightly change Schoener's D.
 #   z <- z$y # * nrow(sp)/sum(z$y)
    Z <- glob1.dens$y # * nrow(glob)/sum(glob1.dens$y)
    
    # Original ecospat code truncated the KDE by the quantile supplied.
    # This functionality is not used here, and setting quantile = 0 leads to bad behaviour,
    # so this section is commented out.
      # glob1r <- sapply(glob1, findInterval, glob1.dens$x)
      # th.env <- quantile(glob1.dens$y[glob1r], th.env)
      # glob1rm <- which(Z < th.env)
      # spr <- sapply(sp, findInterval, z$x)
      # th.sp <- quantile(z$y[spr], th.sp)
      # sprm <- which(z < th.sp)
      # z[sprm] <- 0
      # Z[glob1rm] <- 0
    
    # z.uncor <- z/max(z)
    z.cor <- z$y/Z
    z.cor[is.na(z.cor)] <- 0
    z.cor[z.cor == "Inf"] <- 0
    # z.cor <- z.cor/max(z.cor)
    
    # These outputs are commented out because they aren't needed for later analysis.
#    w <- z.uncor
#    w[w > 0] <- 1
#    l$x <- x
    l$z <- z$y
#    l$z.uncor <- z.uncor
    l$z.cor <- z.cor
    l$Z <- Z
#    l$glob <- glob
#    l$glob1 <- glob1
#    l$sp <- sp
#    l$w <- w
    
    # Add Holland 1995 niche summary metrics to output:
    # peak abundance, preferred environment, & tolerance
    l$pa <- max(z$y)
    pe.pos <- which.max(z$y)
    l$pe <- x[pe.pos]
    l$t <- sd(z$x)
    
    # To approximate the KDE as a function, one also needs x-values:
    l$x <- z$x
    
    return(l)
  }

GSA.ecospat.niche.overlap <- function (z1, z2, cor) {
  if (cor == FALSE) {
    # The original code scales density estimates by the max value,
    # then takes sum as if it were a discrete mass function.
#    p1 <- as.matrix(z1$z.uncor)/sum(as.matrix(z1$z.uncor))
#    p2 <- as.matrix(z2$z.uncor)/sum(as.matrix(z2$z.uncor))
    # However, the KDE approximates a continuous density function.
    # So, find that function and integrate it.
    # integrating the density distribution as a continuous function:
    fun1 <- approxfun(z1$x, z1$z)
    fun2 <- approxfun(z2$x, z2$z)
    # Package warning:
    # 'The value returned by approxfun contains references to the code in the current version of R: 
    # it is not intended to be saved and loaded into a different R session. This is safer for R >= 3.0.0.'
    # Also note: since the functions integrate to < 1, the minimum achievable D value is > 0.
  }
  if (cor == TRUE) {
    fun1 <- approxfun(z1$x, z1$z.cor)
    fun2 <- approxfun(z2$x, z2$z.cor)
#    p1 <- as.matrix(z1$z.cor)/sum(as.matrix(z1$z.cor))
#    p2 <- as.matrix(z2$z.cor)/sum(as.matrix(z2$z.cor))
  }
  
  difFun <- function(x) abs(fun1(x) - fun2(x))
  xmin <- min(z1$x)
  xmax <- max(z1$x)
  int <- integrate(difFun, xmin, xmax)
  D <- 1 - int$value * 0.5
  return(D)
}

# TODO find out why all of the below functions have 'extremely bad integrand behaviour'
# note that subdivisions should be > R

d1 <- function (z1, z2) {
  fun1 <- approxfun(z1$x, z1$z)
  fun2 <- approxfun(z2$x, z2$z)
  difFun <- function(x) abs(fun1(x) - fun2(x))
  a <- z1$x[2] #max(min(z1$x), min(z2$x))
  b <- z1$x[99] #min(max(z1$x), max(z2$x))
  int <- integrate(difFun, a, b, subdivisions = 200)
  1 - int$value * 0.5
}

d2 <- function (z1, z2) {
  if (! identical(z1$x, z2$x)){stop('different axes')}
  diff <- abs(z1$z - z2$z)
  difFun <- approxfun(z1$x, diff)
  a <- z1$x[2] #max(min(z1$x), min(z2$x))
  b <- z1$x[99] #min(max(z1$x), max(z2$x))
  int <- integrate(difFun, a, b, subdivisions = 200)
  1 - int$value * 0.5
}

# Bhattacharyya distance
bc <- function (z1, z2) {
  fun1 <- approxfun(z1$x, z1$z)
  fun2 <- approxfun(z2$x, z2$z)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))
  a <- z1$x[2] #max(min(z1$x), min(z2$x))
  b <- z1$x[99] #min(max(z1$x), max(z2$x))
  int <- integrate(intgrnd, a, b, subdivisions = 200)
  int$value
}

# Hellinger's H
hinv <- function (z1, z2) {
  fun1 <- approxfun(z1$x, z1$z)
  fun2 <- approxfun(z2$x, z2$z)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))
  a <- z1$x[2] #max(min(z1$x), min(z2$x))
  b <- z1$x[99] #min(max(z1$x), max(z2$x))
  int <- integrate(intgrnd, a, b, subdivisions = 200)
  h <- sqrt(1 - int$value)
  1 - h
}

library(pracma)
sumErr <- function(z1, z2){
  xmin <- min(z1$x)
  xmax <- max(z1$x)
  
  # base integrate function, suppressed errors
  fun1lin <- approxfun(z1$x, z1$z)
  fun2lin <- approxfun(z2$x, z2$z)
  errBaseLin <- 2 - (integrate(fun1lin, xmin, xmax, stop.on.error = FALSE)$value +
                       integrate(fun2lin, xmin, xmax, stop.on.error = FALSE)$value)
  
  # base integrate function, from integrand approximated with spline
  fun1spl <- splinefun(z1$x, z1$z)
  fun2spl <- splinefun(z2$x, z2$z)
  errBaseSpl <- 2 - (integrate(fun1spl, xmin, xmax)$value +
                       integrate(fun2spl, xmin, xmax)$value)
  
  # adaptive Simpson approximation will be used since interval is finite.
  errSimpLin <- 2 - (integral(fun1lin, xmin, xmax) +
                       integral(fun2lin, xmin, xmax))
  errSimpSpl <- 2 - (integral(fun1spl, xmin, xmax) +
                       integral(fun2spl, xmin, xmax))
  
  # trapezoidal integral approximation 
  errTrapLin <- 2 - (trapz(z1$x, z1$z) + trapz(z2$x, z2$z))
  
  cbind(errBaseLin, errBaseSpl, errSimpLin, errSimpSpl, errTrapLin)
}
