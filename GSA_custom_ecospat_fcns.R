library(pracma)

kdeNiche <- 
  function (sp, xmax, xmin, R, h.method="nrd0") {
    sp <- as.matrix(sp)
    l <- list()
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
    l$z <- z$y
    
    # Add Holland 1995 niche summary metrics to output:
    # peak abundance, preferred environment, & tolerance
    l$pa <- max(z$y)
    pe.pos <- which.max(z$y)
    l$pe <- z$x[pe.pos]
    l$t <- sd(sp)
    
    # To approximate the KDE as a function, one also needs x-values:
    l$x <- z$x
    
    return(l)
  }

# Hellinger's H
hell <- function (z1, z2) {
  fun1 <- approxfun(z1$x, z1$z)
  fun2 <- approxfun(z2$x, z2$z)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))
  if (! identical(z1$x, z2$x)){warning('different x axes')}
  a <- max(min(z1$x), min(z2$x))
  b <- min(max(z1$x), max(z2$x))
  int <- integral(intgrnd, a, b)
  h <- sqrt(1 - int)
  h
}
