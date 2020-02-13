library(pracma)
library(GoFKernel)

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

# density.reflected is a function from GoFKernel, and is modified here
# to handle the case where bw arg is given as non-integer, e.g. 'nrd0'
density.reflected <- function (x, lower = -Inf, upper = Inf, weights = NULL, ...) {
  mantener <- !is.na(x)
  x <- x[mantener]
  if (upper < max(x)) 
    warning("There are values in the sample higher than the upper limit")
  if (lower > min(x)) 
    warning("There are values in the sample smaller than the lower limit")
  if (sd(x) == 0) {
    dx <- density(c(x, x[1] + .Machine$double.eps, x[1] - 
                      .Machine$double.eps))
  }
  else {
    if (is.null(weights)) {
      pesos <- rep(1/length(x), length(x))
    }
    else {
      pesos <- weights[mantener]
    }
    argumentos <- list(...)
    
    # GSA revision:
    if ("bw" %in% names(argumentos) & is.numeric(argumentos$bw)) {
      #    if ("bw" %in% names(argumentos)) {
      #      if (is.numeric(argumentos$bw)) 
      broad <- 4 * argumentos$bw
    }
    else {
      pesos <- pesos/sum(pesos)
      broad <- 4 * density(x, weights = pesos, ...)$bw
    }
    if (is.infinite(lower) & is.infinite(upper)) {
      dx <- density(x, weights = pesos, ...)
    }
    else if (is.infinite(lower) & is.finite(upper)) {
      reflected <- which(x >= (upper - broad))
      x.reflect <- c(x, 2 * upper - x[reflected])
      p.reflect <- c(pesos, pesos[reflected])
      p.reflect <- p.reflect/sum(p.reflect)
      dx <- density(x.reflect, weights = p.reflect, ...)
      dx$y <- (dx$y[dx$x >= lower & dx$x <= upper])
      dx$x <- (dx$x[dx$x >= lower & dx$x <= upper])
      bw <- dx$x[2] - dx$x[1]
      area.under <- sum(dx$y) * bw
      dx$y <- dx$y/area.under
    }
    else if (is.finite(lower) & is.infinite(upper)) {
      reflected <- which(x <= (lower + broad))
      x.reflect <- c(x, -x[reflected] + 2 * lower)
      p.reflect <- c(pesos, pesos[reflected])
      p.reflect <- p.reflect/sum(p.reflect)
      dx <- density(x.reflect, weights = p.reflect, ...)
      dx$y <- dx$y[dx$x >= lower & dx$x <= upper]
      dx$x <- dx$x[dx$x >= lower & dx$x <= upper]
      bw <- dx$x[2] - dx$x[1]
      area.under <- sum(dx$y) * bw
      dx$y <- dx$y/area.under
    }
    else {
      reflected.inf <- which(x <= (lower + broad))
      reflected.sup <- which(x >= (upper - broad))
      x.reflect <- c(x, -x[reflected.inf] + 2 * lower)
      p.reflect <- c(pesos, pesos[reflected.inf])
      x.reflect <- c(x.reflect, 2 * upper - x[reflected.sup])
      p.reflect <- c(p.reflect, pesos[reflected.sup])
      p.reflect <- p.reflect/sum(p.reflect)
      dx <- density(x.reflect, weights = p.reflect, ...)
      dx$y <- dx$y[dx$x >= lower & dx$x <= upper]
      dx$x <- dx$x[dx$x >= lower & dx$x <= upper]
      bw <- dx$x[2] - dx$x[1]
      area.under <- sum(dx$y) * bw
      dx$y <- dx$y/area.under
    }
  }
  return(dx)
}

# transformation-based density estimation after Barmi & Simonoff 2000
# TODO
# Figure out why the lower and upper bound of transformed x (from inverse) is too high at big n
# Add tweak to uniroot function within inverse, for badly behaved empirical data
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
# TODO
# Use Borrajo's rule of thumb and 2 bootstrap bandwidth estimators
JonesEst <- function(x, w, reflect=FALSE, ...){
  wts <- 1/w(x)
  wts <- wts/sum(wts)
  
  args <- list(...)
  if ('from' %in% names(args)){
    from <- args$from
  } else {
    from <- min(x)
  }
  if ('to' %in% names(args)){
    to <- args$to
  } else {
    to <- max(x)
  }
  
  if (reflect){
    kde <- density.reflected(x, lower = from, upper = to,
                             weights = wts, ...)
  } else {
    kde <- density(x, from = from, to = to, 
                   weights = wts, ...)
  }
  
  f <- approxfun(kde$x, kde$y)
  a <- min(kde$x)
  b <- max(kde$x)
  append(list(f=f, lower=a, upper=b), kde)
}

# Calculate Holland 1995 niche summary metrics:
# peak abundance, preferred environment, & tolerance
nichStats <- function(d, n=1000){
  pa <- max(d$y)
  pePos <- which.max(d$y)
  pe <- d$x[pePos]
#  f <- approxfun(d$x, d$y)
#  obs <- random.function(n=n, f)
  # Error in integrate(f, lower, upper = z): non-finite function value
#  tol <- sd(obs)
  c(pa=pa, pe=pe)
}

# Hellinger's H
hell <- function (d1, d2) {
  fun1 <- approxfun(d1$x, d1$y)
  fun2 <- approxfun(d2$x, d2$y)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))
  if (! identical(d1$x, d2$x)){warning('different x axes')}
  a <- max(min(d1$x), min(d2$x))
  b <- min(max(d1$x), max(d2$x))
  int <- integral(intgrnd, a, b)
  h <- sqrt(1 - int)
  h
}
