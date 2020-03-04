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
# NB: if 'a' and 'b' argumentss given, these should be on the untransformed x scale,
# and one should NOT provide the arguments to/from or lower/upper.
# TODO Add tweak to uniroot function within inverse, for badly behaved empirical data
transformEst <- function(x, w, reflect = FALSE, a=NULL, b=NULL, ...){
  if (is.null(a)){
    test <- c(is.infinite(w(0)), is.na(w(0)))
    if (any(test)){
      a <- min(x)
    } else {
      a <- 0
    }
  }
  
  # transform the observations
  cdf <- function(z){
    integral(w, a, z)
  } 
  Y <- sapply(x, cdf)
  
  # define the boundaries of estimation on the transformed scale
  aYscale <- cdf(a)
  if (is.null(b)){
    b <- max(x)
    bYscale <- max(Y)
  } else {
    bYscale <- cdf(b)
  }
  
  # estimate kernel density, with boundary reflection if specified
  if (reflect){
    kde <- density.reflected(Y, lower = aYscale, upper = bYscale, ...)
  } else {
    kde <- density(Y, from = aYscale, to = bYscale, ...)
  }
  
  # scale to integrate to one
  n <- length(x)
  muHat <- n * sum( w(x)^-1 )^-1
  kde$y <- muHat * kde$y
  
  # Back-transform from the y=W(x) argument to x.
  inv <- inverse(cdf, lower = a, upper = b)
  
  kde$xTrans <- kde$x
  kde$x <- sapply(kde$x, inv)
  
  # if w(x) is undefined at 0, then min and max of X can be off
  # (albeit only slightly) from values at which kernel can estimate 
  f <- approxfun(kde$x, kde$y)
  lwr <- min(kde$x)
  upr <- max(kde$x)
  
  append(list(f=f, lower=lwr, upper=upr), kde)
}

# Borrajo et al. 2017 bandwidth estimation for biased data:
# rule of thumb or bootstrap method
strappy <- function(Y, w, method='brt'){
  if (! method %in% c('brt','rt')){
    stop('method must be brt or rt')
  } 
  
  # estimate of mu sub omega
  n <- length(Y)
  mu <- (sum(1/w(Y))/n)^-1
  
  # estimate of sigma-sq sub omega
  term1 <- sum(w(Y))/n
  sgma <- sqrt(mu*term1 - mu^2)
  
  # estimate of c sub omega
  cHat <- mu * 1/n * sum(1/w(Y)^2)
  
  # estimate rule-of-thumb bandwidth
  K2u <- function(u){
    dnorm(u)^2
  }
  Rk <- integral(K2u, -Inf, Inf)
  u2k <- function(u){
    u^2 * dnorm(u)
  }
  mu2K <- integral(u2k, -Inf, Inf) # 1
  rtNum <- Rk * mu * cHat * 8 * sqrt(pi)
  rtDenom <- n * mu2K * 3
  rt <- (rtNum/rtDenom)^(1/5) * sgma
  
  if (method=='rt'){
    return(rt)
  }
  
  if (method=='brt'){
    # bootstrap estimate of bandwidth, 
    # with rule-of-thumb for pilot
    
    wtsUnscld <- 1/w(Y)
    wts <- wtsUnscld/sum(wtsUnscld)
    fgDens <- density(Y, bw=rt, kernel='gaussian', weights=wts)
    fg <- splinefun(fgDens$x, fgDens$y)
    
    # take 2nd deriv of smooth f estimate, square it, and integrate
    f2g_y <- fg(fgDens$x, deriv=2)
    f2g <- approxfun(fgDens$x, f2g_y)
    f2g_sq <- function(u){
      f2g(u)^2
    } 
    a <- min(fgDens$x)
    b <- max(fgDens$x)
    Rfg <- integral(f2g_sq, a, b)
    
    num <- Rk * mu * cHat
    denom <- n * mu2K * Rfg
    brt <- (num/denom)^(1/5)
    return(brt)
  }
}

# weighted kernel density estimation after Jones 1991
# bw argument can be rt, brt, or a numeric width to use
JonesEst <- function(x, w, bw='brt', reflect=FALSE, a=NULL, b=NULL, ...){
  wts <- 1/w(x)
  wts <- wts/sum(wts)
  
  if (is.null(a)){
    a <- min(x)
  } 
  if (is.null(b)){
    b <- max(x)
  }
  
  # select bandwidth
  if (is.numeric(bw)){
    h <- argmts$bw
  } else {
    if (! bw %in% c('brt','rt')){
      stop('bw method must be brt or rt')
    }
    h <- strappy(x, w, method = bw)
  }
  
  if (reflect){
    kde <- density.reflected(x, bw = h, lower = a, upper = b, weights = wts, ...)
  } else {
    kde <- density(x, bw = h, from = a, to = b, weights = wts, ...)
  }
  
  f <- approxfun(kde$x, kde$y)
  lwr <- min(kde$x)
  upr <- max(kde$x)
  append(list(f=f, lower=lwr, upper=upr), kde)
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
hell <- function (d1, d2, extrap=TRUE) {
  if (extrap==TRUE){
    a <- min(d1$x, d2$x)
    b <- max(d1$x, d2$x)
    
    # ensure the density estimate extends along the full x-axis
    if (min(d1$x) > a){
      d1$x <- c(a, d1$x)
      d1$y <- c(0, d1$y)
    }
    if (min(d2$x) > a){
      d2$x <- c(a, d2$x)
      d2$y <- c(0, d2$y)
    }
    if (max(d1$x) < b){
      d1$x <- c(d1$x, b)
      d1$y <- c(d1$y, 0)
    }
    if (max(d2$x) < b){
      d2$x <- c(d2$x, b)
      d2$y <- c(d2$y, 0)
    }
  }
  
  fun1 <- approxfun(d1$x, d1$y)
  fun2 <- approxfun(d2$x, d2$y)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))
  a <- max(min(d1$x), min(d2$x))
  b <- min(max(d1$x), max(d2$x))
  int <- integral(intgrnd, a, b)
  if (int > 1){
    int <- 1
    warning('estimated Bhattacharyya distance greater than 1')
  }
  h <- sqrt(1 - int)
  h
}
