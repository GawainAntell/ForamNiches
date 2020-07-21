library(kerneval)

# get peak abundance and preferred environment
nichStats <- function(dens){
  pa <- max(dens$y)
  pePos <- which.max(dens$y)
  pe <- dens$x[pePos]
  c(pa = pa, pe = pe)
}

# Output 1 sp's niche overlap (though time), peak abundance, & preferred enviro
nicher <- function(dat, b1, b2, s, env, xmn, xmx,
                   w1 = NULL, w2 = NULL, reflect = FALSE, ...){
  
  sp1rows <- which(dat$sp$species == s & dat$sp$bin == b1)
  sp2rows <- which(dat$sp$species == s & dat$sp$bin == b2)
  sp1 <- dat$sp[sp1rows, envNm]
  sp2 <- dat$sp[sp2rows, envNm]
  
  d1 <- tryCatch(
    transdens(sp1, w = w1, bw = 'nrd0', reflect = reflect, a = xmn, b = xmx, ...),
    #     bw = 'SJ-ste'
    error = function(err){ list() }
  ) 
  d2 <- tryCatch(
    transdens(sp2, w = w2, bw = 'nrd0', reflect = reflect, a = xmn, b = xmx, ...),
    error = function(err){ list() }
  ) 
  
  # The species may be absent in one or both bins. If absent, d is an empty list
  if (length(d1) == 0){
    data.frame(bin = NA, bin2 = NA, sp = NA, h = NA, pa = NA, pe = NA, 
               bw1 = NA, bw2 = NA, n1 = NA, n2 = NA)
  } else {
    d1deets <- nichStats(d1)
    
    if (length(d2) == 0){
      data.frame(bin = b1, bin2 = b2, sp = s, h = NA, t(d1deets),
                 bw1 = d1$bw, bw2 = NA, n1 = length(sp1), n2 = length(sp2))
    } else{
      h <- hell(d1, d2) 
      data.frame(bin = b1, bin2 = b2, sp = s, h = h, t(d1deets),
                 bw1 = d1$bw, bw2 = d2$bw, n1 = length(sp1), n2 = length(sp2))
    }
  }
}

kde <- function(dat, bPair, envNm){
  b1 <- bPair[1]
  b2 <- bPair[2]
  xmn <- min(dat$samp[,envNm])
  xmx <- max(dat$samp[,envNm])
  
  # estimate bias function for each time bin based on sampling distribution
  sampRows1 <- which(dat$samp$bin == b1)
  samp1 <- dat$samp[sampRows1, envNm]
  # Reflecting the sample curve doesn't change the sp KDE much
  # except that the ends turn down a bit more (more convexity).
  # Since it's more complicated, don't do it.
  densSamp1 <- density(samp1, bw='nrd0', from=xmn, to=xmx)
  w1 <- approxfun(densSamp1$x, densSamp1$y)
  
  # estimate bias function for the younger time bin
  if (is.na(b2)){
    # in the most recent time bin (4 ka), there is no subsequent bin
    w2 <- NA
  } else {
    sampRows2 <- which(dat$samp$bin == b2)
    samp2 <- dat$samp[sampRows2, envNm]
    densSamp2 <- density(samp2, bw = 'nrd0', from = xmn, to = xmx)
    w2 <- approxfun(densSamp2$x, densSamp2$y)
  } 
  
  zoneSp <- unique(dat$sp$species)
  sList <- lapply(zoneSp, function(s){
    nicher(dat = dat, b1 = b1, b2 = b2, s = s,
           w1 = w1, w2 = w2, xmn = xmn, xmx = xmx)
  })
  do.call(rbind, sList)
}
