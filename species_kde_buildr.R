
# Output 1 sp's niche overlap (though time), peak abundance, & preferred enviro
nicher <- function(dat, b1, b2, s, env, xmn, xmx,
                   w1 = NULL, w2 = NULL, reflect = FALSE, ...){
  
  sp1rows <- which(dat$sp$species==s & dat$sp$bin==b1)
  sp1 <- dat$sp[sp1rows,envNm]
  
  sp2rows <- which(dat$sp$species==s & dat$sp$bin==b2)
  sp2 <- dat$sp[sp2rows,envNm]
  
  noWeight <- any(is.null(w1), is.null(w2))
  if (! noWeight){
    d1 <- tryCatch(
      transformEst(sp1, w = w1, bw='SJ-ste', reflect = reflect, a = xmn, b = xmx, ...),
      #      JonesEst(sp1, w = w1, bw='brt', reflect = reflect, a = xmn, b = xmx, ...),
      error = function(err){ list() }
    ) 
    d2 <- tryCatch(
      transformEst(sp2, w = w2, bw='SJ-ste', reflect = reflect, a = xmn, b = xmx, ...),
      #      JonesEst(sp2, w = w2, bw='brt', reflect = reflect, a = xmn, b = xmx, ...),
      error = function(err){ list() }
    ) 
  } else {
    # use unweighted sampling, either regular density or reflected:
    if (reflect){
      d1 <- tryCatch(
        density.reflected(sp1, bw='SJ-ste', lower = xmn, upper = xmx, ...),
      )
      d2 <- tryCatch(
        density.reflected(sp2, bw='SJ-ste', lower = xmn, upper = xmx, ...),
      )
    } else {
      d1 <- tryCatch(
        density(sp1, bw='SJ-ste', from = xmn, to = xmx, ...),
      )
      d2 <- tryCatch(
        density(sp2, bw='SJ-ste', from = xmn, to = xmx, ...),
      )
    }
  }
  
  # the species may be absent in one or both bins, in which case d is an empty list
  if (length(d1)==0){
    data.frame(bin=NA, bin2=NA, sp=NA, h=NA, pa=NA, pe=NA)
    
  } else {
    stats <- nichStats(d1)
    
    if (length(d2)==0){
      data.frame(bin=b1, bin2=b2, sp=s, h=NA, t(stats))
    } else{
      h <- hell(d1, d2, extrap = TRUE) 
      data.frame(bin=b1, bin2=b2, sp=s, h=h, t(stats))
    }
  }
}

kde <- function(dat, bPair, envNm){
  b1 <- bPair[1]
  b2 <- bPair[2]
  xmn <- min(dat$samp[,envNm])
  xmx <- max(dat$samp[,envNm])
  
  # estimate bias function for each time bin based on sampling distribution
  sampRows1 <- which(dat$samp$bin==b1)
  samp1 <- dat$samp[sampRows1,envNm]
  # Reflecting the sample curve doesn't change the sp KDE much
  # except that the ends turn down a bit more (more convexity).
  # Since it's more complicated and throws warnings, don't do it.
  densSamp1 <- density(samp1, bw='SJ-ste', from=xmn, to=xmx)
  #  densSamp1 <- density.reflected(samp1, lower=xmn, upper=xmx) 
  w1 <- approxfun(densSamp1$x, densSamp1$y)
  
  # in the most recent time bin, there is no subsequent bin
  if (is.na(b2)){
    w2 <- NA
    
    # redefine axis limits - discretization of samp KDE can shrink them a bit
    xmnNew <- min(densSamp1$x)
    xmxNew <- max(densSamp1$x)
  } else {
    sampRows2 <- which(dat$samp$bin==b2)
    samp2 <- dat$samp[sampRows2,envNm]
    densSamp2 <- density(samp2, bw='SJ-ste', from=xmn, to=xmx)
    #  densSamp2 <- density.reflected(samp2, lower=xmn, upper=xmx) 
    w2 <- approxfun(densSamp2$x, densSamp2$y)
    
    # redefine axis limits - discretization of samp KDE can shrink them a bit
#    lim1 <- range(densSamp1$x)
#    lim2 <- range(densSamp2$x)
#    xmnNew <- max(lim1[1], lim2[1])
#    xmxNew <- min(lim1[2], lim2[2])
  } 
#  if (sum(identical(xmn, xmnNew), identical(xmx, xmxNew))!=2){
#    error('code needed after all')
#  }
  zoneSp <- unique(dat$sp$species)
  sList <- lapply(zoneSp, function(s){
    nicher(dat = dat, b1 = b1, b2 = b2, s = s,
           w1 = w1, w2 = w2, xmn = xmn, xmx = xmx)
  })
  do.call(rbind, sList)
}
