
clipper <- function (t, y, width, n.min, n.sigmas){
  n <- length(t)
  qual <- array(TRUE, dim = n)
  y.mean <-array(NA, dim = n)
  y.sd <- array(NA, dim = n)
  for(i in 1:n) {
    t.0 <- t[i] - width/2
    t.1 <- t[i] + width/2
    mask <- which(t >= t.0 & t < t.1)
    n.sub <- length(mask)
    if(n.sub < n.min) next
    y.sub <- y[mask]
    y.mean[i] <- median(y.sub)
    y.sd[i] <- sd(y.sub)
    y.sub <- y.sub - y.mean[i]
    if (abs(y[i] - y.mean[i]) > n.sigmas * y.sd[i]) {
      qual[i] <- FALSE
    }
  }
  return(list(y.mean=y.mean, y.sd=y.sd, qual=qual))
}  



