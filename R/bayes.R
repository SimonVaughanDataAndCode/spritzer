
# ------------------------------------------------
# define PSD model
#
# Additional parameters that are not to be fitted but
# passed to the model functions are stored in the MOD
# array.

 model <- function(par, x, mod=c(0, 1.0)) {

   a.low <- mod[2]

   if (mod[1] == 0) {
     y.mod <- model0(par, x)
   }

   if (mod[1] == 1) {
     y.mod <- model1(par, x, a.low)
   }

   if (mod[1] == 2) {
	y.mod <- model2(par, x)
   }

     if (mod[1] == 3) {
	y.mod <- model3(par, x, a.low)
   }

 return(y.mod)
}

# ------------------------------------------------
# define PSD model as a power law y = b*x^-a + c

model0 <- function(par, x) {
  logmin <- -100
  a <- par[1]
  b <- par[2]
  c <- par[3]

  ly.mod <- -a*log(x)+b
  ly.mod <- pmax(ly.mod,logmin)
  y.mod <- exp(ly.mod) + exp(c)
  return(y.mod)
}

# ------------------------------------------------
# define PSD model as a bending power law
#
# y(x) = b*x^(-a1) / [1 + (x/x_b)^(a2-a1) ] + c
#
# y(x >> x_b) = A * x^(-a2) * x_b^(a2-a1)
# y(x << x_b) = A * x^(-a1)
#
# This is a power law of slope -a1 at x << x_b
# and a power law of slope -a2 at x >> x_b.
# Between these two extremes is smoothly changes
# from one slope to the other.
# (x_b is the break/bend point.)
#
# The parameters are stored in the one dimensional array par
# par = [a2, x_b, b, c]
#      a1  = power law index at x << x_b [given seperately]
#      a2  = power law index at x >> x_b
#      x_b = break point
#      b   = normalisation (input as log[N])
#      c   = constant
#
# x_b, b and c are passed to the function in log_10 units.

model1 <- function(par, x, a.low=1.0) {

  logmin <- -100

  a1 <- a.low
  a2 <- par[1]
  x_b <- par[2]
  b <- par[3]
  c <- par[4]

# convert x to log units

  logx <- log(x)

# set up arrays

  logq <- x*0
  y <- x*0

# Calculate bending factor q=1+z with z = (x/x_b)^(a2-a1).
# But to avoid overflow/underflow we calculate
#   log[z] = (a2 - a1) * (log[x] - log[x_b])

  logz <- (a2-a1)*(logx - x_b)

# split into regions of low, medium and high z.
# i.e. where x << x_b, x ~ x_b and x >> x_b
# and treat these seperately to avoid underflow/
# overflow errors

  lo <- (logz < -16)
  hi <- (logz > 16)
  me <- (logz > -16) & (logz < 16)

  if (sum(lo, na.rm=TRUE) > 0) { logq[lo] <- log(1.0) }
  if (sum(hi, na.rm=TRUE) > 0) { logq[hi] <- logz[hi] }
  if (sum(me, na.rm=TRUE) > 0) { logq[me] <- log(exp(logz[me])+1) }

# calculate log(y)

  logy <- -a1*logx - logq + b

# watch for very low/high values (catch over/underflow)

  lo <- (logy < logmin)
  hi <- (logy > -logmin)
  me <- (logy > logmin) & (logz < -logmin)

  if (sum(hi, na.rm=TRUE) > 0) { y[hi] <- exp(-logmin) }
  if (sum(me, na.rm=TRUE) > 0) { y[me] <- exp(logy[me]) }

  y <- y + exp(c)

  return(y)

}


# ------------------------------------------------
# define PSD model as a power law y = b*x^-a

model2 <- function(par, x) {
  logmin <- -100
  a <- par[1]
  b <- par[2]

  ly.mod <- -a*log(x)+b
  ly.mod <- pmax(ly.mod,logmin)
  y.mod <- exp(ly.mod)
  return(y.mod)
}

# ------------------------------------------------
# define PSD model as a bending power law
#
# y(x) = b*x^(-a1) / [1 + (x/x_b)^(a2-a1) ]
#
# y(x >> x_b) = A * x^(-a2) * x_b^(a2-a1)
# y(x << x_b) = A * x^(-a1)
#
# This is a power law of slope -a1 at x << x_b
# and a power law of slope -a2 at x >> x_b.
# Between these two extremes is smoothly changes
# from one slope to the other.
# (x_b is the break/bend point.)
#
# The parameters are stored in the one dimensional array par
# par = [a2, x_b, b]
#      a1  = power law index at x << x_b [given seperately]
#      a2  = power law index at x >> x_b
#      x_b = break point
#      b   = normalisation (input as log[N])
#
# x_b, b are passed to the function in log_10 units.

model3 <- function(par, x, a.low=1.0) {

  logmin <- -100

  a1 <- a.low
  a2 <- par[1]
  x_b <- par[2]
  b <- par[3]

# convert x to log units

  logx <- log(x)

# set up arrays

  logq <- x*0
  y <- x*0

# Calculate bending factor q=1+z with z = (x/x_b)^(a2-a1).
# But to avoid overflow/underflow we calculate
#   log[z] = (a2 - a1) * (log[x] - log[x_b])

  logz <- (a2-a1)*(logx - x_b)

# split into regions of low, medium and high z.
# i.e. where x << x_b, x ~ x_b and x >> x_b
# and treat these seperately to avoid underflow/
# overflow errors

  lo <- (logz < -16)
  hi <- (logz > 16)
  me <- (logz > -16) & (logz < 16)

  if (sum(lo, na.rm=TRUE) > 0) { logq[lo] <- log(1.0) }
  if (sum(hi, na.rm=TRUE) > 0) { logq[hi] <- logz[hi] }
  if (sum(me, na.rm=TRUE) > 0) { logq[me] <- log(exp(logz[me])+1) }

# calculate log(y)

  logy <- -a1*logx - logq + b

# watch for very low/high values (catch over/underflow)

  lo <- (logy < logmin)
  hi <- (logy > -logmin)
  me <- (logy > logmin) & (logz < -logmin)

  if (sum(hi, na.rm=TRUE) > 0) { y[hi] <- exp(-logmin) }
  if (sum(me, na.rm=TRUE) > 0) { y[me] <- exp(logy[me]) }

  return(y)

}


# ------------------------------------------------
# define periodogram log likelihood function
# S = -log(likelihood)

mlogl <- function(par, x, y, mod=c(0, 1.0)) {
  mody <- model(par, x, mod=mod)
  l <- sum( log(mody) + y/mody , na.rm=TRUE)
  return(l)
}

# ------------------------------------------------
# define posterior in terms of likelihood * prior
# for the parameters. Actually we work with
# -log[posterior] = -log[likelihood] + -log[prior]

lpost <- function(par, x, y, mod=c(0, 1.0)) {
  ml.post <- mlogl(par, x, y, mod=mod) + mlprior(par, mod=mod)
  return(ml.post)
}


# ------------------------------------------------
# Define the prior densities for the parameter
# of the model. Actually, calculate the minus
# log prior density, which can then be combined
# with the minus log likelihood (MLOGL).

mlprior <- function(par, mod=c(0, 1.0)) {

# smallest allowed log(number) - to avoid underflow
# and -log(0) infinities.

  logmin <- -100

# set allowed range of parameters

  alim <- c(-1,8)

# define prior density at parameters = par
# as the product of the prior densities
# for each parameter

  if ( mod[1] == 0 ) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0
     pr <- pa*pb*pc
  }

  if (mod[1] == 1) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     d <- par[4]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0
     pd <- 1.0

# Alternative priors using Normal densities (for RE J1034)
#     pa <- dnorm(a,mean=2,sd=2)
#     pb <- dnorm(b,mean=-6.9,sd=2.3)
#     pc <- dnorm(c,mean=-4.6,sd=2.3)
#     pd <- dnorm(d,mean=0,sd=2.3)

     pr <- pa*pb*pc*pd
  }

  if (mod[1] == 2) {
     a <- par[1]
     b <- par[2]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pr <- pa*pb
  }

  if (mod[1] == 3) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0

     pr <- pa*pb*pc
  }


# take the log, watching out for zeros

  if (pr > 0) {
    mlp <- -log(pr)
  } else {
    mlp <- -logmin
  }

  return(mlp)

}

# ------------------------------------------------
# function to compute the periodogram from an input
# time series

pgram <- function(x, dt, rms=TRUE, leahy=FALSE, plotq=TRUE) {
  N <- length(x)
  df <- 1/(N*dt)
  f <- seq(1,N/2)*df
  x.mean <- mean(x)

# define the different normalisation options

  norm <- 2.0*dt/N
  if (rms == TRUE && leahy == FALSE) {
    norm <- norm/x.mean^2
    ylab <- "Power density [(rms/mean)^2/Hz]"
  }
  if (leahy == TRUE) {
    norm <- norm/x.mean
    ylab <- "Power density [Leahy norm.]"
  }
  if (leahy == FALSE && rms == FALSE) {
    ylab <- "Power density [abs. norm.]"
  }
  xlab <- "Frequency"

# calculate the FFT, square it and normalise

  p <- fft(x)
  pow <- abs(p[2:(N/2+1)])
  pow <- pow * pow * norm

# make a plot if requested

  if (plotq == TRUE) {
    plot(f, pow, bty="l", type="s", log="xy", bty="l")

    n.filt <- 10
    filt <- array(1, dim=n.filt)/n.filt
    filter.y <- filter(log(pow), filter=filt)
    filter.y <- exp(filter.y + 0.577214)
    lines(f, filter.y, type="l", col="red", lw=5)

    filt <- c(0.25, 0.5, 0.25)
    filter.y <- filter(pow, filter=filt)
    lines(f, filter.y, type="s", col="green", lw=3)

    smooth.per <- lowess(log(f), log(pow), f=0.1)
    smooth.x <- exp(smooth.per$x)
    smooth.y <- exp(smooth.per$y + 0.557214)
    lines(smooth.x, smooth.y, type="l", col="blue", lw=5)
  }

# put the results in a list to return

  per <- list(f=f, p=pow, ylab=ylab, xlab=xlab, p.smooth=filter.y)

  return(per)
}

# ------------------------------------------------
# The MCMC routine. Metropolis-Hastings algorithm
# using a Normal random walk for the proposal
# distribution.
#
#  theta.0 - vector of starting parameter values
#  cov     - covariance matrix for proposal dist.
#  M       - number of parameters in theta
#  N       - number of iterations to run
#  accept  - counter for number of acceptances

mchain <- function(x, y, N, mod=c(0, 1.0), theta.0, cov, discard=floor(N/2)) {

# load the MNORMT package (for multivariate normals RMNORM)

  require(mnormt)

# set up output arrays

  M <- length(theta.0)
  theta <- array(0, dim=c(M,N))
  log.p <- array(0, dim=N)
  accept <- 0

# starting location of chain

  theta[,1] <- theta.0
  log.p[1] <- lpost(theta[,1], x, y, mod=mod)

# loop over the chain

  for (t in 2:N) {

# draw a value from the proposal distribution

#    theta.prop <- rmnorm(1, theta[,t-1], cov)

# or use alternative proposal (Student's T distribution)

    theta.prop <- rmt(1, theta[,t-1], cov/3, df=3)

# compute ratio of posteriors at new and old locations
#   r = p(theta.prop|data)/p(theta[t-1]|data)
# in terms of the minus log posterior function lpost
# this is log[p(theta.old|data)] - log[p(theta.new|data)]

  log.pnew <- lpost(theta.prop, x, y, mod=mod)
  log.r <- log.p[t-1] - log.pnew

# decide whether or not to update theta.
# theta is updated with probability min(r,1)
# otherwise it is left as before.

  log.r <- min(log.r,0)
  r <- exp(log.r)
  update <- sample(c(TRUE,FALSE), size=1 ,prob=c(r,1-r))
  if (update) {
    theta[,t] <- theta.prop
    log.p[t] <- log.pnew
    if (t > discard) { accept <- accept+1 }
  } else {
    theta[,t] <- theta[,t-1]
    log.p[t] <- log.p[t-1]
  }

# end of loop over the chain

# uncomment this line to output iterations in real time
#  cat("-->", t, theta[,t], r, log.p[t], update, fill=TRUE)

  }

# once the chain has run we discard the first
# section - to remove memory of the starting point

  theta <- theta[,(discard+1):N]
  log.p <- log.p[(discard+1):N]
  L <- N - discard
  accept <- accept/L

  mcmc.out <- list(theta=theta,
                   accept=accept,
                   log.p=log.p, L=L)

  return(mcmc.out)
}

# ------------------------------------------------
# Produce some diagnostic plots of MCMC output
# Specifically, plot the time series, autocorrelation
# and histogram for each parameter of each chain

  mcmc.diag <- function(mcmc,theta.map,theta.err,j) {

  M <- dim(mcmc$theta)[1]

  cat(" acceptance rate: ",mcmc$accept,fill=TRUE)

  par(mfrow=c(M,3))

  txt <- paste("Chain",j)
  for (m in 1:M) {

    lab <- paste("theta[",m,"]", sep="")
    plot(mcmc$theta[m,], type="l", ylab=lab, xlab="t", bty="n")
    if (m == 1) { title(txt) }

    hist(mcmc$theta[m,], prob=TRUE, col="blue", xlab=lab, main="")
    xlim <- c(-3,3)*theta.err[m]+theta.map[m]
    x <- seq(xlim[1], xlim[2], length.out=100)
    lines(x, dnorm(x, mean=theta.map[m], sd=theta.err[m]))

    acf(mcmc$theta[m,], lag.max=30, main=lab)

    }

  }

# ------------------------------------------------
# for each parameter theta[1]...theta[M]
# calculate the R.hat statistic (Gelman & Rubin 1992)
# See also Gelman et al. (2004, sect 11.6)

  Rhat <- function(theta, M, L) {

    R.hat <- array(0, dim=M)
    for (m in 1:M) {
      tj <- theta[m,,]
      sj <- apply(tj, 2, var)
      W <- mean(sj)
      mj <- apply(tj, 2, mean)
      B <- var(mj)*L
      v <- (L-1)/L*W + B/L
      R.hat[m] <- sqrt(v/W)

      cat("-- Theta[",m,"] R.hat = ",R.hat[m],sep="")
      if (R.hat[m] > 1.2) {
        cat(" ** High R.hat. Check results. **",fill=TRUE)
      } else {
        cat(" -- Good R.hat.", fill=TRUE)
      }
    }
    return(R.hat)
  }

# ------------------------------------------------
# The master MCMC routine

  mcmc <- function(x, y, N, J, mod=c(0, 1.0), theta.map, cov, discard=floor(N/2)) {

  M <- length(theta.map)
  L <- N - discard
  theta.err <- sqrt(diag(cov))
  theta <- array(0, dim=c(M, L, J))

  cat("----------------------------------------- \n")
  cat("-- MCMC analysis \n")
  cat("-- ",J," chains, each with ",N," iterations (first ",
          discard," discarded) \n", sep="")
  cat("-- n.sims = ", J*L, " iterations saved \n")

# loop over multiple chains j=1,...,J

  for (j in 1:J) {

     theta.0 <- theta.map + sample(c(2,3,-3,-2),replace=TRUE,size=M)*theta.err
     mc.out <- mchain(x, y, N, mod=mod, theta.0, cov, discard)
     L <- mc.out$L
     cat("-- Chain",j)

# run diagnostic plots for each chain

     mcmc.diag(mc.out, theta.map, theta.err, j)

# then collate the chains into one array

     theta[,,j] <- mc.out$theta
  }

# perform checks for convergence of the J chains

  mcmc.conv(theta, M, L, J)

# collate the results from J chains into one array
# return the MCMC draws to the calling routine

  dim(theta) <- c(M, J*L)

# transpose the array from an M*N matrix to N*M
# where M = no parameters and N = no simulations.

  theta <- t(theta)

# label the columns

  colnames(theta) <- paste("theta[",1:M,"]", sep="")

  return(theta)

}

# ------------------------------------------------
# Perform checks for convergence of multiple
# Markov chains. Uses Gelman & Rubin's R.hat
# for each parameter, and also a visual check of
# the 80% regions for each parameter.

  mcmc.conv <- function(theta, M, L, J) {

  layout(t(c(1,2)), widths=c(0.3,0.7))

# Calculate R.hat (Gelman & Rubin 1992) for each
# parameter as a check for convergence of chains

  R.hat <- Rhat(theta, M, L)

# plot the R.hat results

  plot(R.hat, 1:M, xlim=c(1,2), ylim=c(0.5,M+0.5), bty="n", pch=16, yaxp=c(1,M,M-1),
  xlab="R.hat", ylab="Parameter")
  abline(v=1.1, lty=2)
  axis(3)

# for each parameter theta[1]...theta[M]
# calculate and plot the 80% intervals from each chain

  plot(rep(0,M+1), 1:(M+1), ylim=c(0.5,M+0.5), xlim=c(-2,2), type="n", bty="n",
    ylab="Parameter", xlab="80% region (scaled)", yaxp=c(1,M,M-1))
  axis(3)

  for (m in 1:M) {
    tj <- theta[m,,]
    intv <- apply(tj, 2, quantile, prob=c(0.1,0.9))
    ci0 <- intv[1,]
    ci1 <- intv[2,]
    mt <- apply(tj, 2, mean)
    scale <- mean(ci1-ci0)/2
    offset <- mean(mt)
    ci0 <- (ci0 - offset)/scale
    ci1 <- (ci1 - offset)/scale
    for (j in 1:J) {
      x <- m + (j-1)/J/4
      segments(ci0[j], x,ci1[j], x, col=j)
    }
  }

  }



# ------------------------------------------------

infer <- function(theta) {

  cat("----------------------------------------- \n")
  cat("-- MCMC inferences \n")

  N <- dim(theta)[1]
  M <- dim(theta)[2]

# compute the covariance matrix from the simulated draws

  covar <- cov(theta)
  cat("-- Covariance matrix (from simulations) \n")

# print on screen

  print(covar)

# calculate for each parameter its mean and
# equal tail 90% interval. These are the posterior
# means and 90% credible intervals from the MCMC.

  theta.mean <- apply(theta, 2, mean)
  theta.sd <- apply(theta, 2, sd)
  theta.ci <- apply(theta, 2, quantile, prob=c(0.05,0.95))

# collect the output into one array

  result <- cbind(theta.mean, theta.sd, t(theta.ci))

# put some names on the rows/columns to make it
# easier to read.

  colnames(result) <- c("mean", "sd", "5%", "95%")
  rownames(result) <- colnames(theta)

# print the output on screen

  cat(" \n")
  cat("-- Posterior summaries of parameters \n")
  print(result)

# 2D correlation plots of the parameters
#  pairs(theta, gap=0, labels=colnames(theta), pch=".")

#  for (i in 1:M) { labels[i] <- expression(theta[i]) }
  cont.pairs(theta)

# return the output to the calling routine

  return(result)

  }


# ------------------------------------------------
# Model fitting function
# uses NLM to minimise the minus log posterior
# (or likelihood).

  fit.model <- function(x, y, mod=c(0, 1.0), par.0, obs=FALSE, ylab="Power") {

# define the size, resolution of the data
# and number of model parameters

    nx <- length(x)
    dx <- x[1]
    M <- length(par.0)

# before minimising the fit statistic we make a crude
# estimate of the normalisation in order that this is
# of the right order to begin with.
# We renormalise the power law by taking the ratio of
# the observed power to the model power.

   var.obs <- sum(y)
   var.mod <- sum(model(par.0, x, mod=mod))
   renorm <- (var.obs/var.mod)
   if (mod[1] == 0) { par.0[2] <- par.0[2] + log(renorm) }
   if (mod[1] == 1) { par.0[3] <- par.0[3] + log(renorm) }
   if (mod[1] == 2) { par.0[2] <- par.0[2] + log(renorm) }
   if (mod[1] == 3) { par.0[3] <- par.0[3] + log(renorm) }

# minimise the function lpost starting from position
# par.0. If this is the real data (obs=TRUE) then
# also calculate the hessian matrix (second derivatives
# of lpost). This can be inverted to get the covariance
# matrix and the mode.

    result <- nlm(lpost, par.0, hessian=obs, x=x, y=y, mod=mod)

# Store the result - MAP = maximum a posteriori

    theta.map <- result$estimate

# compute the model at the posterior mode

    y.map <- model(theta.map, x, mod=mod)

# compute the deviance (-2 * log likelihood)

    D.obs = 2.0*mlogl(theta.map,x,y,mod=mod)

# and highest outlier R.max = max(2*data/model)

    rat <- 2.0*(y/y.map)
    T.obs <- max(rat)
    j.peak <- which.max(rat)
    f.peak <- x[j.peak]

# compute the KS statistic comparing the data/mode
# residuals to an exponential distribution (~chi^2_2/2)

    ks.out <- ks.test(rat/2, "pexp")
    ks <- ks.out$statistic

# highest outlier from smoothed data.
# i.e. smooth raw periodogram with a filter
# and take ratio of smoothed data to fitted model

    kern <- c(0.25, 0.5, 0.25)
    y.smooth <- filter(y, filter=kern)
    rat.smooth <- 2.0*(y.smooth/y.map)
    srat <- max(rat.smooth, na.rm=TRUE)
    sj.peak <- which.max(rat.smooth)
    sf.peak <- x[sj.peak]

# compute the sum of the square standardised residuals
# (like chi square). See Anderson et al. (1990; ApJ, 364, 699)

    sse <- sum(((y-y.map)/y.map)^2)

# If obs == TRUE (i.e. this is the "real" observation)
# then present the results of the model fitting to the
# user. Otherwise (e.g. if this is a simulation) remain
# silent.

    if (obs == TRUE) {

# plot the MAP model and data/model residuals
# save the existing plotting parameters,
# to be restored once finished

      par.mfrow <- par()$mfrow
      par.mfcol <- par()$mfcol
      par.mar <- par()$mar
      par.oma <- par()$oma

      layout(matrix(c(1,2)), heights=c(1.5,1))
      par(mar=c(0,6,0,0),oma=c(0,2,2,2))

      x.plot <- c(x[1]-dx/2,x)
      y.plot <- model(theta.map, x.plot, mod=mod)
      plot(x-dx/2, y, log="xy", type="s", xlab="", ylab=ylab, xaxt="n", bty="l")
      axis(1, labels=FALSE)
      lines(x.plot, y.plot, lty=1, lwd=4, col="red")

      par(mar=c(6,6,0,0))
      plot(x-dx/2, rat/2, log="xy", type="s", xlab="Frequency", ylab="data/model", bty="l")
      lines(range(x)-dx/2,c(1,1), lwd=4, col="red")

# restore the graphics device parameters

      par(mfcol=par.mfcol, mfrow=par.mfrow, oma=par.oma, mar=par.mar)

# calculate errors based on Normal approximation

      theta.cov <- solve(result$hessian)
      theta.err <- sqrt(diag(theta.cov))

# present results on screen

      cat("----------------------------------------- \n")
      cat("-- Model ",mod," result: \n", sep="")
      for (i in 1:M) {
        cat("--   theta[",i,"] = ", theta.map[i]," +/- ", theta.err[i],
          sep="",fill=TRUE)
      }

# minimisation diagnostics

      if (result$code == 3) {
        cat ("** Last global step in NLM failed to locate a point
          lower than estimate.\n")
      }
      if (result$code == 4) {
        cat("** Iteration limit exceeded in NLM. \n")
      }
      if (result$code == 5) {
        cat("** Maximum stepsize in NLM exceeded five times. \n")
      }

# display the deviance and highest outlier from observation

      cat(" ",fill=TRUE)
      cat("-- Number of frequencies = ", nx, fill=TRUE)
      cat("-- Deviance [-2 log L] D = ", D.obs, fill=TRUE)
      cat("-- Highest data/model outlier 2I/S = ", T.obs)
      cat(" at frequency ", f.peak, fill=TRUE)
      cat("-- Highest smoothed data/model outlier 2Sm[I]/S = ", srat)
      cat(" at frequency ", sf.peak, fill=TRUE)

# calculate the expected summed residual conditional on the model
# begin correct: S=sum(2*data/model) --> E[S] = 2N and V[S] = 4N
# since 2*data/model is chi^2_2 per data point.

      S.obs <- sum(rat)
      S.exp <- 2.0*(nx-M)
      S.sd <- sqrt(4.0*(nx-M))

      cat("-- Summed residual S = ", S.obs, sep="", fill=TRUE)
      cat("-- Expected S ~ ", S.exp," +/- ", S.sd, sep="", fill=TRUE)

      cat("-- KS test p-value (use with caution) p = ", ks.out$p.value,
           sep="", fill=TRUE)
      cat("-- Merit function (SSE) = ", sse, fill=TRUE,sep="")

    } else {

      theta.cov <- 0

    }

# include the fitted parameters, covariance and test
# statistics in the output list

  result <- list(theta.map=theta.map, theta.cov=theta.cov,
                 T=T.obs, D=D.obs, f.peak=j.peak, ks=ks, sse=sse,
                 y.fit=y.map, srat=srat)

  return(result)

  }

# ------------------------------------------------
  pp <- function(theta, x, mod=c(0, 1.0), Q=dim(theta)[2], x0) {

# use the MCMC draws from the posterior to perform
# posterior predictive checks.
# Randomly draw parameters theta.sim, use these to
# produce some simulated data y.sim. Fit the
# data and record the test quantities (as with the true data).

# randomly shuffle the MCMC output to mix the
# individual chains

    L <- dim(theta)[1]
    M <- dim(theta)[2]
    theta <- theta[sample(1:L),]

# Q = number of simulations (Q <= L)

    Q <- min(Q, L)

# initialise arrays for storing the test quantities

    D <- array(0, dim=Q)
    T <- D
    f.peak <- D
    ks <- D
    sse <- D
    y0 <- D
    lrt <- D
    srat <- D

    nx <- length(x)

# loop over Q simulations

    cat("----------------------------------------- \n")
    cat("-- Posterior predictive checks with ",Q," simulations",sep="",fill=TRUE)

    jump <- floor(Q/10)
    for (sim in 1:Q) {

      theta.sim <- theta[sim,]                # draw random paras
      y.sim <- model(theta.sim, x, mod=mod)   # 'true' spectrum
      y.sim <- y.sim*rexp(nx)                 # random data

      fit.sim <- fit.model(x, y.sim, mod=mod, theta.sim) # fit data

# store output

      D[sim] <- fit.sim$D
      T[sim] <- fit.sim$T
      ks[sim] <- fit.sim$ks
      sse[sim] <- fit.sim$sse
      f.peak[sim] <- fit.sim$f.peak
      y0[sim] <- y.sim[x0]
      srat[sim] <- fit.sim$srat

      lrt[sim] <- calc.lrt(x, y.sim, mod=mod, theta.sim, D[sim])

      if (sim %% jump == 0) { cat("-- ",sim/jump*10,"%",sep="",fill=TRUE) }
    }

    pp.out <- list(T.sim=T, D.sim=D, f.sim=f.peak, ks.sim=ks, sse.sim=sse,
                   y0=y0, lrt.sim=lrt, srat.sim=srat)
    return(pp.out)

  }

# ------------------------------------------------
# Function to calculate the Likelihood Ratio Test
# statistic (LRT). Two (nested?) models are fitted
# to the same data, and the difference between the two
# -2*log[likelihood] functions is LRT.
#
# If original model is 0 (or 1) use 1 (or 0) as alternative
# If original model is 2 (or 3) use 3 (or 2) as alternative

calc.lrt <- function(x, y, mod, theta, D) {

# guestimate a bend frequency - in middle of log[freq] range

      f.br <- 0.5*(max(log(x)) + min(log(x)))

      if (mod[1] == 0) {
	   mod.alt <- c(1, mod[2])
         par0 <- c(theta[1], f.br, theta[2], theta[3])
         y.mod <- model(par0, x, mod=mod.alt)
         par0[3] <- par0[3] + log(sum(y)/sum(y.mod))
         fit.alt <- nlm(mlogl, par0, x=x, y=y, mod=mod.alt, fscale=D)
       }

	if (mod[1] == 1) {
	   mod.alt <- c(0, mod[2])
         par0 <- c(theta[1], theta[3], theta[4])
         y.mod <- model(par0, x, mod=mod.alt)
         par0[2] <- par0[2] + log(sum(y)/sum(y.mod))
         fit.alt <- nlm(mlogl, par0, x=x, y=y, mod=mod.alt, fscale=D)
      }

      if (mod[1] == 2) {
	   mod.alt <- c(3, mod[2])
         par0 <- c(theta[1], f.br, theta[2])
         y.mod <- model(par0, x, mod=mod.alt)
         par0[3] <- par0[3] + log(sum(y)/sum(y.mod))
         fit.alt <- nlm(mlogl, par0, x=x, y=y, mod=mod.alt, fscale=D)
      }

	if (mod[1] == 3) {
	   mod.alt <- c(2, mod[2])
         par0 <- c(theta[1], theta[3])
         y.mod <- model(par0, x, mod=mod.alt)
         par0[2] <- par0[2] + log(sum(y)/sum(y.mod))
         fit.alt <- nlm(mlogl, par0, x=x, y=y, mod=mod.alt, fscale=D)
      }

     D.alt <- 2.0*fit.alt$minimum

# calculate LRT = different between the -2*log[likelihoods]

     if (mod[1] == 0 || mod[1] == 2) {
       lrt <- D - D.alt
     } else {
       lrt <- D.alt - D
     }

     return(lrt)
}

# ------------------------------------------------
# function to plot histograms of posterior predictive
# distributions for test quantities

  pp.plot <- function(T.obs, T.sim, name="T.sim", p.T) {

# calculate histogram of T.sim

  hist.T <- hist(T.sim, breaks=20, plot=FALSE)
  xmin <- min(hist.T$breaks, T.obs)
  xmax <- max(hist.T$breaks, T.obs)
  txt.title <- paste("Histogram of", name)
  plot(hist.T, col="blue", xlim=c(xmin, xmax),
     main=txt.title,xlab=name, bty="n", yaxt="n", ylab="")

# mark the observed valuelab

  abline(v=T.obs)

# annotate plot

  xpos <- (xmax-xmin)*0.8 + xmin
  ypos <- max(hist.T$counts)
  txt <- paste("p=", signif(p.T, 2), sep="")
  text(xpos, ypos, txt)

  }

# ------------------------------------------------
# Function to produce a matrix plot showing
# contour plots of pairs of parameters,
# similar to the output of the PAIRS function
# but using contour plots rather than scatter plots.
#
# Input:
#  theta   - N*M array of parameter values
#  n       - number of points to add to scatter plots

cont.pairs <- function(theta, n=2000, cex=1.5,
                       cont=TRUE, labels=colnames(theta)) {

# load the MASS library, which includes the KDE2D function

  require(MASS)

# number of points to plot

  np = max(n,2000)

# get the number of parameters

  M <- dim(theta)[2]

# save the existing plotting parameters,
# to be restored once finished

 par.mfrow <- par()$mfrow
 par.mfcol <- par()$mfcol
 par.mar <- par()$mar
 par.oma <- par()$oma
 par.mgp <- par()$mgp

# define a plotting array comprising M*M regions

  par(mfcol=c(M,M))

# leave some room around the edges of the plot

  par(mar=c(0,0,0,0), oma=c(6,6,6,6), mgp=c(3,3,0))

# number of data points to include on scatter plots

  n <- min(n, dim(theta)[1])

# loop over an M*M square of parameter pairs theta_i, theta_j

  for (i in 1:M) {
    for (j in 1:M) {

# calculate a two dimensional density from the
# values of {theta(,i), theta(,j)} and plot the
# density contours

      if (i != j) {
        z <- kde2d(theta[,i], theta[,j])
        contour(z, nlevels=7, drawlabels=FALSE, xlab="",
                ylab="",main="",axes=FALSE)

# draw a box around the plot region

         box(which="plot")

# add axes labels and tick marks only to outer
# plots (as in the PAIRS function)

         if (i == 1 & j %% 2 == 0) {
            axis(2,label=TRUE)
         }
         if (i == M & j %% 2 == 1) {
            axis(4,label=TRUE)
         }
         if (i %% 2 == 1 & j == M) {
            axis(1,label=TRUE)
         }
         if (i %% 2 == 0 & j == 1) {
            axis(3,label=TRUE)
         }
      }

# if we're on the diagonal do not plot any data, just
# give the name of the ith parameter

      if (j == i) {
        plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
             xlim=c(0,1), ylim=c(0,1))
        text(0.5, 0.5, labels[i], cex=cex)
      }

# if we're in the upper right half of the plot array add
# the points to make a scatter plot.

      if (i > j) {
        if (cont == TRUE) {
          contour(z, nlevels=7, drawlabels=FALSE, add=TRUE, col="grey45", lwd=1)
        }
        points(theta[1:np,i], theta[1:np,j], pch=16, cex=0.2)
      }            # end if (i > j)
    }              # end loop over j
  }                # end loop over i

# restore the graphics device parameters

  par(mfcol=par.mfcol, mfrow=par.mfrow, oma=par.oma, mar=par.mar, mgp=par.mgp)

}


# ------------------------------------------------


# ------------------------------------------------
# ------------------------------------------------
# main function
#
# N <- full length of each MCMC chain
# J <- number of chains
# Q <- number of simulations to use in ppp tests (Q <= N*J/2)
# ps <- TRUE/FALSE: output to PostScript file?
# pois <- TRUE/FALSE are the data in ct/s units subject to Poisson noise?
# prep <- TRUE/FALSE try pre-processing of time series?
# a.low <- [fixed] low-frequency slope in bending spectral models
# mod <- choice of model (integer), or if FALSE then prompt user to select model
# t.range <- time intervals (numeric vector), or if FALSE prompt user to select range
# cols <- columns to load from file, or if FALSE prompt user to choose
# par.init <- initial values for model parameters, or if FALSE prompt user
# pause <- if TRUE then pause between certain steps to view onscreen displays
# file.out <- if TRUE then save the posterior parameter summaries to a text file
#
# packages needed:
#   MASS, MNORMT
#
# ------------------------------------------------

spritz <- function(file = "",
           N = 5000,
           J = 5,
           Q = 1000,
           ps = FALSE,
           pois = FALSE,
           prep = FALSE,
           a.low = 1.0,
           mod = FALSE,
           t.range = FALSE,
           cols = FALSE,
           par.init = FALSE,
           pause = TRUE,
           file.out = FALSE) {

# define eta: the minimum (fractional) allowed spread in sampling
# (considered as numerical jitter)

  eta <- 1.0e-3

# deide to output to screen or PS file (if q.ps == TRUE)

   if (ps == TRUE) {
     postscript("bayes.ps", family="Times", horizontal=TRUE)
     par(cex=1.5, lwd=1.5, ps=15, las=1)
   }
  layout(1)

# ------------------------------------------------
# load time series data

# select the file to load

  if (file == "") {
    file <- file.choose()
  }

# get the 'base' of the filename (exclude the path) and
# extract its root (i.e. before file extention)

  if (file.out == TRUE) {
    file.base <- basename(file)
    file.bits <- unlist(strsplit(file.base, "\\."))
    file.n <- length(file.bits)
    file.root <- paste(file.bits[1:(file.n-1)], collapse=".")
    text.out <- paste(file.root, "_bayes.txt", sep="")
    ps.out <- paste(file.root, "_bayes.ps", sep="")
  } else {
    ps.out <- "bayes.ps"
  }

# load data table from file

  data <- read.table(file, skip=1)

# assumes there are two columns: TIME, VARIABLE
# if there are more than two columns, prompt the user
# to select which to use. Or take values specified upon
# calling (when cols not FALSE)

  dims <- dim(data)
  n.cols <- dims[2]

  if (cols == FALSE) {

    cols <- c(1,2)

    if (n.cols != 2) {
      cat("----------------------------------------- \n")
      cat("-- The chosen file contains ",n.cols," columns", fill=TRUE)
      cat("-- Which columns contain the TIME and VARIABLE? [1, 2]", fill=TRUE)
      str <- readline("--> ")

      if (str != "") {
        cols <-  as( unlist( strsplit(str, split=" ") ), "numeric")
      }
    }
  }

  t <- data[,cols[1]]
  x <- data[,cols[2]]

  dt <- abs(t[2]-t[1])
  n <- length(t)

# ------------------------------------------------
# check for uneven sampling

# calculate first differences

 t.diff <- t[2:n] - t[1:(n-1)]

# calculate sd/mean of first differences

 t.spread <- sd(t.diff) / mean(t.diff)

# is this unacceptable? (compared to numerical error)

  uneven <- FALSE
  if (t.spread > eta) { uneven <- TRUE }
  if (uneven == TRUE) {
    cat("----------------------------------------- \n")
    cat("-- Uneven sampling of time series", fill=TRUE)
    cat("-- Std dev. of sampling rate is ", 100*t.spread, "%", fill=TRUE)
    cat("-- Mean sampling rate is ", mean(t.diff), fill=TRUE)
    cat("-- Will resample to even grid - PLEASE CHECK RESULTS CAREFULLY", fill=TRUE)
    intp <- approx(t, x, n=n)
    t.old <- t
    x.old <- x
    t <- intp$x
    x <- intp$y
    dt <- t[2] - t[1]
  }

# ------------------------------------------------


# Poisson noise level for photon counting data
  if (pois == TRUE) {
    pn <- 2.0/mean(x, na.rm=TRUE)
  }

  cat("----------------------------------------- \n")
  cat("-- Time series: ",file,sep="",fill=TRUE)
  cat("-- Length ",length(x),", sampling interval ",dt,sep="",fill=TRUE)

# ------------------------------------------------
# plot the time series

  plot(t, x, type="s", bty="n")
  if (uneven == TRUE) {
    lines(t.old, x.old, type="b", col="green")
  }

# ------------------------------------------------
# allow user to select a subset of the data
# unless already specified at calling (t.range not FALSE)

  if (t.range == FALSE) {
    print("-- Select the range of data to use.")
    t.range <- c(min(t), max(t))
    cat("-- Current range [",t.range[1], ",", t.range[2], "] (RETURN to accept)", fill=TRUE)
    str <- readline("--> ")
    if (str != "") {
       t.range <-  as( unlist( strsplit(str, split=" ") ), "numeric")
       mask <- (t >= t.range[1] & t <= t.range[2])
       t <- t[mask]
       x <- x[mask]
       plot(t, x, type="s", bty="n")
    }
  }

# ------------------------------------------------
# pre-process data

  if (prep == TRUE) {
    print("-- Choose pre-processing of time series.")
    print("-- 0: none")
    print("-- 1: end-matching")
    print("-- 2: first differences")
    print("-- 3: subtract LOWESS smoothed trend line")
    print("-- 4: log transform")

    pre.method <- readline(prompt = "--> ")

    n.t <- length(t)

    if (pre.method == 1) {
      dy <- x[n.t] - x[1]
      dx <- t[n.t] - t[1]
      mm <- dy/dx
      cc <- x[1] - mm*t[1]
      x.mod <- mm*t + cc
      x <- x - x.mod
      plot(t, x, type="s", bty="l")
    }

    if (pre.method == 2) {
      dy <- x[2:n.t] - x[1:n.t-1]
      dx <- 0.5 * (t[2:n.t] + t[1:n.t-1])
      t <- dx; x <- dy
      plot(t, x, type="s", bty="l")
    }

    if (pre.method == 3) {
      smoo <- lowess(t, x, f=0.1)
      x <- x - smoo$y
      plot(t, x, type="s", bty="l")
    }

    if (pre.method == 4) {
       x <- log10(x)
       plot(t, x, type="s", bty="l")
    }

  }

# ------------------------------------------------
# compute and plot periodogram

  if (pause == TRUE) {
    tmp <- readline(prompt = "-- Press return to see periodogram. \n")
  }

  x.t <- x
  pergram <- pgram(x, dt, rms=TRUE)
  y <- pergram$p
  x <- pergram$f
  y.smooth <- pergram$p.smooth

  nx <- length(x)

# ------------------------------------------------
# choose model

  if (mod == FALSE) {
    print("-- Which model do you want?")
    print("-- 0:  Power law + constant         ")
    print("-- 1:  Bending power law + constant ")
    print("-- 2:  Power law                    ")
    print("-- 3:  Bending power law            ")
    mdl <- readline(prompt="--> ")
  } else {
    mdl <- mod
  }

# log(midpoint) of freq scale. used as initial guess of
# the break/bend frequency.

  f.br <- 0.5*(max(log(x)) + min(log(x)))

  if (mdl == "0") {
    mod <- c(0, a.low)
    p.noise <- log(mean(y[(nx-10):nx]))
    if (pois == TRUE) { p.noise <- pn }
    par.0 <- c(1.5, -5, p.noise)
    y.mod <- model(par.0, x, mod)
    par.0[2] <- par.0[2] + log(sum(y)/sum(y.mod))
  }

  if (mdl == "1") {
    mod <- c(1, a.low)
    p.noise <- log(mean(y[(nx-10):nx]))
    if (pois == TRUE) { p.noise <- pn }
    par.0 <- c(2.8, f.br, -5, p.noise)
    y.mod <- model(par.0, x, mod)
    par.0[3] <- par.0[3] + log(sum(y)/sum(y.mod))
  }

  if (mdl =="2") {
    mod <- c(2, a.low)
    par.0 <- c(1.5, -5)
    y.mod <- model(par.0, x, mod)
    par.0[2] <- par.0[2] + log(sum(y)/sum(y.mod))
  }

  if (mdl == "3") {
    mod <- c(3, a.low)
    par.0 <- c(1.3, f.br, -5)
    y.mod <- model(par.0, x, mod)
    par.0[3] <- par.0[3] + log(sum(y)/sum(y.mod))
  }

# use the user-supplied starting parameters if given (by par.init)
# If par.init == FALSE then prompt the user to update initial guess
# If par.init == TRUE then use initial guess without prompting user

  if (par.init != FALSE && par.init != TRUE) {
     par.0 <- par.init
  }

# calculate model using initial parameter values

  y.mod <- model(par.0, x, mod)
  lines(x, y.mod, type="l", col="purple")

  print("-- Initial parameter values")
  cat("-- ", par.0, fill=TRUE)

  if (par.init == FALSE) {
    print("-- Enter alternative values(or RETURN to accept)")
    par.str <- readline(prompt = "-->")
    while (par.str != "") {
      par.0 <- as( unlist( strsplit(par.str, split=" ") ), "numeric")
      y.mod <- model(par.0, x, mod)
      lines(x, y.mod, type="l", col="purple")
      print("-- Enter alternative values(or RETURN to accept)")
      par.str <- readline(prompt = "-->")
    }
  }

# ------------------------------------------------
# define model spectrum and fit to periodogram data

  fit.obs <- fit.model(x, y, mod=mod, par.0, obs=TRUE, ylab=pergram$ylab)

  theta.map <- fit.obs$theta.map
  theta.cov <- fit.obs$theta.cov
  f.peak=fit.obs$f.peak

# ------------------------------------------------
# calculate LRT statistic

  lrt.obs <- calc.lrt(x, y, mod=mod, theta.map, fit.obs$D)
  cat("-- Likelihood Ratio Test (LRT) = ", lrt.obs, fill=TRUE,sep="")
  p.lrt <- pchisq(lrt.obs, df=1, lower.tail=FALSE)
  cat("-- LRT test p-value = ", p.lrt, fill=TRUE, sep="")

# ------------------------------------------------
# generate the MCMC chains

# set up the covariance matrix for the
# (Normal) proposal distribution used to
# generate the MCMC outp

  cov <- theta.cov * 1.2

# N = (full) length of each chain
# J = number of chains
# call the mcmc routines

  mcmc.out <- mcmc(x, y, N, J, mod=mod, theta.map, cov)
  L <- dim(mcmc.out)[2]

# get parameter inferences from MCMC draws

  mcmc.infer <- infer(mcmc.out)

# save to a text file if requested

  if (file.out == TRUE) {
	write.table(mcmc.infer, file=text.out)
  }

# posterior predictive model checks

  pp.out <- pp(mcmc.out, x, mod=mod, Q, f.peak)

  T.sim <- pp.out$T.sim
  D.sim <- pp.out$D.sim
  ks.sim <- pp.out$ks.sim
  sse.sim <- pp.out$sse.sim
  srat.sim <- pp.out$srat.sim
  y0.sim <- pp.out$y0

  lrt.sim <- pp.out$lrt.sim

# ------------------------------------------------

# plot histograms of the posterior predictive distributions
# of statistics

  T.obs <- fit.obs$T
  D.obs <- fit.obs$D
  ks.obs <- fit.obs$ks
  sse.obs <- fit.obs$sse
  srat.obs <- fit.obs$srat

# calculate p-values of statistics

  p.T <- mean(T.sim > T.obs)
  p.D <- mean(D.sim > D.obs)
  p.ks <- mean(ks.sim > ks.obs)
  p.sse <- mean(sse.sim > sse.obs)
  p.lrt <- mean(lrt.sim > lrt.obs)
  p.srat <- mean(srat.sim > srat.obs)

# and 'errors' on p-values
# (from binomial formula sd = sqrt[p*(1-p)N] )

  p.T.err <- sqrt( p.T * (1-p.T) / length(T.sim) )
  p.D.err <- sqrt( p.D * (1-p.D) / length(D.sim) )
  p.ks.err <- sqrt( p.ks * (1-p.ks) / length(ks.sim) )
  p.sse.err <- sqrt( p.sse * (1-p.sse) / length(sse.sim) )
  p.lrt.err <- sqrt( p.lrt * (1-p.lrt) / length(lrt.sim) )
  p.srat.err <- sqrt( p.srat * (1-p.srat) / length(srat.sim) )

# display the results

  cat("-- Bayesian p-value for T   : p = ", p.T,  " [+/-",signif(p.T.err, digits=3),"]", fill=TRUE, sep="")
  cat("-- Bayesian p-value for D   : p = ", p.D,  " [+/-",signif(p.D.err, digits=3),"]", fill=TRUE, sep="")
  cat("-- Bayesian p-value for KS  : p = ", p.ks, " [+/-",signif(p.ks.err, digits=3),"]", fill=TRUE, sep="")
  cat("-- Bayesian p-value for SSE : p = ", p.sse," [+/-",signif(p.sse.err, digits=3),"]", fill=TRUE, sep="")
  cat("-- Bayesian p-value for Sm.T: p = ", p.srat," [+/-",signif(p.srat.err, digits=3),"]", fill=TRUE, sep="")
  cat("-- Bayesian p-value for LRT : p = ", p.lrt," [+/-",signif(p.lrt.err, digits=3),"]", fill=TRUE, sep="")


# plot histograms of T and D

  png("bayes2.png", width = 960, height = 480)


  layout(rbind(c(1,2,3),c(4,5,6)))

  pp.plot(T.obs, T.sim, name="T.sim", p.T)
  pp.plot(srat.obs, srat.sim, name="srat.sim", p.srat)
  pp.plot(lrt.obs, lrt.sim, name="lrt.sim", p.lrt)
  pp.plot(sse.obs, sse.sim, name="sse.sim", p.sse)
  pp.plot(D.obs, D.sim, name="D.sim", p.D)
  pp.plot(ks.obs, ks.sim, name="ks.sim", p.ks)

  dev.off(which = dev.cur())

  if (ps == TRUE) {
     dev.off(dev.cur())
  }

# -------------------------------------------
# generate 'pretty' output graphical summary


#  postscript(ps.out, family="Times", horizontal=TRUE)
  png("bayes.png", width = 960, height = 960)
  par(cex=1.5, lwd=1.5, ps=15, las=1)

   layout(rbind(c(1,1,1),
               c(2,2,4),
               c(3,3,5)))

# frame 1: time series

  plot(t, x.t, type="l", bty="l")

# frame 2: periodogram + model

  dx <- x[1]
  x.plot <- c(x[1]-dx/2,x)
  y.plot <- model(theta.map, x.plot, mod=mod)
  plot(x-dx/2, y, log="xy", type="s", xlab="Frequency", ylab="Power", bty="l")

  smooth.per <- lowess(log(x), log(y), f=0.2)
  smooth.x <- exp(smooth.per$x)
  smooth.y <- exp(smooth.per$y + 0.557214)
  lines(smooth.x, smooth.y, type="l", col="blue", lw=5)

  lines(x.plot, y.plot, lty=1, lwd=4, col="red")

# frame 3: residuals

  y.map <- model(theta.map, x, mod=mod)
  rat <- y/y.map
  plot(x-dx/2, rat, type="s", bty="l", log="xy", xlab="Frequency",
        ylab="data/model", ylim=c(0.1, max(10,max(rat)*1.2)))

  filt <- c(0.25, 0.5, 0.25)
  filter.y <- filter(y, filter=filt)
  lines(x-dx/2, filter.y/y.map, type="s", col="blue", lw=3)

  lines(range(x)-dx/2, c(1,1), lwd=4, col="red")

# frame 4: histogram of T

  pp.plot(srat.obs, srat.sim, name="srat.sim", p.srat)

# frame 5: histogram of SSE

  pp.plot(sse.obs, sse.sim, name="sse.sim", p.sse)

  dev.off(which = dev.cur())

# -------------------------------------------

# return the periodogram data to the user

  result <- list(f=x, p=y, p.smooth=y.smooth)
  return(result)

# end of BAYES routine

  }

# ------------------------------------------------

