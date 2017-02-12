
t <- seq(0, 100, by = 0.07)
n <- length(t)
y <- cumsum(rnorm(n))*0.3 + rnorm(n)
y[200] <- 20
y[400] <- -40


result <- clipper(t, y, width = 2, n.min = 5, n.sigmas = 4)
plot(t, y, bty = "n", pch = 16)
lines(t, result$y.mean, col = "blue", lwd = 3)
mask <- which(result$qual == FALSE)
points(t[mask], y[mask], col="pink", cex=2, lwd = 2)
lines(t, result$y.mean + result$y.sd, lty = 2)
lines(t, result$y.mean - result$y.sd, lty = 2)

