rm(list=ls())

f = function(x, c) I(x^2) - (x > c)*I((x - c)^2)

getX = function (x, c) {
  cbind( "beta0" = x, "beta1" = f(x,c) )
}

#getX <- function (x, c) {
  x <- x - c
  cbind("beta0" = 1, "beta1" = x, "beta2" = pmin(x, 0) ^ 2)
}

## `x`, `y` give data points; `c` is known break point
est <- function (x, y, c) {
  ## model matrix
  X <- getX(x, c)
  p <- dim(X)[2L]
  ## solve least squares with QR factorization
  fit <- .lm.fit(X, y)
  ## compute Pearson estimate of `sigma ^ 2`
  r <- c(fit$residuals)
  n <- length(r)
  RSS <- c(crossprod(r))
  sig2 <- RSS / (n - p)
  ## coefficients summary table
  beta <- fit$coefficients
  R <- "dimnames<-"(fit$qr[1:p, ], NULL)
  Rinv <- backsolve(R, diag(p))
  se <- sqrt(rowSums(Rinv ^ 2) * sig2)
  tstat <- beta / se
  pval <- 2 * pt(abs(tstat), n - p, lower.tail = FALSE)
  tab <- matrix(c(beta, se, tstat, pval), nrow = p, ncol = 4L,
                dimnames = list(dimnames(X)[[2L]], 
                                c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
  ## 2 * negative log-likelihood
  nega2logLik <- n * log(2 * pi * sig2) + (n - p)
  ## AIC / BIC
  aic <- nega2logLik + 2 * (p + 1)
  bic <- nega2logLik + log(n) * (p + 1)
  ## multiple R-squared and adjusted R-squared
  TSS <- c(crossprod(y - sum(y) / n))
  r.squared <- 1 - RSS / TSS
  adj.r.squared <- 1 - sig2 * (n - 1) / TSS
  ## return
  list(coefficients = beta, residuals = r, fitted.values = c(X %*% beta),
       R = R, sig2 = sig2, coef.table = tab, aic = aic, bic = bic, c = c,
       RSS = RSS, r.squared = r.squared, adj.r.squared = adj.r.squared, p = p, n = n)
}

choose.c <- function (x, y, c.grid) {
  if (is.unsorted(c.grid)) stop("'c.grid' in not increasing")
  ## model list
  lst <- lapply(c.grid, est, x = x, y = y)
  ## RSS trace
  RSS <- sapply(lst, "[[", "RSS")
  ## verbose
  dev.new()
  plot(c.grid, RSS, type = "b", pch = 19)
  ## find `c` / the model minimizing `RSS`
  lst[[which.min(RSS)]]
}

pred <- function (model, x.new, y) {
  ## prediction matrix
  X <- getX(x.new, model$c)
  p <- dim(X)[2L]
  ## predicted mean
  fit <- X %*% model$coefficients
  ## prediction standard error
  Qt <- forwardsolve(t(model$R), t(X))
  se <- sqrt(colSums(Qt ^ 2) * model$sig2)
  ## 95%-confidence interval
  alpha <- qt(0.025, length(model$residuals) - p)
  lwr <- fit + alpha * se
  upr <- fit - alpha * se
  ## 95%-prediction interval
  tval = qt((1-0.95)/2, df = model$n - model$p)
  #lwr = fit + tval*s(y,fit,model$n,model$p)*sqrt( 1/(model$n-model$p) + I((x.new - mean(x.new))^2)/Sxx(x.new) )
  #upr = fit - tval*s(y,fit,model$n,model$p)*sqrt( 1/(model$n-model$p) + I((x.new - mean(x.new))^2)/Sxx(x.new) )
  ## return
  matrix(c(fit, se, lwr, upr), ncol = 4L,
         dimnames = list(NULL, c("fit", "se", "lwr", "upr")))
}

Sxx = function(x) sum( I((x - mean(x))^2) )
s = function(y, yhat, n, p) sqrt( RSS(y, yhat) / (n - p) )
RSS = function(y ,yhat) sum( I((y - yhat)^2) ) 
ESS = function(y, yhat) sum( (yhat - mean(y))^2 ) # SSR or MSS
TSS = function(y, yhat) sum( (y - mean(y))^2 )
R.squared = function(y, yhat) 1 - (RSS(y, yhat)/TSS(y, yhat))
R.squared.adj = function(y, yhat, n, p) 1 - (1 - R.squared(y, yhat))*( (n - 1)/(n - p) ) 


## we first generate a true model
set.seed(17)
x <- runif(200, min = 0, max = 15)  ## sample points on [0, 1]
beta <- -c(0.3389, 0.0783)  ## true coefficients
X <- getX(x, 8)  ## model matrix with true break point at 0.6
y <- X %*% beta + rnorm(200, 0, 1)  ## observations with Gaussian noise
dev.new()
plot(x, exp(y), log = "y")

c.grid <- seq(1, 14, 0.5)
fit <- choose.c(x, y, c.grid)
fit$c

#x.new <- sort(x)
x.new = seq(0, 15, 0.5)
p <- pred(fit, x.new, y)

dev.new()
plot(x, exp(y), cex = 0.5, log = "y")
matlines(x.new, exp(p[,-2]), col = c(1,2,2), lty = c(1,2,2), lwd = 2)


RSS(y, fit$fitted.values)
fit$RSS

fit$r.squared
R.squared(y, fit$fitted.values )

fit$adj.r.squared
R.squared.adj(y, fit$fitted.values, fit$n, fit$p)

fit$sig2
s(y,fit$fitted.values,fit$n,fit$p)^2

tval.mine = qt((1-0.95)/2, df = fit$n - fit$p)
lwr.mine = fit$fitted.values + tval.mine*s(y,fit$fitted.values,fit$n,fit$p)*sqrt( 1/(fit$n) + (x.new - mean(x.new))^2 / Sxx(x.new) )
upr.mine = fit$fitted.values - tval.mine*s(y,fit$fitted.values,fit$n,fit$p)*sqrt( 1/(fit$n) + (x.new - mean(x.new))^2 / Sxx(x.new) )
lwr = p[,3]
upr = p[,4]
