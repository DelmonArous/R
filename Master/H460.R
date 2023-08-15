rm(list=ls())

## Define the LQ and USC functions
LQfunc = function(alpha, beta, D) -(alpha)*D - (beta)*I(D^2)
f = function(D, DT) I(D^2) - (D > DT)*I((D - DT)^2)
USCfunc = function(alpha, beta, SF, D, DT) {
  
  fit = (-(alpha)*D - (beta)*I( f(D, DT) ))
  residuals = SF - fit
  RSS = c(crossprod(residuals))
  
  list(RSS = RSS)
}
  
## Define function for computing RSS, MSS, TSS, R-squared and adjusted R-squared
Sxx = function(x)   sum( I((x - mean(x))^2) )
s = function(y, yhat, n, p) sqrt( RSS(y, yhat) / (n - p) )
sy = function(x, y, yhat, n, p) s(y, yhat, n, p)*sqrt( (1/n) + (x - mean(x))^2/Sxx(x) )
RSS = function(y ,yhat) sum( I((y - yhat)^2) ) 
MSS = function(y, yhat) sum( (yhat - mean(y))^2 )
TSS = function(y, yhat) sum( (y - mean(y))^2 )
R.squared = function(y, yhat) 1 - (RSS(y, yhat)/TSS(y, yhat))
R.squared.adj = function(y, yhat, n, p) 1 - (1 - R.squared(y, yhat))*( (n - 1)/(n - p) ) 
standardresid = function(x, y, yhat, n, p) (y - yhat)/(s(y, yhat, n, p)*sqrt(1 - (1/n) - (x - mean(x))^2/Sxx(x)))

## Computes the regression when given the data points (x,y), and DT is the transition dose DT
computeFit = function (D, SF, DT = NULL) {
  
  ## Model matrix
  if ( is.null(DT) == TRUE ) {
    X = cbind( "alpha" = D, "beta" = I(D^2) )
  } else {
    X = cbind( "alpha" = D, "beta" = f(D,DT) )
  }
  
  ## Number of predictor variables, i.e. the dimension of the model matrix
  p = dim(X)[2L]
  
  ## Use QR factorization to solve least squares
  fit = .lm.fit(X, SF)
  
  ## Compute the residuals, RSS and estimate the standard deviation sigma^2
  residuals = c(fit$residuals)
  n = length(residuals)
  RSS = c(crossprod(residuals))
  sig2 = RSS/(n - p)
  
  ## Compute the regression coefficient summary table
  beta = fit$coef
  R = "dimnames<-"(fit$qr[1:p, ], NULL)
  Rinv = backsolve(R, diag(p))
  se = sqrt(rowSums(Rinv ^ 2) * sig2)
  tstat = beta/se
  pval = 2*pt(abs(tstat), n - p, lower.tail = FALSE)

  ## Compute AIC
  AIC = n*log(2*pi*sig2) + (n - p) + 2*(p + 1) 
  
  ## Compute R-squared and adjusted R-squared
  TSS = c(crossprod(SF - sum(SF) / n))
  r.squared = 1 - RSS/TSS
  adj.r.squared = 1 - sig2*(n - 1)/TSS
  
  ## Return estimated values
  if ( is.null(DT) == TRUE ) {
    list(coef = beta, residuals = residuals, fitted.values = c(X %*% beta),
         R = R, sig2 = sig2, se = se, AIC = AIC,
         RSS = RSS, r.squared = r.squared, adj.r.squared = adj.r.squared, p = p, n = n)
  } else {
    list(coef = beta, residuals = residuals, fitted.values = c(X %*% beta),
         R = R, sig2 = sig2, se = se, AIC = AIC, DT = DT,
         RSS = RSS, r.squared = r.squared, adj.r.squared = adj.r.squared, p = p, n = n)
  }
  
}

## Different vectors containing interval values to aid the following computation  
Doserange = 17
DTrange = c(3, 5)
legend.title = c("Asynchronous H460 cells")
D0.vec = c(0.86)
Dq.vec = c(1.52)

# Extracted data set
D = c(0, 2.105, 5, 7.5, 10, 11.97, 13.95, 15.92)
SFdata = exp(c(0, -0.2, -0.8, -1.9, -3.4, -4.5, -5.65, -6.8))

## Linear and non-linear regression of LQ and MTSH, respectively
LQmod = lm(log(SFdata) ~ -1 + D + I(D^2))
MTSHmod = nls(log(SFdata) ~ log(1 - (1 - exp(-v*D))^m), start = list(v=0.5, m=1))

## Model coefficients
alpha = (-1)*coef(LQmod)[1]
beta = (-1)*coef(LQmod)[2]
v = coef(MTSHmod)[1] 
m = coef(MTSHmod)[2] 
D0 = 1/v
Dq = D0*log(m)
DT = (2*Dq)/(1 - alpha*D0)

## Compute USC curve 
USC = (-(alpha)*D - (beta)*I( f(D, DT) ))

## Determine the transition dose as an independent parameter and find best USC fit

DTspan = seq(5, 11, 0.01)
SF = log(SFdata)

## Regression for each value of DT in span
# lst = lapply(DTspan, computeFit, D = D, SF = SF)
lst = lapply(DTspan, USCfunc, alpha = alpha, beta = beta, SF = SF, D = D)

## Compute RSS for each value of DT
RSS_USC = sapply(lst, "[[", "RSS")

## Plot RSS vs grid of DT
dev.new()
plot(DTspan, RSS_USC,
     ylab = "RSS", xlab = "DT grid [Gy]",
     ylim = c(min(RSS_USC), max(RSS_USC)), xlim = c(min(DTspan), max(DTspan)),
     type = "o", pch = c(16, rep(NA, 24)),
     cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
legend("topright", legend=c("H460 NSCLC cell line"), pch = 16, lty = 1, bty = "o")

# ## Plot dose-response curves for the cell population under aerobic and extremely-hypoxic cells
dev.new()
plot(D, log(SFdata),
     ylab = "log SF",
     las = 1, pch = 16, xaxt = "n", ann = FALSE, 
     cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
axis(3, cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
mtext("Dose [Gy]", side = 3, line = 3, cex = 1.5)
mtext("log SF", side = 2, line = 3, cex = 1.5)
D_range = seq(0, max(D), 0.01)
lines(D_range, (LQfunc(alpha, beta, D_range)), lty = 2)
lines(D_range, ((-(alpha)*D_range - (beta)*I( f(D_range, DT) ))), lty = 1)

#lines(sort(D), (exp(fitted(LQmod))[order(D)]), lty = 2)
#lines(sort(D), (exp(fitted(MTSHmod))[order(D)]), lty = 2)
#lines(sort(D), (exp(USC)[order(D)]), lty = 1)
legend("bottomleft", 
       legend=c("H460 NSCLC cell line", 
                paste0("USC fit, adjusted R-squared = ", signif(R.squared.adj(log(SFdata), USC, length(D), 4), digits = 3)),
                paste0("LQ fit, adjusted R-squared = ", signif(summary(LQmod)$adj.r.squared, digits = 3))),
                lty = c(0, 1, 2), bty = "n", cex = 1.0)

cat("Adjusted R-sqaured: ", summary(LQmod)$adj.r.squared, "\n")
cat("Adjusted R-sqaured: ", R.squared.adj(log(SFdata), USC, length(D), 4), "\n")

fit = computeFit(D, SF, DTspan[which(RSS_USC == min(RSS_USC))])
USC_opt = USCfunc(-0.01048588, 0.03596738, D, DTspan[which(RSS_USC == min(RSS_USC))])
cat("Adjusted R-sqaured: ", R.squared.adj(log(SFdata), USC_opt, length(D), 3), "\n")

SDalpha = coef(summary(LQmod))[1, 2]
SDbeta = coef(summary(LQmod))[2, 2]
SDv = coef(summary(MTSHmod))[1, 2]
SDm = coef(summary(MTSHmod))[2, 2]
SDD0 = D0^2*SDv
SDDq = sqrt( (SDD0*log(m))^2 + (SDm*D0/m)^2 )
SDDT = (DT/Dq)*sqrt( SDDq^2 + (D0*DT*SDalpha/2)^2 + (alpha*DT*SDD0/2)^2 )

