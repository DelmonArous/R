setwd("C:/Users/Delmon/Dropbox/Master/R")
getwd()

## Clear and detach all data
rm(list=ls())
detach(SFdata)

## Read in cell survival observations and attach the object SFdata to R search path
SFdata = read.table(file = "NHIK3025_v1.csv", header = TRUE, sep = ";")
attach(SFdata)

## Define the LQ and USC functions
LQfunc = function(alpha, beta, D) -(alpha)*D - (beta)*I(D^2)
f = function(D, DT) I(D^2) - (D > DT)*I((D - DT)^2)
USCfunc = function(alpha, beta, D, DT) -(alpha)*D - (beta)*I( f(D, DT) )

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
  AIC = n*log(2*pi*sig2) + (n - p) + 2*p
  
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

## Compute the 95% CI for the fitted y values
CI = function (model, D.new) {
  
  ## Prediction matrix
  if ( is.null(model$DT) == TRUE ) {
    X = cbind( "alpha" = D.new, "beta" = D.new^2 )
    p = dim(X)[2L]
  } else {
    X = cbind( "alpha" = D.new, "beta" = f(D.new,model$DT) )
    p = 4
  }
  
  ## Fitted values or expected mean
  fit = X %*% model$coef
  
  ## Prediction standard error
  Qt = forwardsolve(t(model$R), t(X))
  se = sqrt(colSums(Qt^2)*model$sig2)
  
  ## 95% CI
  alpha = qt(0.025, length(model$residuals) - p)
  lwr = fit + alpha * se
  upr = fit - alpha * se

  ## Return estimated values
  matrix(c(fit, se, lwr, upr), ncol = 4L,
         dimnames = list(NULL, c("fit", "se", "lwr", "upr")))
  
}

## Different vectors containing interval values to aid the following computation  
stdave.O2 = c(0.07092, 0.17892, 0.2731, 0.08199)
stdave.N2 = c(0.07504, 0.06343, 0.07921, 0.03519)
Doserange.O2 = c( 15, 25, 25, 25)
Doserange.N2 = c( 30, 40, 40, 40)
DTrange.O2 = c( c(1.5, 3), c(3.5, 6), c(3, 5), c(3, 5))
DTrange.N2 = c( c(5, 8), c(13, 17), c(11, 15), c(15, 20))
legend.title = c("Mitosis (0 hours)", "G1 (5 hours)", "S (13 hours)", "Asynchronous cells")
v.O2.vec = c(1, 0.5, 0.5, 0.5)
m.O2.vec = c(1, 2, 1, 1)
v.N2.vec = c(1, 0.2, 0.2, 0.2)
m.N2.vec = c(1, 2, 2, 2)
D0.O2.vec = c(0.84, 1.08, 0.82, 0.86)
D0.N2.vec = c(0.5, 3.46, 2.55, 2.36)
Dq.O2.vec = c(0.29, 2.30, 1.66, 1.52)
Dq.N2.vec = c(3, 6.16, 5.90, 8.57)

## Loop through the different cell populations 
k = 0
for (i in c("M", "G1", "S", "Asynch") ) {
  
  k = k + 1
  
  ## Extract dose-response observations
  SF.O2 = SFdata[Cellcycle == i & Hypoxia == "No", 1]
  SF.N2 = SFdata[Cellcycle == i & Hypoxia == "Yes", 1]
  D.O2 = SFdata[Cellcycle == i & Hypoxia == "No", 2]
  D.N2 = SFdata[Cellcycle == i & Hypoxia == "Yes", 2]
  
  ## Linear and non-linear regression of LQ and MTSH, respectively
  LQmod.O2 = lm(log(SF.O2) ~ -1 + D.O2 + I(D.O2^2))
  LQmod.N2 = lm(log(SF.N2) ~ -1 + D.N2 + I(D.N2^2))
  MTSHmod.O2 = nls(log(SF.O2) ~ log(1 - (1 - exp(-v*D.O2))^m), start = list(v=v.O2.vec[k], m=m.O2.vec[k]))
  MTSHmod.N2 = nls(log(SF.N2) ~ log(1 - (1 - exp(-v*D.N2))^m), start = list(v=v.N2.vec[k], m=m.N2.vec[k]))
  
  ## Model coefficients
  alpha.N2 = (-1)*coef(LQmod.N2)[1]
  beta.N2 = (-1)*coef(LQmod.N2)[2]
  v.N2 = coef(MTSHmod.N2)[1] 
  m.N2 = coef(MTSHmod.N2)[2] 
  D0.N2 = 1/v.N2
  Dq.N2 = D0.N2*log(m.N2)
  DT.N2 = (2*Dq.N2)/(1 - alpha.N2*D0.N2)
  
  alpha.O2 = (-1)*coef(LQmod.O2)[1]
  beta.O2 = (-1)*coef(LQmod.O2)[2]
  v.O2 = coef(MTSHmod.O2)[1] 
  m.O2 = coef(MTSHmod.O2)[2] 
  D0.O2 = 1/v.O2
  Dq.O2 = D0.O2*log(m.O2)
  DT.O2 = (2*Dq.O2)/(1 - alpha.O2*D0.O2)
  
  ## Compute USC curve 
  USC.O2 = USCfunc(alpha.O2, beta.O2, D.O2, DT.O2)
  USC.N2 = USCfunc(alpha.N2, beta.N2, D.N2, DT.N2)
  
  # ## Plot dose-response curves for the cell population under aerobic and extremely-hypoxic cells
  dev.new()
  plot(D.O2, SF.O2,
       ylab = "Surviving fraction", xlab = "Dose [Gy]",
       ylim = c(min(SF.O2), 1), xlim = c(0, max(D.N2)),
       las = 1, pch = 1, yaxt = "n", log = "y",
       cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
  points(D.N2, SF.N2, pch = 16)
  ticks = seq(-4, 0, by=1)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(2, at=c(0.0001, 0.001, 0.01, 0.1, 1), labels=labels,
       cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)

  D_aer = seq(0, max(D.O2), 0.01)
  D_hyp = seq(0, max(D.N2), 0.01)
  lines(D_aer, exp(LQfunc(alpha.O2, beta.O2, D_aer)), lty = 1)
  lines(D_hyp, exp(LQfunc(alpha.N2, beta.N2, D_hyp)), lty = 1)
  lines(D_aer, (1 - (1 - exp(-v.O2*D_aer))^m.O2), lty = 2)
  lines(D_hyp, (1 - (1 - exp(-v.N2*D_hyp))^m.N2), lty = 2)
  lines(D_aer, exp(USCfunc(alpha.O2, beta.O2, D_aer, DT.O2)), lty = 3)
  lines(D_hyp, exp(USCfunc(alpha.N2, beta.N2, D_hyp, DT.N2)), lty = 3)
  legend("topright", legend=c(legend.title[k], "Aerobic", "Extremely hypoxic", "LQ", "MTSH", "USC"),
         lty = c(0, 0, 0, 1, 2, 3), pch = c(NA, 1, 16, NA, NA, NA), bty = "o")

    # AIC_MTSH_O2 = 2*2 - 2*logLik(MTSHmod.O2) # length(D.O2)*log(2*pi*sig2_MTSH_O2) + (length(D.O2) - 2) + 2*2 
    # cat("MTSH O2, AIC = ", AIC_MTSH_O2, "\n")
    # 
    # AIC_MTSH_N2 = 2*2 - 2*logLik(MTSHmod.N2) # length(D.N2)*log(2*pi*sig2_MTSH_N2) + (length(D.N2) - 2) + 2*2 
    # cat("MTSH N2, AIC = ", AIC_MTSH_N2, "\n")
    # 
    # AIC_LQ_O2 = 2*2 - 2*logLik(LQmod.O2) # length(D.O2)*log(2*pi*sig2_LQ_O2) + (length(D.O2) - 2) + 2*2 
    # cat("LQ O2, AIC = ", AIC_LQ_O2, "\n")
    # 
    # AIC_LQ_N2 = 2*2 - 2*logLik(LQmod.N2) # length(D.N2)*log(2*pi*sig2_LQ_N2) + (length(D.N2) - 2) + 2*2 
    # cat("LQ N2, AIC = ", AIC_LQ_N2, "\n")
    # 
    # sig2_USC_O2 = RSS(log(SF.O2), USC.O2)/(length(D.O2) - 4)
    # AIC_USC_O2 = length(D.O2)*log(2*pi*sig2_USC_O2) + (length(D.O2) - 4) + 2*4 
    # cat("USC O2, AIC = ", AIC_USC_O2, "\n")
    # 
    # sig2_USC_N2 = RSS(log(SF.N2), USC.N2)/(length(D.N2) - 4)
    # AIC_USC_N2 = length(D.N2)*log(2*pi*sig2_USC_N2) + (length(D.N2) - 4) + 2*4 
    # cat("USC N2, AIC = ", AIC_USC_N2, "\n")
  
  ## Compute model residuals
  LQresid.O2 = resid(LQmod.O2)
  LQresid.N2 = resid(LQmod.N2)
  MTSHresid.O2 = resid(MTSHmod.O2) 
  MTSHresid.N2 = resid(MTSHmod.N2) 
  USCresid.O2 = log(SF.O2) - USC.O2 
  USCresid.N2 = log(SF.N2) - USC.N2 

  ## Aerobic residual plot
  dev.new()
  plot(exp(fitted(LQmod.O2)), LQresid.O2,
       ylab = "Residuals", xlab = "Estimated surviving fraction (fitted values)",
       ylim = c(-2, 2), xlim = c(min(SF.O2), 1),
       log = "x", las = 1, pch = 1, xaxt = "n",
       cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
  points(exp(fitted(MTSHmod.O2)), MTSHresid.O2, pch = 2)
  points(exp(USC.O2), USCresid.O2, pch = 3)
  abline(0, 0, lty = 2)
  ticks = seq(-4, 0, by=1)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=c(0.0001, 0.001, 0.01, 0.1, 1), labels=labels,
       cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)

  ## Add trend curves in the aerobic residual plot for all three models
  LOESS.MTSHresid.O2 = loess(MTSHresid.O2 ~ fitted(MTSHmod.O2) )
  LOESS.LQresid.O2 = loess(LQresid.O2 ~ (fitted(LQmod.O2)) )
  LOESS.USCresid.O2 = loess(USCresid.O2 ~ USC.O2)

  lines(sort(exp(fitted(LQmod.O2))), fitted(LOESS.LQresid.O2)[order(exp(fitted(LQmod.O2)))], lty = 1)
  lines(sort(exp(fitted(MTSHmod.O2))), fitted(LOESS.MTSHresid.O2)[order(exp(fitted(MTSHmod.O2)))], lty = 2)
  lines(sort(exp(USC.O2)), fitted(LOESS.USCresid.O2)[order(exp(USC.O2))], lty = 3)
  legend("topright", legend = c(paste(legend.title[k], "- Aerobic"), "LQ", "MTSH", "USC"),
         lty = c(NA, 1, 2, 3), pch = c(NA, 1, 2, 3), bty = "o")

  ## Extremely-hypoxic residual plot
  dev.new()
  plot(exp(fitted(LQmod.N2)), LQresid.N2,
       ylab = "Residuals", xlab = "Estimated surviving fraction (fitted values)",
       ylim = c(-2, 2), xlim = c(min(SF.N2), 1),
       log = "x", las = 1, pch = 1, xaxt = "n",
       cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
  points(exp(fitted(MTSHmod.N2)), MTSHresid.N2, pch = 2)
  points(exp(USC.N2), USCresid.N2, pch = 3)
  abline(0, 0, lty = 2)
  ticks = seq(-4, 0, by=1)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=c(0.0001, 0.001, 0.01, 0.1, 1), labels=labels,
       cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)

  ## Add trend curves in the extremely-hypoxic residual plot for all models
  LOESS.LQresid.N2 = loess(LQresid.N2 ~ (fitted(LQmod.N2)) )
  LOESS.MTSHresid.N2 = loess(MTSHresid.N2 ~ fitted(MTSHmod.N2) )
  LOESS.USCresid.N2 = loess(USCresid.N2 ~ USC.N2)

  lines(sort(exp(fitted(LQmod.N2))), fitted(LOESS.LQresid.N2)[order(exp(fitted(LQmod.N2)))], lty = 1)
  lines(sort(exp(fitted(MTSHmod.N2))), fitted(LOESS.MTSHresid.N2)[order(exp(fitted(MTSHmod.N2)))], lty = 2)
  lines(sort(exp(USC.N2)), fitted(LOESS.USCresid.N2)[order(exp(USC.N2))], lty = 3)
  legend("topright", legend = c(paste(legend.title[k], "- Extremely hypoxic"), "LQ", "MTSH", "USC"),
         lty = c(NA, 1, 2, 3), pch = c(NA, 1, 2, 3), bty = "o")

  ## Determine the transition dose as an independent parameter and find best USC fit

  DTspan.O2 = seq(DTrange.O2[2*k-1], DTrange.O2[2*k], 0.01)
  D = D.O2
  SF = log(SF.O2)

  ## Regression for each value of DT in span
  lst.O2 = lapply(DTspan.O2, computeFit, D = D.O2, SF = SF)

  ## Compute RSS for each value of DT
  RSS.O2 = sapply(lst.O2, "[[", "RSS")

  DTspan.N2 = seq(DTrange.N2[2*k-1], DTrange.N2[2*k], 0.01)
  D = D.N2
  SF = log(SF.N2)

  ## Regression for each value of DT in span
  lst.N2 = lapply(DTspan.N2, computeFit, D = D, SF = SF)

  ## Compute RSS for each value of DT
  RSS.N2 = sapply(lst.N2, "[[", "RSS")

  ## Plot RSS vs grid of DT
  dev.new()
  
  ## Add extra space to right margin of plot within same frame
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  ## Plot first set of data and draw its axis
  plot(DTspan.O2, RSS.O2,
       ylab = "", xlab = "",
       ylim = c(min(RSS.O2), max(RSS.O2)),
       xlim = c(min(DTspan.O2), max(DTspan.O2)),
       axes = FALSE, type = "o", col = "blue", pch = c(1, rep(NA, 24)),
       cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
  axis(2, ylim = c(min(RSS.O2), max(RSS.O2)), 
       col = "blue", col.axis = "blue")
  mtext("RSS", side = 2, col = "blue", line = 2.5)
  box()
  
  ## Plot the second plot and put axis scale on the right
  par(new = TRUE)
  plot(DTspan.N2, RSS.N2, 
       ylab = "", xlab = "",
       ylim = c(min(RSS.N2), max(RSS.N2)),
       xlim = c(min(DTspan.N2), max(DTspan.N2)),
       axes = FALSE, type = "o", col = "red", pch = c(16, rep(NA, 24)))
  axis(4, ylim = c(min(RSS.N2), max(RSS.N2)), 
       col = "red", col.axis = "red")
  mtext("RSS", side = 4, col = "red", line = 2.5)
  
  ## Draw DT axis
  axis(1, pretty(c(min(min(DTspan.O2), min(DTspan.N2)), 
                   max(max(DTspan.O2), max(DTspan.N2))), 2))
  mtext("DT grid [Gy]", side = 1, col = "black", line = 2.5)
  
  ## Add legend
  legend("topright", 
         legend=c(legend.title[k], "Aerobic", "Extremely hypoxic"),
         col = c(NA, "blue", "red"), pch = c(NA, 1, 16), 
         lty = c(NA, 1, 1), bty = "o")

#   ### Extrapolation and prediction interval
# 
#   ## Sample dose points on the same range as the original aerobic experimental data
#   set.seed(k)
#   Dnew.O2 = runif(150, min = 0, max = max(D.O2))
# 
#   ## Generate observation points with Gaussian noise around the best-fit curve
#   ## predicted by the LQ and USC model in the dose range of Dnew.O2
#   set.seed(k)
#   SFnew.LQO2 = LQfunc(alpha.O2, beta.O2, Dnew.O2) + rnorm(150, 0, sd(SF.O2))
#   set.seed(k)
#   SFnew.USCO2 = USCfunc(alpha.O2, beta.O2, Dnew.O2, DT.O2) + rnorm(150, 0, sd(SF.O2))
# 
#   ## Sample dose points on the same range as the original extremely-hypoxic experimental data
#   set.seed(k)
#   Dnew.N2 = runif(150, min = 0, max = max(D.N2))
# 
#   ## Generate observation points with Gaussian noise around the best-fit curve
#   ## predicted by the LQ and USC model in the dose range of Dnew.N2
#   set.seed(k)
#   SFnew.LQN2 = LQfunc(alpha.N2, beta.N2, Dnew.N2) + rnorm(150, 0, sd(SF.N2))
#   set.seed(k)
#   SFnew.USCN2 = USCfunc(alpha.N2, beta.N2, Dnew.N2, DT.N2) + rnorm(150, 0, sd(SF.N2))
# 
#   ## Solve least squares and find the best fit for both LQ and USC
#   newLQfit.O2 = computeFit(Dnew.O2, SFnew.LQO2)
#   newLQfit.N2 = computeFit(Dnew.N2, SFnew.LQN2)
#   newUSCfit.O2 = computeFit(Dnew.O2, SFnew.USCO2, DT.O2) # find.DT(Dnew.O2, SFnew.USCO2, seq(DTrange.O2[2*k-1], DTrange.O2[2*k], 0.01) )
#   newUSCfit.N2 = computeFit(Dnew.N2, SFnew.USCN2, DT.N2) # find.DT(Dnew.N2, SFnew.USCN2, seq(DTrange.N2[2*k-1], DTrange.N2[2*k], 0.01) )
# 
#   ## Compute 95% CI for the extended dose range and plot the predictions
#   Dtest.O2 = seq(from = 0, to = Doserange.O2[k], by = 0.5)
#   Dtest.N2 = seq(from = 0, to = Doserange.N2[k], by = 0.5)
#   p.LQO2 = CI(newLQfit.O2, Dtest.O2)
#   p.LQN2 = CI(newLQfit.N2, Dtest.N2)
#   p.USCO2 = CI(newUSCfit.O2, Dtest.O2)
#   p.USCN2 = CI(newUSCfit.N2, Dtest.N2)
# 
#   dev.new()
#   matplot(Dtest.N2, exp(p.LQN2[,-2]),
#           ylab = "Surviving fraction", xlab = "Dose [Gy]",
#           xlim = c(0, Doserange.N2[k]), ylim = c(0.0000001, 1),
#           log = "y", las = 1, col = c(1,4,4), lty = c(1,2,2), type = "l", lwd = 1, yaxt = "n",
#           cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
#   ticks = seq(-7, 0, by=1)
#   labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
#   axis(2, at=c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1), labels=labels,
#        cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
#   matlines( Dtest.O2, exp(p.LQO2[,-2]), col = c(1,4,4), lty = c(1,2,2), type = "l", lwd = 1 )
#   matlines( Dtest.N2, exp(p.USCN2[,-2]), col = c(1,2,2), lty = c(3,2,2), type = "l", lwd = 1 )
#   matlines( Dtest.O2, exp(p.USCO2[,-2]), col = c(1,2,2), lty = c(3,2,2), type = "l", lwd = 1 )
#   points(Dnew.O2, exp(SFnew.LQO2), pch = 1, cex = 0.5)
#   points(Dnew.N2, exp(SFnew.LQN2), pch = 16, cex = 0.5)
#   points(Dnew.O2, exp(SFnew.USCO2), pch = 1, cex = 0.5)
#   points(Dnew.N2, exp(SFnew.USCN2), pch = 16, cex = 0.5)
#   legend("topright", legend=c(legend.title[k], "Aerobic", "Extremely hypoxic", "LQ", "USC", "LQ 95% CI", "USC 95% CI"),
#          lty = c(NA, 0, 0, 1, 3, 2, 2), pch = c(NA, 1, 16, NA, NA, NA, NA), col = c(NA, 1, 1, 1, 1, 4, 2), lwd = 2, bty = "o")
# 
#   ## Print model coeffisients and their corresponding SD
#   SDalpha.O2 = coef(summary(LQmod.O2))[1, 2]
#   SDalpha.N2 = coef(summary(LQmod.N2))[1, 2]
#   SDbeta.O2 = coef(summary(LQmod.O2))[2, 2]
#   SDbeta.N2 = coef(summary(LQmod.N2))[2, 2]
#   SDv.O2 = coef(summary(MTSHmod.O2))[1, 2]
#   SDv.N2 = coef(summary(MTSHmod.N2))[1, 2]
#   SDm.O2 = coef(summary(MTSHmod.O2))[2, 2]
#   SDm.N2 = coef(summary(MTSHmod.N2))[2, 2]
#   SDD0.O2 = D0.O2^2*SDv.O2
#   SDDq.O2 = sqrt( (SDD0.O2*log(m.O2))^2 + (SDm.O2*D0.O2/m.O2)^2 )
#   SDD0.N2 = D0.N2^2*SDv.N2
#   SDDq.N2 = sqrt( (SDD0.N2*log(m.N2))^2 + (SDm.N2*D0.N2/m.N2)^2 )
#   SDDT.O2 = (DT.O2/Dq.O2)*sqrt( SDDq.O2^2 + (D0.O2*DT.O2*SDalpha.O2/2)^2 + (alpha.O2*DT.O2*SDD0.O2/2)^2 )
#   SDDT.N2 = (DT.N2/Dq.N2)*sqrt( SDDq.N2^2 + (D0.N2*DT.N2*SDalpha.N2/2)^2 + (alpha.N2*DT.N2*SDD0.N2/2)^2 )
# 
#   #SDD0.O2 = coef(summary(MTSHmod.O2))[1, 2]
#   #SDD0.N2 = coef(summary(MTSHmod.N2))[1, 2]
#   #SDDq.O2 = coef(summary(MTSHmod.O2))[2, 2]
#   #SDDq.N2 = coef(summary(MTSHmod.N2))[2, 2]
# 
#   alpha.ratio = alpha.O2/alpha.N2
#   beta.ratio = sqrt(beta.O2/beta.N2)
#   alphabeta.ratio.O2 = alpha.O2/beta.O2
#   alphabeta.ratio.N2 = alpha.N2/beta.N2
# 
#   SDalpha.ratio = alpha.ratio*sqrt( (SDalpha.O2/alpha.O2)^2 + (SDalpha.N2/alpha.N2)^2 )
#   SDbeta.ratio = 0.5*beta.ratio*sqrt( (SDbeta.O2/beta.O2)^2 + (SDbeta.N2/beta.N2)^2 )
#   SDalphabeta.ratio.O2 = alphabeta.ratio.O2*sqrt( (SDalpha.O2/alpha.O2)^2 + (SDbeta.O2/beta.O2)^2 )
#   SDalphabeta.ratio.N2 = alphabeta.ratio.N2*sqrt( (SDalpha.N2/alpha.N2)^2 + (SDbeta.N2/beta.N2)^2 )
#   
#   ntot = length(D.O2) + length(D.N2)
#   t_value_OER1 = abs(alpha.ratio - beta.ratio)/sqrt(SDalpha.ratio^2/ntot + SDbeta.ratio^2/ntot)
#   df_OER1 = (SDalpha.ratio^2/ntot + SDbeta.ratio^2/ntot)^2/( (SDalpha.ratio^2 /ntot)^2/(ntot-1) + (SDbeta.ratio^2 /ntot)^2/(ntot-1) ) 
#   cat("t-value: ", t_value_OER1, "\n")
#   cat("df: ", df_OER1, "\n")
#   cat("n: ", ntot, "\n")
#   
#   write(legend.title[k], file = "NHIK3025.txt", append = T)
#   write(c("alpha (O2)", "beta (O2)", "alpha (N2)", "beta (N2)",
#           "D0 (O2)", "Dq (O2)", "D0 (N2)", "Dq (N2)",
#           "alpha_O2/alpha_N2", "sqrt(beta_O2/beta_N2)",
#           "DT (O2)", "DT (N2)", "alpha/beta (O2)", "alpha/beta (N2)"),
#         file = "NHIK3025.txt", ncolumns = 14, sep = "\t\t\t", append = T)
#   write(c( (-1)*coef(summary(LQmod.O2))[, 1], (-1)*coef(summary(LQmod.N2))[, 1],
#            D0.O2, Dq.O2, D0.N2, Dq.N2,
#            alpha.ratio, beta.ratio, DT.O2, DT.N2, alphabeta.ratio.O2, alphabeta.ratio.N2),
#         file = "NHIK3025.txt", ncolumns = 14, sep = "\t\t\t", append = T)
#   write(c( coef(summary(LQmod.O2))[, 2], coef(summary(LQmod.N2))[, 2],
#            SDD0.O2, SDDq.O2, SDD0.N2, SDDq.N2,
#            SDalpha.ratio, SDbeta.ratio, SDDT.O2, SDDT.N2, SDalphabeta.ratio.O2, SDalphabeta.ratio.N2),
#         file = "NHIK3025.txt", ncolumns = 14, sep = "\t\t\t", append = T)
#   write( c( "Adjusted R-sqaured (O2): ", summary(LQmod.O2)$adj.r.squared, R.squared.adj(log(SF.O2), USC.O2, length(D.O2), 4) ),
#          file = "NHIK3025.txt", ncolumns = 3, sep = "\t", append = T)
#   write( c( "Adjusted R-squared (N2): ", summary(LQmod.N2)$adj.r.squared, R.squared.adj(log(SF.N2), USC.N2, length(D.N2), 4) ),
#          file = "NHIK3025.txt", ncolumns = 3, sep = "\t", append = T)
#   write( c( "RSS (O2): ", RSS(log(SF.O2), fitted(MTSHmod.O2)), RSS(log(SF.O2), fitted(LQmod.O2)), RSS(log(SF.O2), USC.O2)),
#          file = "NHIK3025.txt", ncolumns = 4, sep = "\t\t\t", append = T)
#   write( c( "RSS (N2): ", RSS(log(SF.N2), fitted(MTSHmod.N2)), RSS(log(SF.N2), fitted(LQmod.N2)), RSS(log(SF.N2), USC.N2)),
#          file = "NHIK3025.txt", ncolumns = 4, sep = "\t\t\t", append = T)
#   write("-------------------------------------------------------------------------------------------------------------",
#         file = "NHIK3025.txt", append = T)
  
  
  # if (i == "M") {
  #   
  #   dev.new()
  #   matplot(Dtest.O2, exp(p.LQO2[,-2]),
  #           ylab = "Surviving fraction", xlab = "Dose [Gy]",
  #           xlim = c(0, Doserange.O2[k]), ylim = c(0.0000001, 1),
  #           log = "y", las = 1, col = c(1,4,4), lty = c(1,2,2), type = "l", lwd = 1, yaxt = "n",
  #           cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
  #   ticks = seq(-7, 0, by=1)
  #   labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  #   axis(2, at=c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1), labels=labels,
  #        cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
  #   matlines( Dtest.O2, exp(p.LQO2[,-2]), col = c(1,4,4), lty = c(1,2,2), type = "l", lwd = 1 )
  #   matlines( Dtest.O2, exp(p.USCO2[,-2]), col = c(1,2,2), lty = c(3,2,2), type = "l", lwd = 1 )
  #   points(Dnew.O2, exp(SFnew.LQO2), pch = 1, cex = 0.5)
  #   points(Dnew.O2, exp(SFnew.USCO2), pch = 1, cex = 0.5)
  #   legend("topright", legend=c(legend.title[k], "Aerobic", "LQ", "USC", "LQ 95% CI", "USC 95% CI"),
  #          lty = c(NA, 0, 1, 3, 2, 2), pch = c(NA, 1, NA, NA, NA, NA), col = c(NA, 1, 1, 1, 4, 2), lwd = 2, bty = "o")
  #   
  # }
  # 
  
  
  
}


## All-encompassing model

SF = SFdata[Cellcycle == "S" & Hypoxia == "Yes", 1]
D = SFdata[Cellcycle == "S" & Hypoxia == "Yes", 2]

GLM = glm(SF ~ -1 + Dose..Gy. + I(Dose..Gy.^2) + Dose..Gy.*Cellcycle + Dose..Gy.*Hypoxia  + Cellcycle:Hypoxia, family = quasipoisson(link = "log"))


model = nls(SF ~ exp(-(OMF*alpha*D + beta*OMF^2*I(D^2))), start = list(OMF = 1., alpha = 0.15, beta = 0.001))



