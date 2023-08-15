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
  tab = matrix(c(beta, se, tstat, pval), nrow = p, ncol = 4L,
                dimnames = list(dimnames(X)[[2L]], 
                                c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
  
  ## 2 * negative log-likelihood
  # nega2logLik = n * log(2 * pi * sig2) + (n - p)
  ## Compute AIC. k = p + 1 for the USC since the model is determined by alpha, beta and DT.
  AIC = n*log(2*pi*sig2) + (n - p) + 2*(p + 1) 
  
  ## Compute R-squared and adjusted R-squared
  TSS = c(crossprod(SF - sum(SF) / n))
  r.squared = 1 - RSS/TSS
  adj.r.squared = 1 - sig2*(n - 1)/TSS
  
  ## Return estimated values
  if ( is.null(DT) == TRUE ) {
    list(coef = beta, residuals = residuals, fitted.values = c(X %*% beta),
         R = R, sig2 = sig2, se = se, coef.table = tab, AIC = AIC,
         RSS = RSS, r.squared = r.squared, adj.r.squared = adj.r.squared, p = p, n = n)
  } else {
    list(coef = beta, residuals = residuals, fitted.values = c(X %*% beta),
         R = R, sig2 = sig2, se = se, coef.table = tab, AIC = AIC, DT = DT,
         RSS = RSS, r.squared = r.squared, adj.r.squared = adj.r.squared, p = p, n = n)
  }
  
}

## Compute the transition dose DT as an independent parameter
find.DT = function (D, SF, DTspan) {
  
  # Regression for each value of DT in span 
  lst = lapply(DTspan, computeFit, D = D, SF = SF)
  
  # Compute RSS for each value of DT and plot RSS vs grid of DT
  RSS = sapply(lst, "[[", "RSS")
  dev.new()
  plot(DTspan, RSS, type = "b", pch = 1)
  
  # Find DT that minimizes RSS
  # minindex = lst[[which( RSS == min(RSS) )]]
  # print( DTspan[[length(minindex)]]  )
  lst[[which.min(RSS)]] # lst[[length(minindex)]]
  
}

## Compute the 95% CI for the fitted y values
CI = function (model, D.new) {
  
  ## Prediction matrix
  if ( is.null(model$DT) == TRUE ) {
    X = cbind( "alpha" = D.new, "beta" = D.new^2 )
  } else {
    X = cbind( "alpha" = D.new, "beta" = f(D.new,model$DT) )
  }
  
  ## Number of predictor variables, i.e. the dimension of the model matrix
  p = dim(X)[2L]
  
  ## Fitted values or expected mean
  fit = X %*% model$coef
  
  ## Prediction standard error
  Qt = forwardsolve(t(model$R), t(X))
  se = sqrt(colSums(Qt^2)*model$sig2)
  
  ## 95% CI
  alpha = qt(0.025, length(model$residuals) - p)
  lwr = fit + alpha * se
  upr = fit - alpha * se
  ## 95%-prediction interval
  #lwr = fit + alpha * sqrt(model$sig2 + I(se^2))
  #upr = fit + alpha - sqrt(model$sig2 + I(se^2))
  #lwr = fit + tval*s(y,fit,model$n,model$p)*sqrt( 1/(model$n-model$p) + I((x.new - mean(x.new))^2)/Sxx(x.new) )
  #upr = fit - tval*s(y,fit,model$n,model$p)*sqrt( 1/(model$n-model$p) + I((x.new - mean(x.new))^2)/Sxx(x.new) )
  
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
v.N2.vec = c(0.5, 0.2, 0.2, 0.2)
m.N2.vec = c(2, 2, 2, 2)
D0.O2.vec = c(0.84, 1.08, 0.82, 0.86)
D0.N2.vec = c(0.5, 3.46, 2.55, 2.36)
Dq.O2.vec = c(0.29, 2.30, 1.66, 1.52)
Dq.N2.vec = c(3, 6.16, 5.90, 8.57)

#Cellcycle = relevel(Cellcycle, ref = "M") # change referanse/baseline level 

## Initiate the vectors
Dlist.O2 = Dlist.N2 = SFlist.O2 = SFlist.N2 = vector("list", 4)
vecDT.O2 = vecDT.N2 = vecalpha.O2 = vecbeta.O2 = vecalpha.N2 = vecbeta.N2 = vector()

## Loop through the different cell populations 
k = 0
for (i in c("M", "G1", "S", "Asynch") ) {
  
  k = k + 1
  
  ## Extract dose-response observations
  SF.O2 = SFdata[Cellcycle == i & Hypoxia == "No", 1]
  SF.N2 = SFdata[Cellcycle == i & Hypoxia == "Yes", 1]
  D.O2 = SFdata[Cellcycle == i & Hypoxia == "No", 2]
  D.N2 = SFdata[Cellcycle == i & Hypoxia == "Yes", 2]
  SD.O2 = SFdata[Cellcycle == i & Hypoxia == "No", 5]
  SD.N2 = SFdata[Cellcycle == i & Hypoxia == "Yes", 5]
  SFmean.O2 = SFdata[Cellcycle == i & Hypoxia == "No", 6]
  SFmean.N2 = SFdata[Cellcycle == i & Hypoxia == "Yes", 6]

  ## Linear and non-linear regression of LQ and MTSH, respectively
  LQmod.O2 = lm(log(SF.O2) ~ -1 + D.O2 + I(D.O2^2))
  LQmod.N2 = lm(log(SF.N2) ~ -1 + D.N2 + I(D.N2^2))
  MTSHmod.O2 = nls(log(SF.O2) ~ log(1 - (1 - exp(-(1/D0)*D.O2))^(exp(Dq/D0)) ), start = list(D0=D0.O2.vec[k] , Dq=Dq.O2.vec[k]))
  MTSHmod.N2 = nls(log(SF.N2) ~ log(1 - (1 - exp(-(1/D0)*D.N2))^(exp(Dq/D0))  ), start = list(D0=D0.N2.vec[k], Dq=Dq.N2.vec[k]))
  
  ## Model coefficients
  alpha.N2 = (-1)*coef(LQmod.N2)[1]
  beta.N2 = (-1)*coef(LQmod.N2)[2]
  D0.N2 = coef(MTSHmod.N2)[1]
  Dq.N2 = coef(MTSHmod.N2)[2]
  DT.N2 = (2*Dq.N2)/(1 - alpha.N2*D0.N2)
  
  alpha.O2 = (-1)*coef(LQmod.O2)[1]
  beta.O2 = (-1)*coef(LQmod.O2)[2]
  D0.O2 = coef(MTSHmod.O2)[1]
  Dq.O2 = coef(MTSHmod.O2)[2]
  DT.O2 = (2*Dq.O2)/(1 - alpha.O2*D0.O2)
  
  ## Compute USC curve 
  USC.O2 = USCfunc(alpha.O2, beta.O2, D.O2, DT.O2)
  USC.N2 = USCfunc(alpha.N2, beta.N2, D.N2, DT.N2)
 
  # ## Plot dose-response curves for the cell population under O2obic and extremely-N2oxic cells
  # dev.new()
  # plot(D.O2, SF.O2,
  #      ylab = "Surviving fraction", xlab = "Dose [Gy]",
  #      ylim = c(min(SF.O2), 1), xlim = c(0, max(D.N2)),
  #      las = 1, pch = 1, log = "y")
  # points(D.N2, SF.N2, pch = 16)
  # lines(sort(D.O2), exp(fitted(LQmod.O2))[order(D.O2)], lty = 1)
  # lines(sort(D.N2), exp(fitted(LQmod.N2))[order(D.N2)], lty = 1)
  # lines(sort(D.O2), exp(fitted(MTSHmod.O2))[order(D.O2)], lty = 2)
  # lines(sort(D.N2), exp(fitted(MTSHmod.N2))[order(D.N2)], lty = 2)
  # lines(sort(D.O2), exp(USC.O2)[order(D.O2)], lty = 3)
  # lines(sort(D.N2), exp(USC.N2)[order(D.N2)], lty = 3)
  # legend("topright", legend=c(legend.title[k], "Aerobic", "Extremely hypoxic", "LQ", "MTSH", "USC"),
  #        lty = c(0, 0, 0, 1, 2, 3), pch = c(NA, 1, 16, NA, NA, NA), bty = "o")
  # 
  # ## Compute model residuals
  # LQresid.O2 = resid(LQmod.O2)
  # LQresid.N2 = resid(LQmod.N2)
  # MTSHresid.O2 = resid(MTSHmod.O2) # standardresid(D.O2, log(SF.O2), fitted(MTSHmod.O2), length(resid(MTSHmod.O2)), 2)
  # MTSHresid.N2 = resid(MTSHmod.N2) # standardresid(D.N2, log(SF.N2), fitted(MTSHmod.N2), length(resid(MTSHmod.N2)), 2)
  # USCresid.O2 = log(SF.O2) - USC.O2 # standardresid(D.O2, log(SF.O2), USC.O2, length(USC.O2), 2) # log(SF.O2) - USC.O2
  # USCresid.N2 = log(SF.N2) - USC.N2 # standardresid(D.N2, log(SF.N2), USC.N2, length(USC.N2), 2)
  # 
  # ## Aerobic residual plot
  # dev.new()
  # plot(exp(fitted(LQmod.O2)), LQresid.O2,
  #      ylab = "Residuals", xlab = "Estimated surviving fraction (fitted values)",
  #      ylim = c(-2, 2), xlim = c(min(SF.O2), 1),
  #      log = "x", las = 1, pch = 1)
  # points(exp(fitted(MTSHmod.O2)), MTSHresid.O2, pch = 2)
  # points(exp(USC.O2), USCresid.O2, pch = 3)
  # abline(0, 0, lty = 2)
  # 
  # ## Add trend curves in the aerobic residual plot for all three models
  # LOESS.MTSHresid.O2 = loess(MTSHresid.O2 ~ fitted(MTSHmod.O2) )
  # LOESS.LQresid.O2 = loess(LQresid.O2 ~ (fitted(LQmod.O2)) )
  # LOESS.USCresid.O2 = loess(USCresid.O2 ~ USC.O2)
  # 
  # lines(sort(exp(fitted(LQmod.O2))), fitted(LOESS.LQresid.O2)[order(exp(fitted(LQmod.O2)))], lty = 1)
  # lines(sort(exp(fitted(MTSHmod.O2))), fitted(LOESS.MTSHresid.O2)[order(exp(fitted(MTSHmod.O2)))], lty = 2)
  # lines(sort(exp(USC.O2)), fitted(LOESS.USCresid.O2)[order(exp(USC.O2))], lty = 3)
  # legend("topright", legend = c(paste(legend.title[k], "- Aerobic"), "LQ", "MTSH", "USC"),
  #        lty = c(NA, 1, 2, 3), pch = c(NA, 1, 2, 3), bty = "o")
  #  
  # ## Extremely-hypoxic residual plot
  # dev.new()
  # plot(exp(fitted(LQmod.N2)), LQresid.N2,
  #      ylab = "Residuals", xlab = "Estimated surviving fraction (fitted values)",
  #      ylim = c(-2, 2), xlim = c(min(SF.N2), 1),
  #      log = "x", las = 1, pch = 1)
  # points(exp(fitted(MTSHmod.N2)), MTSHresid.N2, pch = 2)
  # points(exp(USC.N2), USCresid.N2, pch = 3)
  # abline(0, 0, lty = 2)
  # 
  # ## Add trend curves in the extremely-hypoxic residual plot for all models
  # LOESS.LQresid.N2 = loess(LQresid.N2 ~ (fitted(LQmod.N2)) )
  # LOESS.MTSHresid.N2 = loess(MTSHresid.N2 ~ fitted(MTSHmod.N2) )
  # LOESS.USCresid.N2 = loess(USCresid.N2 ~ USC.N2)
  # 
  # lines(sort(exp(fitted(LQmod.N2))), fitted(LOESS.LQresid.N2)[order(exp(fitted(LQmod.N2)))], lty = 1)
  # lines(sort(exp(fitted(MTSHmod.N2))), fitted(LOESS.MTSHresid.N2)[order(exp(fitted(MTSHmod.N2)))], lty = 2)
  # lines(sort(exp(USC.N2)), fitted(LOESS.USCresid.N2)[order(exp(USC.N2))], lty = 3)
  # legend("topright", legend = c(paste(legend.title[k], "- Extremely hypoxic"), "LQ", "MTSH", "USC"),
  #        lty = c(NA, 1, 2, 3), pch = c(NA, 1, 2, 3), bty = "o")
  # 
  
  
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
  plot(DTspan.O2, RSS.O2, 
       ylab = "RSS", xlab = "DT grid [Gy]", 
       ylim = c(min(min(RSS.O2), min(RSS.N2)), max(max(RSS.O2), max(RSS.N2))), 
       xlim = c(min(min(DTspan.O2), min(DTspan.N2)), max(max(DTspan.O2), max(DTspan.N2))),
       type = "o", pch = c(1, rep(NA, 24)) )
  lines(DTspan.N2, RSS.N2, type = "o", pch = c(16, rep(NA, 24)))
  legend("topright", legend=c(legend.title[k], "Aerobic", "Extremely hypoxic"),
         pch = c(NA, 1, 16), bty = "o")
  
  
  ### Extrapolation and prediction interval
  
  ## Sample dose points on the same range as the original aerobic experimental data
  set.seed(k)
  Dnew.O2 = runif(150, min = 0, max = max(D.O2))

  ## Generate observation points with Gaussian noise around the best-fit curve
  ## predicted by the LQ and USC model in the dose range of Dnew.O2
  set.seed(k)
  SFnew.LQO2 = LQfunc(alpha.O2, beta.O2, Dnew.O2) + rnorm(150, 0, 0.2*sd(log(SF.O2)))
  set.seed(k)
  SFnew.USCO2 = USCfunc(alpha.O2, beta.O2, Dnew.O2, DT.O2) + rnorm(150, 0, 0.2*sd(log(SF.O2)))
  
  ## Sample dose points on the same range as the original extremely-hypoxic experimental data
  set.seed(k)
  Dnew.N2 = runif(150, min = 0, max = max(D.N2))

  ## Generate observation points with Gaussian noise around the best-fit curve
  ## predicted by the LQ and USC model in the dose range of Dnew.N2
  set.seed(k)
  SFnew.LQN2 = LQfunc(alpha.N2, beta.N2, Dnew.N2) + rnorm(150, 0, 0.2*sd(log(SF.N2)))
  set.seed(k)
  SFnew.USCN2 = USCfunc(alpha.N2, beta.N2, Dnew.N2, DT.N2) + rnorm(150, 0, 0.2*sd(log(SF.N2)))

  ## Solve least squares and find the best fit for both LQ and USC
  newLQfit.O2 = computeFit(Dnew.O2, SFnew.LQO2)
  newLQfit.N2 = computeFit(Dnew.N2, SFnew.LQN2)
  newUSCfit.O2 = computeFit(Dnew.O2, SFnew.USCO2, DT.O2) # find.DT(Dnew.O2, SFnew.USCO2, seq(DTrange.O2[2*k-1], DTrange.O2[2*k], 0.01) )
  newUSCfit.N2 = computeFit(Dnew.N2, SFnew.USCN2, DT.N2) # find.DT(Dnew.N2, SFnew.USCN2, seq(DTrange.N2[2*k-1], DTrange.N2[2*k], 0.01) )
  
  ## Compute 95% CI for the extended dose range and plot the predictions
  Dtest.O2 = seq(from = 0, to = Doserange.O2[k], by = 0.5)
  Dtest.N2 = seq(from = 0, to = Doserange.N2[k], by = 0.5)
  p.LQO2 = CI(newLQfit.O2, Dtest.O2)
  p.LQN2 = CI(newLQfit.N2, Dtest.N2)
  p.USCO2 = CI(newUSCfit.O2, Dtest.O2)
  p.USCN2 = CI(newUSCfit.N2, Dtest.N2)
  
  dev.new()
  matplot(Dtest.N2, exp(p.LQN2[,-2]),
          ylab = "Surviving fraction", xlab = "Dose [Gy]",
          xlim = c(0, Doserange.N2[k]), ylim = c(0.0000001, 1),
          log = "y", las = 1, col = c(1,2,2), lty = c(1,2,2), type = "l", lwd = 1 )
  matlines( Dtest.O2, exp(p.LQO2[,-2]), col = c(1,2,2), lty = c(1,2,2), type = "l", lwd = 1 )
  matlines( Dtest.N2, exp(p.USCN2[,-2]), col = c(1,4,4), lty = c(3,2,2), type = "l", lwd = 1 )
  matlines( Dtest.O2, exp(p.USCO2[,-2]), col = c(1,4,4), lty = c(3,2,2), type = "l", lwd = 1 )
  points(Dnew.O2, exp(SFnew.LQO2), pch = 1, cex = 0.5)
  points(Dnew.N2, exp(SFnew.LQN2), pch = 16, cex = 0.5)
  points(Dnew.O2, exp(SFnew.USCO2), pch = 1, cex = 0.5)
  points(Dnew.N2, exp(SFnew.USCN2), pch = 16, cex = 0.5)
  legend("topright", legend=c(legend.title[k], "LQ", "USC", "LQ 95% CI", "USC 95% CI"),
         lty = c(NA, 1, 3, 2, 2), col = c(NA, 1, 1, 2, 4), lwd = 2, bty = "o")
  #
  # Dlist.O2[[k]] = D.O2
  # Dlist.N2[[k]] = D.N2
  # SFlist.O2[[k]] = SF.O2
  # SFlist.N2[[k]] = SF.N2
  # 
  # vecalpha.O2 = c(vecalpha.O2, alpha.O2)
  # vecalpha.N2 = c(vecalpha.N2, alpha.N2)
  # vecbeta.O2 = c(vecbeta.O2, beta.O2)
  # vecbeta.N2 = c(vecbeta.N2, beta.N2)
  # vecDT.O2 = c(vecDT.O2, DT.O2)
  # vecDT.N2 = c(vecDT.N2, DT.N2)

  
   ## Print model coeffisients and their corresponding SD
  # SDalpha.O2 = coef(summary(LQmod.O2))[1, 2]
  # SDalpha.N2 = coef(summary(LQmod.N2))[1, 2]
  # SDbeta.O2 = coef(summary(LQmod.O2))[2, 2]
  # SDbeta.N2 = coef(summary(LQmod.N2))[2, 2]
  # #SDv.O2 = coef(summary(MTSHmod.O2))[1, 2]
  # #SDv.N2 = coef(summary(MTSHmod.N2))[1, 2]
  # #SDm.O2 = coef(summary(MTSHmod.O2))[2, 2]
  # #SDm.N2 = coef(summary(MTSHmod.N2))[2, 2]
  # SDD0.O2 = coef(summary(MTSHmod.O2))[1, 2]
  # SDD0.N2 = coef(summary(MTSHmod.N2))[1, 2]
  # SDDq.O2 = coef(summary(MTSHmod.O2))[2, 2]
  # SDDq.N2 = coef(summary(MTSHmod.N2))[2, 2]
  # 
  # alpha.ratio = alpha.O2/alpha.N2
  # beta.ratio = sqrt(beta.O2/beta.N2)
  # alphabeta.ratio.O2 = alpha.O2/beta.O2
  # alphabeta.ratio.N2 = alpha.N2/beta.N2
  # 
  # SDalpha.ratio = alpha.ratio*sqrt( (SDalpha.O2/alpha.O2)^2 + (SDalpha.N2/alpha.N2)^2 )
  # SDbeta.ratio = 0.5*beta.ratio*sqrt( (SDbeta.O2/beta.O2)^2 + (SDbeta.N2/beta.N2)^2 )
  # #SDD0.O2 = D0.O2^2*SDv.O2
  # #SDDq.O2 = sqrt( (SDD0.O2*log(m.O2))^2 + (SDm.O2*D0.O2/m.O2)^2 )
  # SDDT.O2 = (DT.O2/Dq.O2)*sqrt( SDDq.O2^2 + (D0.O2*DT.O2*SDalpha.O2/2)^2 + (alpha.O2*DT.O2*SDD0.O2/2)^2 )
  # #SDD0.N2 = D0.N2^2*SDv.N2
  # #SDDq.N2 = sqrt( (SDD0.N2*log(m.N2))^2 + (SDm.N2*D0.N2/m.N2)^2 )
  # SDDT.N2 = (DT.N2/Dq.N2)*sqrt( SDDq.N2^2 + (D0.N2*DT.N2*SDalpha.N2/2)^2 + (alpha.N2*DT.N2*SDD0.N2/2)^2 )
  # SDalphabeta.ratio.O2 = alphabeta.ratio.O2*sqrt( (SDalpha.O2/alpha.O2)^2 + (SDbeta.O2/beta.O2)^2 )
  # SDalphabeta.ratio.N2 = alphabeta.ratio.N2*sqrt( (SDalpha.N2/alpha.N2)^2 + (SDbeta.N2/beta.N2)^2 )
  # 
  # write(legend.title[k], file = "NHIK3025.txt", append = T)
  # write(c("alpha (O2)", "beta (O2)", "alpha (N2)", "beta (N2)",
  #         "D0 (O2)", "Dq (O2)", "D0 (N2)", "Dq (N2)",
  #         "alpha_O2/alpha_N2", "sqrt(beta_O2/beta_N2)",
  #         "DT (O2)", "DT (N2)", "alpha/beta (O2)", "alpha/beta (N2)"),
  #       file = "NHIK3025.txt", ncolumns = 14, sep = "\t\t\t", append = T)
  # write(c( (-1)*coef(summary(LQmod.O2))[, 1], (-1)*coef(summary(LQmod.N2))[, 1],
  #          coef(summary(MTSHmod.O2))[, 1], coef(summary(MTSHmod.N2))[, 1],
  #          alpha.ratio, beta.ratio, DT.O2, DT.N2, alphabeta.ratio.O2, alphabeta.ratio.N2),
  #       file = "NHIK3025.txt", ncolumns = 14, sep = "\t\t\t", append = T)
  # write(c( coef(summary(LQmod.O2))[, 2], coef(summary(LQmod.N2))[, 2],
  #          coef(summary(MTSHmod.O2))[, 2], coef(summary(MTSHmod.N2))[, 2],
  #          SDalpha.ratio, SDbeta.ratio, SDDT.O2, SDDT.N2, SDalphabeta.ratio.O2, SDalphabeta.ratio.N2),
  #       file = "NHIK3025.txt", ncolumns = 14, sep = "\t\t\t", append = T)
  # write( c( "Adjusted R-sqaured (O2): ", summary(LQmod.O2)$adj.r.squared, R.squared.adj(log(SF.O2), USC.O2, length(D.O2), 2) ),
  #        file = "NHIK3025.txt", ncolumns = 3, sep = "\t", append = T)
  # write( c( "Adjusted R-squared (N2): ", summary(LQmod.N2)$adj.r.squared, R.squared.adj(log(SF.N2), USC.N2, length(D.N2), 2) ),
  #        file = "NHIK3025.txt", ncolumns = 3, sep = "\t", append = T)
  # write( c( "RSS (O2): ", RSS(log(SF.O2), fitted(MTSHmod.O2)), RSS(log(SF.O2), fitted(LQmod.O2)), RSS(log(SF.O2), USC.O2)),
  #        file = "NHIK3025.txt", ncolumns = 4, sep = "\t\t\t", append = T)
  # write( c( "RSS (N2): ", RSS(log(SF.N2), fitted(MTSHmod.N2)), RSS(log(SF.N2), fitted(LQmod.N2)), RSS(log(SF.N2), USC.N2)),
  #        file = "NHIK3025.txt", ncolumns = 4, sep = "\t\t\t", append = T)
  # write("-------------------------------------------------------------------------------------------------------------",
  #       file = "NHIK3025.txt", append = T)

  ## Isoeffect curves
  # nref = 45
  # dref = 1.8
  # Dref = nref*dref # EQD2
  # n = 3
  # d = seq(from = 2, to = 30, by = 1)
  # 
  # Diso.LQ.O2 = Dref*(dref + alpha.O2/beta.O2 )/( d + alpha.O2/beta.O2 ) # + rnorm(length(d), sd = 1)
  # Diso.LQ.N2 = Dref*(dref + alpha.N2/beta.N2 )/( d + alpha.N2/beta.N2 ) # + rnorm(length(d), sd = 1)
  # Diso.USC.O2 = Dref*( alpha.O2 + (beta.O2/dref)*I(f(dref, DT.O2)) )/( alpha.O2 + (beta.O2/d)*I(f(d, DT.O2)) )
  # Diso.USC.N2 = Dref*( alpha.N2 + (beta.N2/dref)*I(f(dref, DT.N2)) )/( alpha.N2 + (beta.N2/d)*I(f(d, DT.N2)) )
  # 
  # dev.new()
  # plot(d, Diso.LQ.O2,
  #      ylab = "Total dose [Gy]", xlab = "Dose/fraction [Gy]",
  #      ylim = c( min(Diso.LQ.O2), max(Diso.LQ.N2)), xlim = rev(range(d)),
  #      las = 1, log = "x", type = "l", lty = 1, col = "Red")
  # points(d, Diso.LQ.N2, type = "l", lty = 1, col = "Blue")
  # points(d, Diso.USC.O2, type = "l", lty = 2, col = "Red")
  # points(d, Diso.USC.N2, type = "l", lty = 2, col = "Blue")
  # legend("topleft", legend=c(legend.title[k], "LQ, O2obic", "USC, O2obic", "LQ, extremely N2oxic", "USC, extremely N2oxic"),
  #        lty = c(0, 1, 2, 1, 2), col = c(0, "Red", "Red", "Blue", "Blue"), bty = "o")

}

# -------------------------------------------------------------------------------------------------------------- #
# G2
SF.G2.N2 = SFdata[Cellcycle == "G2" & N2oxia == "Yes", 1]
D.G2.N2 = SFdata[Cellcycle == "G2" & N2oxia == "Yes", 2]

LQmod.G2.N2 = lm(log(SF.G2.N2) ~ -1 + D.G2.N2 + I(D.G2.N2^2))
MTSHmod.G2.N2 = nls(log(SF.G2.N2) ~ log(1 - (1 - exp(-v*D.G2.N2) )^m ), start = list(v = 1, m = 1))

# USC
alpha.G2.N2 = (-1)*coef(LQmod.G2.N2)[1]
beta.G2.N2 = (-1)*coef(LQmod.G2.N2)[2]
D0.G2.N2 = 1/coef(MTSHmod.G2.N2)[1]
Dq.G2.N2 = log(coef(MTSHmod.G2.N2)[2])*D0.G2.N2
DT.G2.N2 = (2*Dq.G2.N2)/(1 - alpha.G2.N2*D0.G2.N2)

USC.G2.N2 = -alpha.G2.N2*D.G2.N2 - beta.G2.N2*I( f(D.G2.N2, DT.G2.N2) )

#USC.G2.N2 = lm(log(SF.G2.N2) ~ -1 + D.G2.N2 + I( f(D.G2.N2, DT.G2.N2) ))
# USCalpha.G2.N2 = coef(USC.G2.N2)[1]
# USCbeta.G2.N2 = coef(USC.G2.N2)[2]

dev.new()
plot(D.G2.N2, SF.G2.N2,
     ylab="Survival fraction", xlab="Dose [Gy]",
     ylim=c(0.001, 1), xlim=c(0,22),
     las=1, pch=16, log="y")

lines(smooth.spline(D.G2.N2, exp(fitted(LQmod.G2.N2))), lty = 1)
lines(smooth.spline(D.G2.N2, exp(fitted(MTSHmod.G2.N2))), lty = 2)
lines(sort(D.G2.N2), exp(USC.G2.N2)[order(D.G2.N2)], lty = 3)
#lines(smooth.spline(Dose.Gy.[Cellcycle == "G2" & N2oxia == "Yes"], (fitted(GLMmod)[Cellcycle == "G2" & N2oxia == "Yes"]) ), lty = 4)
legend("topright", legend=c("G2 (15 hours)", "Extremely N2oxic", "LQ", "MTSH", "USC"), 
       lty=c(0,0,1,2,3), pch=c(NA,16,NA,NA,NA), bty="o")

write( "G2 (15 hours)", file = "output.txt", append = T )
write( c("alpha (N2)", "beta (N2)", "DT (N2)"), file = "output.txt", ncolumns = 3, sep = "\t\t", append = T )
write( c( (-1)*coef(summary(LQmod.G2.N2))[, 1], DT.G2.N2), file = "output.txt", ncolumns = 3, sep = "\t\t", append = T )
write( c( coef(summary(LQmod.G2.N2))[, 2], "comming"), file = "output.txt", ncolumns = 3, sep = "\t\t", append = T )
write("--------------------------------------------------------------------------------", file = "output.txt", append = T )

#summary(LQmod.G2.N2)
#confint(LQmod.G2.N2, conf.level=0.95)
#summary(LQmod.G2.N2)$adj.r.squared
#summary(MTSHmod.G2.N2)


# Diagnostics
LQresid.G2.N2 = resid(LQmod.G2.N2)
MTSHresid.G2.N2 = resid(MTSHmod.G2.N2)
USCresid.G2.N2 = log(SF.G2.N2) - USC.G2.N2

dev.new()
plot(exp(fitted(LQmod.G2.N2)), LQresid.G2.N2,
     ylab = "Residuals", xlab = "Estimated survival fraction (fitted values)",
     ylim = c(-1,1), xlim = c(0, 1), las = 1, pch = 1)
points(exp(fitted(MTSHmod.G2.N2)), MTSHresid.G2.N2, pch = 2)
points(exp(USC.G2.N2), USCresid.G2.N2, pch = 3)
abline(0, 0, lty=2)

LOESS.LQresid.G2.N2 = loess(LQresid.G2.N2 ~ fitted(LQmod.G2.N2) )
LOESS.MTSHresid.G2.N2 = loess(MTSHresid.G2.N2 ~ fitted(MTSHmod.G2.N2) )
LOESS.USCresid.G2.N2 = loess(USCresid.G2.N2 ~ USC.G2.N2)

lines(sort(exp(fitted(LQmod.G2.N2))), fitted(LOESS.LQresid.G2.N2)[order(exp(fitted(LQmod.G2.N2)))], lty = 1)
lines(sort(exp(fitted(MTSHmod.G2.N2))), fitted(LOESS.MTSHresid.G2.N2)[order(exp(fitted(MTSHmod.G2.N2)))], lty = 2)
lines(sort(exp(fitted(MTSHmod.G2.N2))), fitted(LOESS.USCresid.G2.N2)[order(exp(USC.G2.N2))], lty = 3)

legend("topright", legend = c("G2 (15 hours); extremely N2oxic", "LQ", "MTSH", "USC"), 
       lty = c(NA, 1, 2, 3), pch = c(NA, 1, 2, 3), bty = "o")

# -------------------------------------------------------------------------------------------------------------------------- #
# GLM
GLM.FullMod = lm(log(ï..SF) ~ -1 + Dose..Gy. + I(Dose..Gy.^2) + Dose..Gy.*Cellcycle + Dose..Gy.*N2oxia  + Cellcycle:N2oxia )
GLM.ReducedMod = glm(ï..SF ~ -1 + Dose..Gy. + I(Dose..Gy.^2) + Dose..Gy.*Cellcycle + Dose..Gy.*N2oxia + Cellcycle:N2oxia, family = quasipoisson(link = "log")  )
GLMresid = rstandard(GLMmod)
LOESS.GLMresid = loess(GLMresid ~ fitted(GLMmod))

dev.new()
plot(fitted(GLMmod), GLMresid,
     ylab = "Standardized residuals", xlab = "Estimated survival fraction (fitted values)",
     las = 1, pch = 1)
abline(0, 0, lty = 2)
lines(sort(fitted(GLMmod)), fitted(LOESS.GLMresid)[order(fitted(GLMmod))], lty=1)
legend("topright", legend = c("LOESS curve"), lty = c(1), pch = c(1), bty = "o")

ANOVA = anova(GLMmod)
anova(GLMmod, test = "Chisq")
summary(ANOVA)
summary(GLMmod)
confint(GLMmod) # 95% CI for the coef

dev.new()
LOESS.GLMresid.dose = loess(GLMresid ~ Dose.Gy.)
plot(Dose.Gy., GLMresid, 
     ylab = "Standardized residuals", xlab = "Dose [Gy]",
     las = 1, pch = 1)
abline(0, 0, lty = 2)
lines(sort(Dose.Gy.), fitted(LOESS.GLMresid.dose)[order(Dose.Gy.)] , lty=1)
legend("topright", legend = c("Residual versus predictor plot", "LOESS curve"), lty = c(NA, 1), pch = c(NA, 1), bty = "o")

dev.new()
boxplot(GLMresid ~ N2oxia*Cellcycle, ylab = "Standardized residuals", las = 2, ylim = c(-2.5, 2.5) ,col = c(2,5))
abline(0, 0, lty = 2)
legend("topright", legend = c("Residual versus predictor plot"), bty = "o")

GLMmod1 = glm(SF[Cellcycle=="G1" & N2oxia=="Yes"] ~ -1 + Dose.Gy.[Cellcycle=="G1" & N2oxia=="Yes"] + I(Dose.Gy.[Cellcycle=="G1" & N2oxia=="Yes"]^2), family = quasipoisson(link = "log"))

plot(Dose.Gy.[Cellcycle == "G1" & N2oxia == "Yes"], SF[Cellcycle == "G1" & N2oxia == "Yes"], 
     ylab = "Survival fraction", xlab = "Dose [Gy]", 
     ylim = c(5e-05, 5), xlim = c(0, 30),
     las = 1, pch = 1, log = "y")
lines(smooth.spline(Dose.Gy.[Cellcycle=="G1" & N2oxia=="Yes"], fitted(GLMmod1)), lty=1)
lines(smooth.spline(Dose.Gy.[Cellcycle=="G1" & N2oxia=="Yes"], fitted(GLMmod)[Cellcycle=="G1" & N2oxia=="Yes"] ), lty=2)

# ---------------------------------------------------------------------------------------------- #
# Test 
rho.G1 = 0.481
rho.S = 0.368
rho.M = 0.151
HF = c(0.01, 0.05, 0.1, 0.2)

D.test = c(Dlist.N2[[2]], Dlist.O2[[1]], Dlist.O2[[2]], Dlist.O2[[3]])
SF.test = c(HF*SFlist.N2[[2]], rho.M*(1-HF)*SFlist.O2[[1]], rho.G1*(1-HF)*SFlist.O2[[2]], rho.S*(1-HF)*SFlist.O2[[3]])
#D.test = c(Dlist.N2[[4]], Dlist.O2[[4]])
#SF.test = c(f*SFlist.N2[[4]], (1-f)*SFlist.O2[[4]])

LQmod.test = lm(log(SF.test) ~ -1 + D.test + I(D.test^2))
MTSHmod.test = nls(log(SF.test) ~ log(1 - (1 - exp(-v*D.test) )^m ), start = list(v=1, m=1))

#LOESS.test = loess(log(SF.test) ~ D.test)

# USC
alpha.test = (-1)*coef(LQmod.test)[1]
beta.test = (-1)*coef(LQmod.test)[2]
D0.test = 1/coef(MTSHmod.test)[1]
Dq.test = log(coef(MTSHmod.test)[2])*D0.test
DT.test = (2*Dq.test)/(1 - alpha.test*D0.test)

USC.test = USCfunc(alpha.test, beta.test, D.test, DT.test)

# LQ.test = computeFit(D.test, SF.test)
# USC.test = computeFit(D.test, SF.test, DT.test)

dev.new()
plot(D.test, SF.test,
     ylab = "Survival fraction", xlab = "Dose [Gy]",
     las = 1, pch = 1, log = "y")
lines(sort(D.test), exp(fitted(LQmod.test))[order(D.test)],  lty = 1)

lines(sort(D.test), exp(USC.test)[order(D.test)], lty = 3)
#lines(sort(D.test), exp(fitted(LOESS.test))[order(D.test)], lty = 1)

legend("topright", legend=c("LQ", "USC"), lty = c(1, 3), bty = "o")


D.test.N2 = c(D.M.N2, D.G1.N2, D.S.N2, D.G2.N2)
SF.test.N2 = c((rho.G2.M/2)*SF.M.N2, rho.G1*SF.G1.N2, rho.S*SF.S.N2, (rho.G2.M/2)*SF.G2.N2)
mod.test.N2 = lm(log(SF.test.N2) ~ -1 + D.test.N2 + I(D.test.N2^2))
LOESS.N2 = loess(log(SF.test.N2) ~ D.test.N2)

dev.new()
plot(D.test.N2,  exp(fitted(mod.test.N2)),
     ylab = "Survival fraction", xlab = "Dose [Gy]",
     ylim = c(0.0001, 1), xlim = c(0, 27),
     log = "y", las = 1, pch = 1)
points(D.Asynch.N2, exp(fitted(LQmod.Asynch.N2)), pch=2)

#lines(smooth.spline(D.test, exp(fitted(mod.test))), lty=1)
lines(sort(D.test.N2), exp(fitted(LOESS.N2))[order(D.test.N2)], lty=1)
lines(smooth.spline(D.Asynch.N2, exp(fitted(LQmod.Asynch.N2))), lty=2)
lines(smooth.spline(D.test.N2, exp(fitted(mod.test.N2))), lty=3)

legend("topright", legend=c("Synchronous; LOESS", "Asynchronous; LQ"), 
       lty=c(1,2), pch=c(1,2), bty="o")


D.test.O2 = c(D.M.O2, D.G1.O2, D.S.O2)
SF.test.O2 = c(rho.G2.M*SF.M.O2, rho.G1*SF.G1.O2, rho.S*SF.S.O2)
SFhat.O2 = c(rho.G2.M*exp(fitted(LQmod.M.O2)), rho.G1*exp(fitted(LQmod.G1.O2)), rho.S*exp(fitted(LQmod.S.O2)) )

mod.test.O2 = lm(log(SFhat.O2) ~ -1 + D.test.O2 + I(D.test.O2^2))
LOESS.O2 = loess(log(SFhat.O2) ~ D.test.O2)

dev.new()
plot(D.test.O2,  SFhat.O2,
     ylab = "Survival fraction", xlab = "Dose [Gy]",
     ylim = c(0.0001, 1), xlim = c(0, 10),
     log = "y", las = 1, pch = 1)
points(D.Asynch.O2, SF.Asynch.O2, pch=2)

#lines(smooth.spline(D.test, exp(fitted(mod.test))), lty=1)
lines(sort(D.test.O2), exp(fitted(LOESS.O2))[order(D.test.O2)], lty=1)
lines(smooth.spline(D.Asynch.O2, exp(fitted(LQmod.Asynch.O2))), lty=2)
lines(smooth.spline(D.test.O2, exp(fitted(mod.test.O2))), lty=3)

# ---------------------------------------------------------------------------------------------- #
alpha.beta.normal = 3
alpha.normal = 0.206
beta.normal = alpha.normal/alpha.beta.normal
DT.normal = 5.8

alpha.beta.tumour = 10
alpha.tumour = 0.3446
beta.tumour = alpha.tumour/alpha.beta.tumour
DT.tumour = 6.61

nref = 60
dref = 1.8 
Dref = nref*dref # EQD2 
n = 3
d = rep(2:30, times = 1)

set.seed(1)
Diso.LQ.normal = Dref*(dref + alpha.beta.normal )/( d + alpha.beta.normal )
Diso.LQ.tumour = Dref*(dref + alpha.beta.tumour )/( d + alpha.beta.tumour )
Diso.USC.normal = Dref*( alpha.normal + (beta.normal/dref)*I(f(dref, DT.normal)) )/( alpha.normal + (beta.normal/d)*I(f(d, DT.normal)) )
Diso.USC.tumour = Dref*( alpha.tumour + (beta.tumour/dref)*I(f(dref, DT.tumour)) )/( alpha.tumour + (beta.tumour/d)*I(f(d, DT.tumour)) ) 

plot(d, Diso.LQ.normal,
     ylab = "Total dose [Gy]", xlab = "Dose/fraction [Gy]", xlim = rev(range(d)),
     las = 1, log = "x", type = "l", lty = 1, col = "Red")
points(d, Diso.LQ.tumour, type = "l", lty = 1, col = "Blue")
points(d, Diso.USC.normal, type = "l", lty = 2, col = "Red")
points(d, Diso.USC.tumour, type = "l", lty = 2, col = "Blue")
legend("topleft", legend=c("LQ normal", "USC normal", "LQ tumour", "USC tumour"), 
       lty = c(1, 2, 1, 2), col = c("Red", "Red", "Blue", "Blue"), bty = "o")

# ------------------------------------------------------------------------------------------------------------------------ #
alpha.test = 0.33
beta.test = 0.0384
DT.test = 6.2
D0.test = 1.25
Dq.test = 1.8
HF = 0.2

D.test = seq(from = 0, to = 35, length.out = 50)

LQ.test.N2 = LQfunc(alpha.test, beta.test, D.test, 3)
USC.test.N2 = USCfunc(alpha.test, beta.test, D.test, DT.test, 3)
LQ.test.O2 = LQfunc(alpha.test, beta.test, D.test, 1)
USC.test.O2 = USCfunc(alpha.test, beta.test, D.test, DT.test, 1)
SFmix.LQ = (1 - HF)*exp(LQ.test.O2) + HF*exp(LQ.test.N2)
SFmix.USC = (1 - HF)*exp(USC.test.O2) + HF*exp(USC.test.N2)
 
dev.new()
plot(D.test, exp(LQ.test.N2), 
     ylab = "Survival fraction", xlab = "Dose [Gy]",
     las = 1, pch = 1, log = "y", type = "l", lty = 1, col = "Red")
points(D.test, exp(USC.test.N2), type = "l", lty = 2, col = "Red")
points(D.test, exp(LQ.test.O2), type = "l", lty = 1, col = "Blue")
points(D.test, exp(USC.test.O2), type = "l", lty = 2, col = "Blue")
points(D.test, (SFmix.LQ), type = "l", lty = 1, col = "Green")
points(D.test, (SFmix.USC), type = "l", lty = 2, col = "Green")
legend("topright", legend=c("LQ, extremely N2oxic ", "USC, extremely N2oxic", 
                            "LQ, O2obic", "USC, O2obic", 
                            "LQ, partly N2oxic", "USC, partly N2oxic"), 
       lty = c(1, 2, 1, 2, 1, 2), col = c("Red", "Red", "Blue", "Blue", "Green", "Green"), bty = "o")


# -------------------------------------------------------------------------------------------------------- #
alpha.test = 0.33
beta.test = 0.0384
DT.test = 6.2
D0.test = 1.25
Dq.test = 1.8
HF = 0.2

SF.NSCLC = 10^(c(0, -0.2, -0.8, -1.9, -3.5, -4.6, -5.7, -6.8))
D.NSCLC = c(0, 2, 5, 7.7, 10, 12, 14.1, 16)

#LQ.NSCLC.test = LQfunc(alpha.test, beta.test, D.NSCLC, 1) 
#USC.NSCLC.test = USCfunc(alpha.test, beta.test, D.NSCLC.test, DT.test, 1.) 

LQmod.NSCLC = lm(log(SF.NSCLC) ~ -1 + D.NSCLC + I(D.NSCLC^2))
MTSHmod.NSCLC = nls(log(SF.NSCLC) ~ log(1 - (1 - exp(-v*D.NSCLC) )^m ), start = list(v=1, m=1))

# Model coefficients
alpha.NSCLC = (-1)*coef(LQmod.NSCLC)[1]
beta.NSCLC = (-1)*coef(LQmod.NSCLC)[2]
# v.NSCLC = coef(MTSHmod.NSCLC)[1]
# m.NSCLC = coef(MTSHmod.NSCLC)[2]
# D0.NSCLC = 1/v.NSCLC
# Dq.NSCLC = log(m.NSCLC)*D0.NSCLC
# DT.NSCLC = (2*Dq.NSCLC)/(1 - alpha.NSCLC*D0.NSCLC)

#DT.grid = seq(1, 15, 0.05)
#USCmod.NSCLC = find.DT(D.NSCLC, SF.NSCLC, DT.grid)
#USCmod.NSCLC$c

# Compute USC curve
USC.NSCLC = USCfunc(alpha.NSCLC, beta.NSCLC, D.NSCLC, DT.test, 1)

dev.new()
plot(D.NSCLC, (SF.NSCLC), 
     ylab = "Survival fraction", xlab = "Dose [Gy]",
     las = 1, pch = 1, log = "y")
points(D.NSCLC, exp(fitted(LQmod.NSCLC)), type = "l", lty = 1)
points(D.NSCLC, exp(USC.NSCLC), type = "l", lty = 2)
legend("topright", legend=c("LQ", "USC"), lty = c(1, 2), bty = "o")



