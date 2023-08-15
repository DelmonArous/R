setwd("C:/Users/Delmon/Dropbox/Master/R")
getwd()

## Clear and detach all data
rm(list=ls())
detach(SFdata)

## Read in cell survival observations and attach the object SFdata to R search path
SFdata = read.table(file = "NHIK1922_v2.csv", header = TRUE, sep = ";")
attach(SFdata)

## Define the LQ and USC functions
LQfunc = function(alpha, beta, D, OER) -(alpha/OER)*D - (beta/OER^2)*I(D^2)
f = function(D, DT) I(D^2) - (D > DT)*I((D - DT)^2)
USCfunc = function(alpha, beta, D, DT, OER) -(alpha/OER)*D - (beta/OER^2)*I( f(D, DT) )

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
  ## 95%-prediction interval
  #lwr = fit + alpha * sqrt(model$sig2 + I(se^2))
  #upr = fit + alpha - sqrt(model$sig2 + I(se^2))
  #lwr = fit + tval*s(y,fit,model$n,model$p)*sqrt( 1/(model$n-model$p) + I((x.new - mean(x.new))^2)/Sxx(x.new) )
  #upr = fit - tval*s(y,fit,model$n,model$p)*sqrt( 1/(model$n-model$p) + I((x.new - mean(x.new))^2)/Sxx(x.new) )
  
  ## Return estimated values
  matrix(c(fit, se, lwr, upr), ncol = 4L,
         dimnames = list(NULL, c("fit", "se", "lwr", "upr")))
  
}

Doserange.O2 = c(30, 50)
Doserange.N2 = c(50, 50)
DTrange.O2 = c( c(2, 9), c(5, 25) )
DTrange.N2 = c( c(10.5, 28.5), c(10, 25) )
legend.title = c("Irradiated in vitro (suspension)", "Irradiated in vivo")  
v.vec = c(0.5, 0.2, 0.2, 0.2)
m.vec = c(1, 2, 2, 2)
D0.O2.vec = c(0.644, 0.268)
D0.N2.vec = c(0.223, 0.226)
Dq.O2.vec = c(2.534, 1.589)
Dq.N2.vec = c(1.942, 3.650)

k = 0
for (i in c("in vitro", "in vivo")) {
  
  k = k + 1
  
  if (i == "in vitro") {
    
    ## Extract dose-response observations
    SF.O2 = SFdata[Procedure == i & Hypoxia == "No", 1]
    SF.N2 = SFdata[Procedure == i & Hypoxia == "Yes", 1]
    D.O2 = SFdata[Procedure == i & Hypoxia == "No", 2]
    D.N2 = SFdata[Procedure == i & Hypoxia == "Yes", 2]
  
    ## Linear and non-linear regression of LQ and MTSH, respectively
    LQmod.O2 = lm(log(SF.O2) ~  -1 + D.O2 + I(D.O2^2))
    LQmod.N2 = lm(log(SF.N2) ~  -1 + D.N2 + I(D.N2^2))
    #MTSHmod.O2 = nlsLM(log(SF.O2) ~ log(1 - (1 - exp(-(1/D0)*D.O2))^(exp(Dq/D0)) ), start = list(D0=D0.O2.vec[k] , Dq=Dq.O2.vec[k]))
    #MTSHmod.N2 = nlsLM(log(SF.N2) ~ log(1 - (1 - exp(-(1/D0)*D.N2))^(exp(Dq/D0))  ), start = list(D0=D0.N2.vec[k], Dq=Dq.N2.vec[k]))
    MTSHmod.O2 = nls(log(SF.O2) ~ log(1 - (1 - exp(-v*D.O2) )^m ), start = list(v=v.vec[2*k-1], m=m.vec[2*k-1]))
    MTSHmod.N2 = nls(log(SF.N2) ~ log(1 - (1 - exp(-v*D.N2) )^m ), start = list(v=v.vec[2*k], m=m.vec[2*k]))
  
    ## Model coefficients
    alpha.O2 = (-1)*coef(LQmod.O2)[1]
    beta.O2 = (-1)*coef(LQmod.O2)[2]
    v.O2 = coef(MTSHmod.O2)[1] 
    m.O2 = coef(MTSHmod.O2)[2] 
    D0.O2 = 1/v.O2
    Dq.O2 = D0.O2*log(m.O2)
    DT.O2 = (2*Dq.O2)/(1 - alpha.O2*D0.O2)
    
    alpha.N2 = (-1)*coef(LQmod.N2)[1]
    beta.N2 = (-1)*coef(LQmod.N2)[2]
    v.N2 = coef(MTSHmod.N2)[1]
    m.N2 = coef(MTSHmod.N2)[2]
    D0.N2 = 1/v.N2
    Dq.N2 = D0.N2*log(m.N2)
    DT.N2 = (2*Dq.N2)/(1 - alpha.N2*D0.N2)
    
    ## Compute USC curve
    USC.O2 = USCfunc(alpha.O2, beta.O2, D.O2, DT.O2, 1)
    USC.N2 = USCfunc(alpha.N2, beta.N2, D.N2, DT.N2, 1)
    
    ## Plot dose-response curves for various oxygenation
    dev.new()
    plot(D.O2, SF.O2,
         ylab = "Surviving fraction", xlab = "Dose [Gy]",
         ylim = c(min(SF.O2), 1), xlim = c(0, max(D.N2)),
         las = 1, pch = 1, log = "y", yaxt = "n", 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    points(D.N2, SF.N2, pch = 16)
    ticks = seq(-4, 0, by=1)
    labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
    axis(2, at=c(0.0001, 0.001, 0.01, 0.1, 1), labels=labels, 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    
    D_aer = seq(0, max(D.O2), 0.01)
    D_hyp = seq(0, max(D.N2), 0.01)
    lines(D_aer, exp(LQfunc(alpha.O2, beta.O2, D_aer, 1)), lty = 1)
    lines(D_hyp, exp(LQfunc(alpha.N2, beta.N2, D_hyp, 1)), lty = 1)
    lines(D_aer, (1 - (1 - exp(-v.O2*D_aer))^m.O2), lty = 2)
    lines(D_hyp, (1 - (1 - exp(-v.N2*D_hyp))^m.N2), lty = 2)
    lines(D_aer, exp(USCfunc(alpha.O2, beta.O2, D_aer, DT.O2, 1)), lty = 3)
    lines(D_hyp, exp(USCfunc(alpha.N2, beta.N2, D_hyp, DT.N2, 1)), lty = 3)
    legend("topright", legend=c(legend.title[k], "Aerobic", "Extremely hypoxic", "LQ", "MTSH", "USC"),
           lty = c(0, 0, 0, 1, 2, 3), pch = c(NA, 1, 16, NA, NA, NA), bty = "o")
    
    AIC_MTSH_O2 = 2*2 - 2*logLik(MTSHmod.O2)
    cat("MTSH O2, AIC = ", AIC_MTSH_O2, "\n")

    AIC_LQ_O2 = 2*2 - 2*logLik(LQmod.O2)
    cat("LQ O2, AIC = ", AIC_LQ_O2, "\n")
    
    sig2_USC_O2 = RSS(log(SF.O2), USC.O2)/(length(D.O2) - 4)
    AIC_USC_O2 = length(D.O2)*log(2*pi*sig2_USC_O2) + (length(D.O2) - 4) + 2*4 
    cat("USC O2, AIC = ", AIC_USC_O2, "\n")
    
    AIC_MTSH_N2 = 2*2 - 2*logLik(MTSHmod.N2) 
    cat("MTSH N2, AIC = ", AIC_MTSH_N2, "\n")
    
    AIC_LQ_N2 = 2*2 - 2*logLik(LQmod.N2) 
    cat("LQ N2, AIC = ", AIC_LQ_N2, "\n")
    
    sig2_USC_N2 = RSS(log(SF.N2), USC.N2)/(length(D.N2) - 4)
    AIC_USC_N2 = length(D.N2)*log(2*pi*sig2_USC_N2) + (length(D.N2) - 4) + 2*4 
    cat("USC N2, AIC = ", AIC_USC_N2, "\n")
    
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
         ylim = c(-2, 2) , xlim = c(min(SF.O2) , 1),
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
    LOESS.LQresid.O2 = loess(LQresid.O2 ~ (fitted(LQmod.O2)) )
    LOESS.MTSHresid.O2 = loess(MTSHresid.O2 ~ fitted(MTSHmod.O2) )
    LOESS.USCresid.O2 = loess(USCresid.O2 ~ USC.O2)
    
    lines(sort(exp(fitted(LQmod.O2))), fitted(LOESS.LQresid.O2)[order(exp(fitted(LQmod.O2)))], lty = 1)
    lines(sort(exp(fitted(MTSHmod.O2))), fitted(LOESS.MTSHresid.O2)[order(exp(fitted(MTSHmod.O2)))], lty = 2)
    lines(sort(exp(USC.O2)), fitted(LOESS.USCresid.O2)[order(exp(USC.O2))], lty = 3)
    legend("topright", legend = c(paste(legend.title[k], "- aerobic"), "LQ", "MTSH", "USC"),
           lty = c(NA, 1, 2, 3), pch = c(NA, 1, 2, 3), bty = "o")
    
    ## Extremely-hypoxic residual plot
    dev.new()
    plot(exp(fitted(LQmod.N2)), LQresid.N2,
         ylab = "Residuals", xlab = "Estimated surviving fraction (fitted values)",
         ylim = c(-2, 2) , xlim = c(min(SF.N2), 1),
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
    legend("topright", legend = c(paste(legend.title[k], "- extremely hypoxic"), "LQ", "MTSH", "USC"),
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
    plot(DTspan.O2, RSS.O2, 
         ylab = "RSS", xlab = "DT grid [Gy]", 
         ylim = c(min(min(RSS.O2), min(RSS.N2)), max(max(RSS.O2), max(RSS.N2))), 
         xlim = c(min(min(DTspan.O2), min(DTspan.N2)), max(max(DTspan.O2), max(DTspan.N2))),
         type = "o", pch = c(1, rep(NA, 44)), 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    lines(DTspan.N2, RSS.N2, type = "o", pch = c(16, rep(NA, 44)))
    legend("topright", legend=c(legend.title[k], "Aerobic", "Extremely hypoxic"),
           pch = c(NA, 1, 16), lty = c(NA, 1, 1), bty = "o")
    
    ### Extrapolation and prediction interval
    
    ## Sample dose points on the same range as the original aerobic experimental data
    set.seed(k)
    Dnew.O2 = runif(150, min = 0, max = max(D.O2))
    
    ## Generate observation points with Gaussian noise around the best-fit curve
    ## predicted by the LQ and USC model in the dose range of Dnew.O2
    set.seed(k)
    SFnew.LQO2 = LQfunc(alpha.O2, beta.O2, Dnew.O2, 1) + rnorm(150, 0, sd(SF.O2))
    set.seed(k)
    SFnew.USCO2 = USCfunc(alpha.O2, beta.O2, Dnew.O2, DT.O2, 1) + rnorm(150, 0, sd(SF.O2))
    
    ## Sample dose points on the same range as the original extremely-hypoxic experimental data
    set.seed(k)
    Dnew.N2 = runif(150, min = 0, max = max(D.N2))
    
    ## Generate observation points with Gaussian noise around the best-fit curve
    ## predicted by the LQ and USC model in the dose range of Dnew.N2
    set.seed(k)
    SFnew.LQN2 = LQfunc(alpha.N2, beta.N2, Dnew.N2, 1) + rnorm(150, 0, sd(SF.N2))
    set.seed(k)
    SFnew.USCN2 = USCfunc(alpha.N2, beta.N2, Dnew.N2, DT.N2, 1) + rnorm(150, 0, sd(SF.N2))
    
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
            log = "y", las = 1, col = c(1,4,4), lty = c(1,2,2), type = "l", lwd = 1, yaxt = "n",
            cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    ticks = seq(-7, 0, by=1)
    labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
    axis(2, at=c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1), labels=labels, 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    matlines( Dtest.O2, exp(p.LQO2[,-2]), col = c(1,4,4), lty = c(1,2,2), type = "l", lwd = 1 )
    matlines( Dtest.N2, exp(p.USCN2[,-2]), col = c(1,2,2), lty = c(3,2,2), type = "l", lwd = 1 )
    matlines( Dtest.O2, exp(p.USCO2[,-2]), col = c(1,2,2), lty = c(3,2,2), type = "l", lwd = 1 )
    points(Dnew.O2, exp(SFnew.LQO2), pch = 1, cex = 0.5)
    points(Dnew.N2, exp(SFnew.LQN2), pch = 16, cex = 0.5)
    points(Dnew.O2, exp(SFnew.USCO2), pch = 1, cex = 0.5)
    points(Dnew.N2, exp(SFnew.USCN2), pch = 16, cex = 0.5)
    legend("topright", legend=c(legend.title[k], "Aerobic", "Extremely hypoxic", "LQ", "USC", "LQ 95% CI", "USC 95% CI"),
           lty = c(NA, 0, 0, 1, 3, 2, 2), pch = c(NA, 1, 16, NA, NA, NA, NA), col = c(NA, 1, 1, 1, 1, 4, 2), lwd = 2, bty = "o")
    
    ## Print model parameters and their corresponding SD
    SDalpha.O2 = coef(summary(LQmod.O2))[1, 2]
    SDalpha.N2 = coef(summary(LQmod.N2))[1, 2]
    SDbeta.O2 = coef(summary(LQmod.O2))[2, 2]
    SDbeta.N2 = coef(summary(LQmod.N2))[2, 2]
    SDv.O2 = coef(summary(MTSHmod.O2))[1, 2]
    SDv.N2 = coef(summary(MTSHmod.N2))[1, 2]
    SDm.O2 = coef(summary(MTSHmod.O2))[2, 2]
    SDm.N2 = coef(summary(MTSHmod.N2))[2, 2]
    SDD0.O2 = D0.O2^2*SDv.O2
    SDDq.O2 = sqrt( (SDD0.O2*log(m.O2))^2 + (SDm.O2*D0.O2/m.O2)^2 )
    SDD0.N2 = D0.N2^2*SDv.N2
    SDDq.N2 = sqrt( (SDD0.N2*log(m.N2))^2 + (SDm.N2*D0.N2/m.N2)^2 )
    SDDT.O2 = (DT.O2/Dq.O2)*sqrt( SDDq.O2^2 + (D0.O2*DT.O2*SDalpha.O2/2)^2 + (alpha.O2*DT.O2*SDD0.O2/2)^2 )
    SDDT.N2 = (DT.N2/Dq.N2)*sqrt( SDDq.N2^2 + (D0.N2*DT.N2*SDalpha.N2/2)^2 + (alpha.N2*DT.N2*SDD0.N2/2)^2 )
    
    #SDD0.O2 = coef(summary(MTSHmod.O2))[1, 2]
    #SDD0.N2 = coef(summary(MTSHmod.N2))[1, 2]
    #SDDq.O2 = coef(summary(MTSHmod.O2))[2, 2]
    #SDDq.N2 = coef(summary(MTSHmod.N2))[2, 2]
    
    alpha.ratio = alpha.O2/alpha.N2
    beta.ratio = sqrt(beta.O2/beta.N2)
    alphabeta.ratio.O2 = alpha.O2/beta.O2
    alphabeta.ratio.N2 = alpha.N2/beta.N2
    SDalpha.ratio = alpha.ratio*sqrt( (SDalpha.O2/alpha.O2)^2 + (SDalpha.N2/alpha.N2)^2 )
    SDbeta.ratio = 0.5*beta.ratio*sqrt( (SDbeta.O2/beta.O2)^2 + (SDbeta.N2/beta.N2)^2 )
    SDalphabeta.ratio.O2 = alphabeta.ratio.O2*sqrt( (SDalpha.O2/alpha.O2)^2 + (SDbeta.O2/beta.O2)^2 )
    SDalphabeta.ratio.N2 = alphabeta.ratio.N2*sqrt( (SDalpha.N2/alpha.N2)^2 + (SDbeta.N2/beta.N2)^2 )
    
    write(legend.title[k], file = "NHIK1922.txt", append = T)
    write(c("alpha (O2)", "beta (O2)", "alpha (N2)", "beta (N2)",
            "D0 (O2)", "Dq (O2)", "D0 (N2)", "Dq (N2)",
            "alpha_O2/alpha_N2", "sqrt(beta_O2/beta_N2)",
            "DT (O2)", "DT (N2)", "alpha/beta (O2)", "alpha/beta (N2)"),
          file = "NHIK1922.txt", ncolumns = 14, sep = "\t\t\t", append = T)
    write(c( (-1)*coef(summary(LQmod.O2))[, 1], (-1)*coef(summary(LQmod.N2))[, 1],
             D0.O2, Dq.O2, D0.N2, Dq.N2,
             alpha.ratio, beta.ratio, DT.O2, DT.N2, alphabeta.ratio.O2, alphabeta.ratio.N2), 
          file = "NHIK1922.txt", ncolumns = 14, sep = "\t\t\t", append = T)
    write(c( coef(summary(LQmod.O2))[, 2], coef(summary(LQmod.N2))[, 2],
             SDD0.O2, SDDq.O2, SDD0.N2, SDDq.N2,
             SDalpha.ratio, SDbeta.ratio, SDDT.O2, SDDT.N2, SDalphabeta.ratio.O2, SDalphabeta.ratio.N2),
          file = "NHIK1922.txt", ncolumns = 14, sep = "\t\t\t", append = T)
    write( c( "Adjusted R-sqaured (O2): ", summary(LQmod.O2)$adj.r.squared, R.squared.adj(log(SF.O2), USC.O2, length(D.O2), 4) ),
           file = "NHIK1922.txt", ncolumns = 3, sep = "\t", append = T)
    write( c( "Adjusted R-squared (N2): ", summary(LQmod.N2)$adj.r.squared, R.squared.adj(log(SF.N2), USC.N2, length(D.N2), 4) ),
           file = "NHIK1922.txt", ncolumns = 3, sep = "\t", append = T)
    write( c( "RSS (O2): ", RSS(log(SF.O2), fitted(MTSHmod.O2)), RSS(log(SF.O2), fitted(LQmod.O2)), RSS(log(SF.O2), USC.O2)),
           file = "NHIK1922.txt", ncolumns = 4, sep = "\t\t\t", append = T)
    write( c( "RSS (N2): ", RSS(log(SF.N2), fitted(MTSHmod.N2)), RSS(log(SF.N2), fitted(LQmod.N2)), RSS(log(SF.N2), USC.N2)),
           file = "NHIK1922.txt", ncolumns = 4, sep = "\t\t\t", append = T)
    write("-------------------------------------------------------------------------------------------------------------",
          file = "NHIK1922.txt", append = T)
    
  } else {
    
    ## Extract dose-response observations
    SF.N2 = SFdata[Procedure == i & Hypoxia == "Yes", 1]
    D.N2 = SFdata[Procedure == i & Hypoxia == "Yes", 2]
    
    ## Linear and non-linear regression of LQ and MTSH, respectively
    LQmod.N2 = lm(log(SF.N2) ~  -1 + D.N2 + I(D.N2^2))
    #MTSHmod.N2 = nlsLM(log(SF.N2) ~ log(1 - (1 - exp(-(1/D0)*D.N2))^(exp(Dq/D0))  ), start = list(D0=D0.N2.vec[k], Dq=Dq.N2.vec[k]))
    MTSHmod.N2 = nls(log(SF.N2) ~ log(1 - (1 - exp(-v*D.N2) )^m ), start = list(v=v.vec[2*k], m=m.vec[2*k]))
    
    ## Model coefficients
    alpha.N2 = (-1)*coef(LQmod.N2)[1]
    beta.N2 = (-1)*coef(LQmod.N2)[2]
    v.N2 = coef(MTSHmod.N2)[1]
    m.N2 = coef(MTSHmod.N2)[2]
    D0.N2 = 1/v.N2
    Dq.N2 = D0.N2*log(m.N2)
    DT.N2 = (2*Dq.N2)/(1 - alpha.N2*D0.N2)
    
    ## Compute USC curve
    USC.N2 = USCfunc(alpha.N2, beta.N2, D.N2, DT.N2, 1)
    
    ## Plot dose-response curves for various oxygenation
    dev.new()
    plot(D.N2, SF.N2,
         ylab = "Surviving fraction", xlab = "Dose [Gy]",
         ylim = c(min(SF.N2), 1), xlim = c(0, max(D.N2)),
         las = 1, pch = 1, log = "y", yaxt = "n", 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    points(D.N2, SF.N2, pch = 16)
    ticks = seq(-4, 0, by=1)
    labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
    axis(2, at=c(0.0001, 0.001, 0.01, 0.1, 1), labels=labels, 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    
    D_hyp = seq(0, max(D.N2), 0.01)
    lines(D_hyp, exp(LQfunc(alpha.N2, beta.N2, D_hyp, 1)), lty = 1)
    lines(D_hyp, (1 - (1 - exp(-v.N2*D_hyp))^m.N2), lty = 2)
    lines(D_hyp, exp(USCfunc(alpha.N2, beta.N2, D_hyp, DT.N2, 1)), lty = 3)
    legend("topright", legend=c(legend.title[k], "Extremely hypoxic", "LQ", "MTSH", "USC"),
           lty = c(0, 0, 1, 2, 3), pch = c(NA, 16, NA, NA, NA), bty = "o")
    
    ## Compute model residuals
    LQresid.N2 = resid(LQmod.N2)
    MTSHresid.N2 = resid(MTSHmod.N2)
    USCresid.N2 = log(SF.N2) - USC.N2 
    
    ## Extremely-hypoxic residual plot
    dev.new()
    plot(exp(fitted(LQmod.N2)), LQresid.N2,
         ylab = "Residuals", xlab = "Estimated surviving fraction (fitted values)",
         ylim = c(-2, 2) , xlim = c(min(SF.N2), 1),
         log = "x", las = 1, pch = 1, xaxt = "n", 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    points(exp(fitted(MTSHmod.N2)), MTSHresid.N2, pch = 2)
    points(exp(USC.N2), USCresid.N2, pch = 3)
    abline(0, 0, lty = 2)
    ticks = seq(-4, 0, by=1)
    labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
    axis(1, at=c(0.0001, 0.001, 0.01, 0.1, 1), labels=labels, 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    
    AIC_MTSH_N2 = 2*2 - 2*logLik(MTSHmod.N2)
    cat("MTSH N2, AIC = ", AIC_MTSH_N2, "\n")
    
    AIC_LQ_N2 = 2*2 - 2*logLik(LQmod.N2) 
    cat("LQ N2, AIC = ", AIC_LQ_N2, "\n")
    
    sig2_USC_N2 = RSS(log(SF.N2), USC.N2)/(length(D.N2) - 4)
    AIC_USC_N2 = length(D.N2)*log(2*pi*sig2_USC_N2) + (length(D.N2) - 4) + 2*4 
    cat("USC N2, AIC = ", AIC_USC_N2, "\n")
    
    ## Add trend curves in the extremely-hypoxic residual plot for all models
    LOESS.LQresid.N2 = loess(LQresid.N2 ~ (fitted(LQmod.N2)) )
    LOESS.MTSHresid.N2 = loess(MTSHresid.N2 ~ fitted(MTSHmod.N2) )
    LOESS.USCresid.N2 = loess(USCresid.N2 ~ USC.N2)
    
    lines(sort(exp(fitted(LQmod.N2))), fitted(LOESS.LQresid.N2)[order(exp(fitted(LQmod.N2)))], lty = 1)
    lines(sort(exp(fitted(MTSHmod.N2))), fitted(LOESS.MTSHresid.N2)[order(exp(fitted(MTSHmod.N2)))], lty = 2)
    lines(sort(exp(USC.N2)), fitted(LOESS.USCresid.N2)[order(exp(USC.N2))], lty = 3)
    legend("topright", legend = c(paste(legend.title[k], "- extremely hypoxic"), "LQ", "MTSH", "USC"),
           lty = c(NA, 1, 2, 3), pch = c(NA, 1, 2, 3), bty = "o")
    
    DTspan.N2 = seq(DTrange.N2[2*k-1], DTrange.N2[2*k], 0.01)
    D = D.N2
    SF = log(SF.N2)
    
    ## Regression for each value of DT in span 
    lst.N2 = lapply(DTspan.N2, computeFit, D = D, SF = SF)
    
    ## Compute RSS for each value of DT
    RSS.N2 = sapply(lst.N2, "[[", "RSS")
    
    ## Plot RSS vs grid of DT
    dev.new()
    plot(DTspan.N2, RSS.N2, 
         ylab = "RSS", xlab = "DT grid [Gy]", 
         ylim = c(min(RSS.N2), max(RSS.N2)), 
         xlim = c(min(DTspan.N2), max(DTspan.N2)),
         type = "o", pch = c(16, rep(NA, 44)), 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    legend("topright", legend=c(legend.title[k], "Extremely hypoxic"),
           pch = c(NA, 16), lty = c(NA, 1), bty = "o")
    
    ### Extrapolation and prediction interval
    
    ## Sample dose points on the same range as the original extremely-hypoxic experimental data
    set.seed(k)
    Dnew.N2 = runif(150, min = 0, max = max(D.N2))
    
    ## Generate observation points with Gaussian noise around the best-fit curve
    ## predicted by the LQ and USC model in the dose range of Dnew.N2
    set.seed(k)
    SFnew.LQN2 = LQfunc(alpha.N2, beta.N2, Dnew.N2, 1) + rnorm(150, 0, sd(SF.N2))
    set.seed(k)
    SFnew.USCN2 = USCfunc(alpha.N2, beta.N2, Dnew.N2, DT.N2, 1) + rnorm(150, 0, sd(SF.N2))
    
    ## Solve least squares and find the best fit for both LQ and USC
    newLQfit.N2 = computeFit(Dnew.N2, SFnew.LQN2)
    newUSCfit.N2 = computeFit(Dnew.N2, SFnew.USCN2, DT.N2)
    
    ## Compute 95% CI for the extended dose range and plot the predictions
    Dtest.N2 = seq(from = 0, to = Doserange.N2[k], by = 0.5)
    p.LQN2 = CI(newLQfit.N2, Dtest.N2)
    p.USCN2 = CI(newUSCfit.N2, Dtest.N2)
    
    dev.new()
    matplot(Dtest.N2, exp(p.LQN2[,-2]),
            ylab = "Surviving fraction", xlab = "Dose [Gy]",
            xlim = c(0, Doserange.N2[k]), ylim = c(0.0000001, 1),
            log = "y", las = 1, col = c(1,4,4), lty = c(1,2,2), type = "l", lwd = 1, yaxt = "n", 
            cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    ticks = seq(-7, 0, by=1)
    labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
    axis(2, at=c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1), labels=labels, 
         cex.lab = 1.5, cex.axis = 1.5, cex.sub = 1.5)
    matlines( Dtest.N2, exp(p.USCN2[,-2]), col = c(1,2,2), lty = c(3,2,2), type = "l", lwd = 1 )
    points(Dnew.N2, exp(SFnew.LQN2), pch = 16, cex = 0.5)
    points(Dnew.N2, exp(SFnew.USCN2), pch = 16, cex = 0.5)
    legend("topright", legend=c(legend.title[k], "Extremely hypoxic", "LQ", "USC", "LQ 95% CI", "USC 95% CI"),
           lty = c(NA, 0, 1, 3, 2, 2), pch = c(NA, 16, NA, NA, NA, NA), col = c(NA, 1, 1, 1, 4, 2), lwd = 2, bty = "o")
    
    ## Print model parameters and their corresponding SD
    SDalpha.N2 = coef(summary(LQmod.N2))[1, 2]
    SDbeta.N2 = coef(summary(LQmod.N2))[2, 2]
    SDv.N2 = coef(summary(MTSHmod.N2))[1, 2]
    SDm.N2 = coef(summary(MTSHmod.N2))[2, 2]
    SDD0.N2 = D0.N2^2*SDv.N2
    SDDq.N2 = sqrt( (SDD0.N2*log(m.N2))^2 + (SDm.N2*D0.N2/m.N2)^2 )
    SDDT.N2 = (DT.N2/Dq.N2)*sqrt( SDDq.N2^2 + (D0.N2*DT.N2*SDalpha.N2/2)^2 + (alpha.N2*DT.N2*SDD0.N2/2)^2 )
    
    #SDD0.N2 = coef(summary(MTSHmod.N2))[1, 2]
    #SDDq.N2 = coef(summary(MTSHmod.N2))[2, 2]
    
    alphabeta.ratio.N2 = alpha.N2/beta.N2
    SDalphabeta.ratio.N2 = alphabeta.ratio.N2*sqrt( (SDalpha.N2/alpha.N2)^2 + (SDbeta.N2/beta.N2)^2 )
    
    write(legend.title[k], file = "NHIK1922.txt", append = T)
    write(c("alpha (N2)", "beta (N2)", "D0 (N2)", "Dq (N2)", "DT (N2)", "alpha/beta (N2)"),
          file = "NHIK1922.txt", ncolumns = 6, sep = "\t\t\t", append = T)
    write(c( (-1)*coef(summary(LQmod.N2))[, 1], D0.N2, Dq.N2, DT.N2, alphabeta.ratio.N2), 
          file = "NHIK1922.txt", ncolumns = 6, sep = "\t\t\t", append = T)
    write(c( coef(summary(LQmod.N2))[, 2], SDD0.N2, SDDq.N2, SDDT.N2, SDalphabeta.ratio.N2),
          file = "NHIK1922.txt", ncolumns = 6, sep = "\t\t\t", append = T)
    write( c( "Adjusted R-squared (N2): ", summary(LQmod.N2)$adj.r.squared, R.squared.adj(log(SF.N2), USC.N2, length(D.N2), 4) ),
           file = "NHIK1922.txt", ncolumns = 3, sep = "\t", append = T)
    write( c( "RSS (N2): ", RSS(log(SF.N2), fitted(MTSHmod.N2)), RSS(log(SF.N2), fitted(LQmod.N2)), RSS(log(SF.N2), USC.N2)),
           file = "NHIK1922.txt", ncolumns = 4, sep = "\t\t\t", append = T)
    write("-------------------------------------------------------------------------------------------------------------",
          file = "NHIK1922.txt", append = T)
      
  }
  
}
  
# --------------------------------------------------------------------------------------------------------------------------- #

SF.test = SFdata[Procedure == "in vivo" & Hypoxia == "No" , 1]
D.test = SFdata[Procedure == "in vivo" & Hypoxia == "No" , 2]

model = nls(SF.test ~ HF*exp(-vecalpha.N2[2]*D.test - vecbeta.N2[2]*I(D.test^2)) + 
                (1 - HF)*exp(-a*D.test - b*I(D.test^2)),
            start = list(HF = 0.3, a = 0.05, b = 0.008), 
            control = nls.lm.control(maxiter = 1024))


# -------------------------------------------------------------------------------------------------- #

SF.vivo.O2 = SFdata[Procedure == "in vivo" & Hypoxia == "No", 1]
SF.vivo.N2 = SFdata[Procedure == "in vivo" & Hypoxia == "Yes", 1]
D.vivo.O2 = SFdata[Procedure == "in vivo" & Hypoxia == "No", 2]/100
D.vivo.N2 = SFdata[Procedure == "in vivo" & Hypoxia == "Yes", 2]/100

## Aerobic model
LQmod.vivo.O2 = HFopt*exp(LQfunc(alpha.vivo.N2, beta.vivo.N2, D.vivo.O2, 1)) + 
                  (1 - HFopt)*exp(LQfunc(alpha.vivo.O2, beta.vivo.O2, D.vivo.O2, 1))
USCmod.vivo.O2 = HFopt*exp(USCfunc(alpha.vivo.N2, beta.vivo.N2, D.vivo.O2, 18.57297, 1)) + 
                  (1 - HFopt)*exp(USCfunc(alpha.vivo.O2, beta.vivo.O2, D.vivo.O2, 18.57297, 1))

## Hypoxic model
LQmod.vivo.N2 = lm(log(SF.vivo.N2) ~ -1 + D.vivo.N2 + I(D.vivo.N2^2))
MTSHmod.vivo.N2 = nls(log(SF.vivo.N2) ~ log(1 - (1 - exp(-v*D.vivo.N2) )^m ), start = list(v=0.2, m=2))

alpha.vivo.N2 = (-1)*coef(LQmod.vivo.N2)[1]
beta.vivo.N2 = (-1)*coef(LQmod.vivo.N2)[2]
v.vivo.N2 = coef(MTSHmod.vivo.N2)[1]
m.vivo.N2 = coef(MTSHmod.vivo.N2)[2]
D0.vivo.N2 = 1/v.vivo.N2
Dq.vivo.N2 = log(m.vivo.N2)*D0.vivo.N2
DT.vivo.N2 = (2*Dq.vivo.N2)/(1 - alpha.vivo.N2*D0.vivo.N2)

USCmod.vivo.N2 = USCfunc(alpha.vivo.N2, beta.vivo.N2, D.vivo.N2, DT.vivo.N2, 1)

## Plot
dev.new()
plot(D.vivo.O2, SF.vivo.O2,
     ylab = "Surviving fraction", xlab = "Dose [Gy]",
     ylim = c(min(SF.vivo.O2), 1), xlim = c(0, max(D.vivo.N2)),
     las = 1, pch = 1, log = "y")
points(D.vivo.N2, SF.vivo.N2, pch = 16)
lines(sort(D.vivo.O2), LQmod.vivo.O2[order(D.vivo.O2)], lty = 1)
lines(sort(D.vivo.N2), exp(fitted(LQmod.vivo.N2))[order(D.vivo.N2)], lty = 1)
lines(sort(D.vivo.O2), USCmod.vivo.O2[order(D.vivo.O2)], lty = 3)
lines(sort(D.vivo.N2), exp(USCmod.vivo.N2)[order(D.vivo.N2)], lty = 3)
legend("topright", legend=c(legend.title[2], "Aerobic", "Extremely hypoxic", "LQ", "USC"),
       lty = c(0, 0, 0, 1, 3), pch = c(NA, 1, 16, NA, NA), bty = "o")

# -------------------------------------------------------------------------------------------------- #

SF.vivo.O2 = SFdata[Procedure == "in vivo" & Hypoxia == "No", 1] # & Hypoxia == "No"
D.vivo.O2 = SFdata[Procedure == "in vivo" & Hypoxia == "No", 2]
SF.vivo.N2 = SFdata[Procedure == "in vivo" & Hypoxia == "Yes", 1]
D.vivo.N2 = SFdata[Procedure == "in vivo" & Hypoxia == "Yes", 2]

## Hypoxic model
LQmod.vivo.N2 = lm(log(SF.vivo.N2) ~ -1 + D.vivo.N2 + I(D.vivo.N2^2))
MTSHmod.vivo.N2 = nls(log(SF.vivo.N2) ~ log(1 - (1 - exp(-v*D.vivo.N2) )^m ), start = list(v=0.2, m=2))

alpha.vivo.N2 = (-1)*coef(LQmod.vivo.N2)[1]
beta.vivo.N2 = (-1)*coef(LQmod.vivo.N2)[2]
v.vivo.N2 = coef(MTSHmod.vivo.N2)[1]
m.vivo.N2 = coef(MTSHmod.vivo.N2)[2]
D0.vivo.N2 = 1/v.vivo.N2
Dq.vivo.N2 = log(m.vivo.N2)*D0.vivo.N2
DT.vivo.N2 = (2*Dq.vivo.N2)/(1 - alpha.vivo.N2*D0.vivo.N2)

LQmod.vivo.O2 = nls(SF.vivo.O2 ~ HF*exp(-(alpha.vivo.N2)*D.vivo.O2 - (beta.vivo.N2)*I(D.vivo.O2^2)) 
                    + (1-HF)*exp(-(alpha)*D.vivo.O2 - (beta)*I(D.vivo.O2^2)), 
                    start = list(HF = 0.2, alpha = 0.1, beta = 0.02))

dev.new()
plot(D.vivo.O2, SF.vivo.O2,
     ylab = "Surviving fraction", xlab = "Dose [Gy]",
     ylim = c(min(SF.vivo.O2), 1), xlim = c(0, max(D.vivo.N2)),
     las = 1, pch = 1, log = "y")
points(D.vivo.N2, SF.vivo.N2, pch = 16)
lines(sort(D.vivo.O2), fitted(LQmod.vivo.O2)[order(D.vivo.O2)], lty = 1)
lines(sort(D.vivo.N2), exp(fitted(LQmod.vivo.N2))[order(D.vivo.N2)], lty = 1)
#lines(sort(D.vivo.O2), USCmod.vivo.O2[order(D.vivo.O2)], lty = 3)
#lines(sort(D.vivo.N2), exp(USCmod.vivo.N2)[order(D.vivo.N2)], lty = 3)
legend("topright", legend=c(legend.title[2], "Aerobic", "Extremely hypoxic", "LQ", "USC"),
       lty = c(0, 0, 0, 1, 3), pch = c(NA, 1, 16, NA, NA), bty = "o")


HF = seq(0, 1, length.out = 100)
alpha.O2.span = seq(0.01, 0.5, length.out = 100)
beta.O2.span = seq(0.001, 0.5, length.out = 100)
sumsq = array(rep(NaN, 100*100*100), dim = c(100, 100, 100));
#OER = seq(1, 5, length.out = 100)
#sumsq = matrix(NA, nrow = 100, ncol = 100)   # array(NA, dim = c(100,100,1))

i = 0
j = 0
k = 0
for (hf in HF) {
  
  i = i + 1
  
  for (alpha in alpha.O2.span) {
    
    j = j + 1
    
    for (beta in beta.O2.span) {
      
      k = k + 1
      
      #SFmix = (1-hf)*exp(-(vecalpha.N2[2])*oer*D.test - (vecbeta.N2[2]*oer^2)*I(D.test^2) ) + 
      #  (hf)*exp( -vecalpha.N2[2]*D.test - vecbeta.N2[2]*I(D.test^2) )
      
      SFtest = hf*exp(-(alpha.vivo.N2)*D.vivo.O2 - (beta.vivo.N2)*I(D.vivo.O2^2)) 
      + (1-hf)*exp(-(alpha)*D.vivo.O2 - (beta)*I(D.vivo.O2^2))
      
      sumsq[i,j,k] = RSS(SF.vivo.O2, SFtest)
      
    }
    
    k = 0
    
  }
  
  j = 0
  
}

opt = which(sumsq == min(sumsq))
#HFopt = HF[opt[1]]
#OERopt = OER[opt[2]]
#alpha.vivo.O2 = OERopt*vecalpha.N2[2]
#beta.vivo.O2 = OERopt*vecbeta.N2[2]







