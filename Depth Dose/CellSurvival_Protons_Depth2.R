setwd("C:/Users/Delmon/Dropbox/R/Anne Marit")
getwd()

rm(list=ls())

SFdata = read.table(file = "SF_Protons_Depth.csv", header = TRUE, sep = ";")
attach(SFdata)

Sxx = function(x)   sum( I((x - mean(x))^2) )
s = function(y, yhat, n, p) sqrt( RSS(y, yhat) / (n - p) )
sy = function(x, y, yhat, n, p) s(y, yhat, n, p)*sqrt( (1/n) + (x - mean(x))^2/Sxx(x) )
RSS = function(y ,yhat) sum( I((y - yhat)^2) ) 
ESS = function(y, yhat) sum( (yhat - mean(y))^2 ) # SSR or MSS
TSS = function(y, yhat) sum( (y - mean(y))^2 )
R.squared = function(y, yhat) 1 - (RSS(y, yhat)/TSS(y, yhat))
R.squared.adj = function(y, yhat, n, p) 1 - (1 - R.squared(y, yhat))*( (n - 1)/(n - p) ) 

f = function(D, DT) I(D^2) - (D > DT)*I((D - DT)^2)
USCfunc = function(alpha, beta, D, DT, OER) -(alpha/OER)*D - (beta/OER^2)*I( f(D, DT) )
LQfunc = function(alpha, beta, D, OER) -(alpha/OER)*D - (beta/OER^2)*I(D^2)

SF = SFdata[, 1] # P == "P9" | P == "P10"
D = SFdata[, 2]
P = SFdata[, 3]

LQmod = lm(log(SF) ~ -1 + D + I(D^2) )
MTSHmod = nls(log(SF) ~ log(1 - (1 - exp(-v*D) )^m ), start = list(v=1, m=1))

# USC
alpha = (-1)*coefficients(LQmod)[1]
beta = (-1)*coefficients(LQmod)[2]
D0 = 1/coefficients(MTSHmod)[1]
Dq = log(coefficients(MTSHmod)[2])*D0
DT = (2*Dq)/(1 - alpha*D0)

USCmod = USCfunc(alpha, beta, D, DT, 1)

dev.new()
plot(D, SF,
     ylab = "Survival fraction", xlab = "Dose [Gy]",
     las = 1, pch = 1, log = "y")
lines(sort(D), exp(fitted(LQmod))[order(D)], lty = 1)
lines(sort(D), exp(fitted(MTSHmod))[order(D)], lty = 2)
lines(sort(D), exp(USCmod)[order(D)], lty = 3)
legend("topright", legend=c("LQ", "MTSH", "USC"), 
       lty = c(1, 2, 3), pch = c(1, 1, 1), bty = "o")

R.squared.adj.USC = R.squared.adj(log(SF), USCmod, length(D), 2)

# LQresid = resid(LQmod)
# MTSHresid = resid(MTSHmod)
# USCresid = log(SF) - USC
# 
# dev.new()
# plot(exp(fitted(LQmod)), LQresid,
#      ylab = "Residuals", xlab = "Estimated survival fraction (fitted values)",
#      ylim = c(-0.9, 0.4), las = 1, pch = 1)
# points(exp(fitted(MTSHmod)), MTSHresid, pch = 2)
# points(exp(USC), USCresid, pch = 3)
# abline(0, 0, lty = 2)
# 
# LOESS.LQresid = loess(LQresid ~ exp(fitted(LQmod)) )
# LOESS.MTSHresid = loess(MTSHresid ~ exp(fitted(MTSHmod)) )
# LOESS.USCresid = loess(USCresid ~ exp(USC) )
# 
# lines(sort(exp(fitted(LQmod))), fitted(LOESS.LQresid)[order(exp(fitted(LQmod)))], lty = 1)
# lines(sort(exp(fitted(MTSHmod))), fitted(LOESS.MTSHresid)[order(exp(fitted(MTSHmod)))], lty = 2)
# lines(sort(exp(USC)), fitted(LOESS.USCresid)[order(exp(USC))], lty = 3)
# legend("topright", legend = c("LQ", "MTSH", "USC"), lty = c(1, 2, 3), pch = c(1, 2, 3), bty = "o")


