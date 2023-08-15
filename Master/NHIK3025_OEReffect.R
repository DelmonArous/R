setwd("C:/Users/Delmon/Dropbox/R/Master")
getwd()

# Import .csv-data from Excel 
sfData_MitoseAerobRoom = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=",", nrow=5)
sfData_MitoseAerobIce = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=",", nrow=5, skip=6)
sfData_G1AerobRoom = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=",", nrow=8, skip=13)
sfData_SAerobRoom = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=",", nrow=8, skip=23)
sfData_AsynchAerob = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=",", nrow=7, skip=33)
sfData_G2_15hHyp = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=",", nrow=7, skip=42)
sfData_G1_5hHyp = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=";", nrow=9, skip=51)
sfData_S_13hHyp = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=";", nrow=7, skip=62)
sfData_Mitose_0hHyp = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=";", nrow=5, skip=71)
sfData_expAsynchHyp = read.table(file="NHIK3025_synkronisert_v2.csv", header=TRUE, sep=";", nrow=7, skip=78)

sfData_MitoseAerobRoom
sfData_MitoseAerobIce
sfData_G1AerobRoom
sfData_SAerobRoom
sfData_AsynchAerob
sfData_G2_15hHyp
sfData_G1_5hHyp
sfData_S_13hHyp
sfData_Mitose_0hHyp
sfData_expAsynchHyp

rm(list=ls())

# Scatterplott 
# Appropriate for examining the relationship between 2 numeric variables
# --------------------------------------------------------------------- #
# sf is the outcome variable (Y), dependent variable
# Dose.Gy. is the numerical explanatory variable (X)
# cex: point size 

# ----------------------------------------------------------------------#
# Mitosis
attach(sfData_MitoseAerobRoom)
newDose_aer = rep(Dose.Gy., 3)
newF_aer = c(f1, f2, f3)
par(mfrow=c(1,1))
plot(newDose_aer[1:5], f1, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim = c(0.001, 1), xlim = c(0,10), 
     las=1, pch=1, log="y")
points(newDose_aer[1:5], f2, pch=1)
points(newDose_aer[1:5], f3, pch=1)
detach(sfData_MitoseAerobRoom)

attach(sfData_Mitose_0hHyp)
newDose_hyp = rep(Dose.Gy., 2)
newF_hyp = c(F3823,F3824)
points(newDose_hyp[1:5], F3823, pch=16)
points(newDose_hyp[1:5], F3824, pch=16)
detach(sfData_Mitose_0hHyp)

mod_aer = lm(log(newF_aer) ~ newDose_aer + I(newDose_aer^2))
mod_hyp = lm(log(newF_hyp) ~ newDose_hyp + I(newDose_hyp^2))
targetmod_aer = nls(log(newF_aer) ~ log(1 - (1 - exp(-v*newDose_aer) )^m ), start = list(v=1, m=1), na.action = na.exclude)
targetmod_hyp = nls(log(newF_hyp) ~ log(1 - (1 - exp(-v*newDose_hyp) )^m ), start = list(v=1, m=1), na.action = na.exclude)

lines(smooth.spline(newDose_aer[1:5], I(exp(predict(mod_aer)[1:5]))), lty=1)
lines(smooth.spline(newDose_hyp[1:5], I(exp(predict(mod_hyp)[1:5]))), lty=1)
lines(smooth.spline(newDose_aer[1:5], I(exp(predict(targetmod_aer)[1:5]))), lty=2)
lines(smooth.spline(newDose_hyp[1:5], I(exp(predict(targetmod_hyp)[1:5]))), lty=2)
legend("topright", legend=c("Mitosis (0 hours)", "Aerobic", "Extremely hypoxic", "Multitarget-singel hit"), 
       lty=c(0,1,1,2), pch=c(NA,1,16,NA), bty="o")

#prd = predict(mod_aer, level=0.95, interval = "confidence")
#prd = predict(mod_hyp, level=0.95, interval = "confidence")
#lines(sort(newDose), exp(prd[1:5, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[1:5, 3][order(newDose)]), col=2, lty=2)
#paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
#paste("R-Sq =", round(summary(mod)$r.squared, 4)) 

# Residual standard error: standard deviations of the errors/residuals  
# average/typical sized error/residual variations of observations variasjoner av observasjoner rundt linjen, sqrt(MSE)
summary(mod_aer)
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_aer)$adj.r.squared
summary(mod_hyp)
confint(mod_hyp, conf.level=0.95)
summary(mod_hyp)$adj.r.squared
summary(targetmod_aer)
summary(targetmod_hyp)

# Diagnostics
resids_aer = resid(mod_aer)
resids_hyp = resid(mod_hyp)
targetresids_aer = resid(targetmod_aer)
targetresids_hyp = resid(targetmod_hyp)

# Quantile normal (Q-Q) plot
par(mfrow=c(1,2))
qqnorm(resids_aer, las=1, ylab="Standardized residuals e*", pch=1)
qqline(resids_aer)
legend("topleft", legend=c("Mitosis (0 hours)", "Aerobic"), lty=c(0,1), pch=c(0,1), bty="o")
qqnorm(resids_hyp, las=1, ylab="Standardized residuals e*", pch=16)
qqline(resids_hyp)
legend("topleft", legend=c("Mitosis (0 hours)", "Extremely hypoxic"), lty=c(0,1), pch=c(0,16), bty="o")

# Residual plot
par(mfrow=c(1,1))
matplot(exp(fitted(mod_aer)), resids_aer,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        xlim=c(0, 0.6), las=1, pch=1)
points(exp(fitted(targetmod_aer)), targetresids_aer, pch=2)
abline(0, 0, lty=2)

LOESSresids_aer = loess(resids_aer ~ fitted(mod_aer) )
LOESStargetresids_aer = loess(targetresids_aer ~ fitted(targetmod_aer) )

lines(sort(exp(fitted(mod_aer))), predict(LOESSresids_aer)[order(exp(fitted(mod_aer)))], lty=1)
lines(sort(exp(fitted(targetmod_aer))), predict(LOESStargetresids_aer)[order(exp(fitted(targetmod_aer)))], lty=3)
legend("topright", legend=c("Mitosis (0 hours); aerobic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

matplot(exp(fitted(mod_hyp)), resids_hyp,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_hyp)), targetresids_hyp, pch=2)
abline(0, 0, lty=2)

LOESSresids_hyp = loess(resids_hyp ~ fitted(mod_hyp) )
LOESStargetresids_hyp = loess(targetresids_hyp ~ fitted(targetmod_hyp) )

lines(sort(exp(fitted(mod_hyp))), predict(LOESSresids_hyp)[order(exp(fitted(mod_hyp)))], lty=1)
lines(sort(exp(fitted(targetmod_hyp))), predict(LOESStargetresids_hyp)[order(exp(fitted(targetmod_hyp)))], lty=3)
legend("topright", legend=c("Mitosis (0 hours); extremely hypoxic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

#S_aer = sqrt(deviance(mod_aer)/df.residual(mod_aer))
#targetS_aer = sqrt(deviance(targetmod_aer)/df.residual(targetmod_aer))
#S_hyp = sqrt(deviance(mod_hyp)/df.residual(mod_hyp))
#targetS_hyp = sqrt(deviance(targetmod_hyp)/df.residual(targetmod_hyp))

# ----------------------------------------------------------------------#
# G1
attach(sfData_G1AerobRoom)
newDose_aer = rep(Dose.Gy., 4)
newF_aer = c(F2725,F2828,F2731,F2747)
plot(newDose_aer[9:16], F2725, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim=c(0.0001, 1), xlim=c(0,27),
     las=1, pch=1, log="y")
points(newDose_aer[9:16], F2828, pch=1)
points(newDose_aer[9:16], F2731, pch=1)
points(newDose_aer[9:16], F2747, pch=1)
detach(sfData_G1AerobRoom)

attach(sfData_G1_5hHyp)
newDose_hyp = rep(Dose.Gy., 5)
newF_hyp = c(F3797,F3799,F3804,F3848,F3852)
points(newDose_hyp[1:9], F3797, pch=16)
points(newDose_hyp[1:9], F3799, pch=16)
points(newDose_hyp[1:9], F3804, pch=16)
points(newDose_hyp[1:9], F3848, pch=16)
points(newDose_hyp[1:9], F3852, pch=16)
detach(sfData_G1_5hHyp)

mod_aer = lm(log(newF_aer) ~ newDose_aer + I(newDose_aer^2), na.action = na.exclude)
mod_hyp = lm(log(newF_hyp) ~ newDose_hyp + I(newDose_hyp^2), na.action = na.exclude)
tempvec_aer = predict(mod_aer)[9:16]
tempvec_hyp = c(predict(mod_hyp)[1:4], predict(mod_hyp)[23:27])

targetmod_aer = nls(log(newF_aer) ~ log(1 - (1 - exp(-v*newDose_aer) )^m ), start = list(v=1, m=1), na.action = na.exclude)
targetmod_hyp = nls(log(newF_hyp) ~ log(1 - (1 - exp(-v*newDose_hyp) )^m ), start = list(v=1, m=1), na.action = na.exclude)
targettempvec_aer = predict(targetmod_aer)[9:16]
targettempvec_hyp = c(predict(targetmod_hyp)[1:4], predict(targetmod_hyp)[23:27])

lines(smooth.spline(newDose_aer[9:16], exp(tempvec_aer)), lty=1)
lines(smooth.spline(newDose_hyp[1:9], exp(tempvec_hyp)), lty=1)
lines(smooth.spline(newDose_aer[9:16], exp(targettempvec_aer)), lty=2)
lines(smooth.spline(newDose_hyp[1:9], exp(targettempvec_hyp)), lty=2)
legend("topright", legend=c("G1 (5 hours)", "Aerobic", "Extremely hypoxic", "Multitarget-singel hit"),
       lty=c(NA,1,1,2), pch=c(NA,1,16,NA), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")[9:16]
#lines(sort(newDose), exp(prd[9:16, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[9:16, 3][order(newDose)]), col=2, lty=2)

summary(mod_aer)
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_aer)$adj.r.squared
summary(mod_hyp)
confint(mod_hyp, conf.level=0.95)
summary(mod_hyp)$adj.r.squared
summary(targetmod_aer)
summary(targetmod_hyp)

# Diagnostics
resids_aer = resid(mod_aer)
resids_hyp = resid(mod_hyp)
targetresids_aer = resid(targetmod_aer)
targetresids_hyp = resid(targetmod_hyp)

# Quantile normal (Q-Q) plot for LQ
par(mfrow=c(1,2))
qqnorm(resids_aer, las=1, ylab="Standardized residuals e*", pch=1)
qqline(resids_aer)
legend("topleft", legend=c("G1 (5 hours)", "Aerobic"), lty=c(0,1), pch=c(0,1), bty="o")
qqnorm(resids_hyp, las=1, ylab="Standardized residuals e*", pch=16)
qqline(resids_hyp)
legend("topleft", legend=c("G1 (5 hours)", "Extremely hypoxic"), lty=c(0,1), pch=c(0,16), bty="o")

# Residual plot
par(mfrow=c(1,1))
matplot(exp(fitted(mod_aer)), resids_aer,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        ylim=c(-1.5, 1), xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_aer)), targetresids_aer, pch=2)
abline(0, 0, lty=2)

LOESSresids_aer = loess(resids_aer[!is.na(resids_aer)] ~ fitted(mod_aer)[!is.na(resids_aer)] )
LOESStargetresids_aer = loess(targetresids_aer[!is.na(targetresids_aer)] ~ fitted(targetmod_aer)[!is.na(targetresids_aer)])

lines(sort(exp(fitted(mod_aer))[!is.na(resids_aer)]), 
      predict(LOESSresids_aer)[order(exp(fitted(mod_aer))[!is.na(resids_aer)])], lty=1)
lines(sort(exp(fitted(targetmod_aer))[!is.na(targetresids_aer)]), 
      predict(LOESStargetresids_aer)[order(exp(fitted(targetmod_aer))[!is.na(targetresids_aer)] )], lty=3)
legend("topright", legend=c("G1 (5 hours); aerobic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

matplot(exp(fitted(mod_hyp)), resids_hyp,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        ylim=c(-1.5, 1.5), xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_hyp)), targetresids_hyp, pch=2)
abline(0, 0, lty=2)

LOESSresids_hyp = loess(resids_hyp[!is.na(resids_hyp)] ~ fitted(mod_hyp)[!is.na(resids_hyp)] )
LOESStargetresids_hyp = loess(targetresids_hyp[!is.na(targetresids_hyp)] ~ fitted(targetmod_hyp)[!is.na(targetresids_hyp)])

lines(sort(exp(fitted(mod_hyp))[!is.na(resids_hyp)]), 
      predict(LOESSresids_hyp)[order(exp(fitted(mod_hyp))[!is.na(resids_hyp)])], lty=1)
lines(sort(exp(fitted(targetmod_hyp))[!is.na(targetresids_hyp)]), 
      predict(LOESStargetresids_hyp)[order(exp(fitted(targetmod_hyp))[!is.na(targetresids_hyp)] )], lty=3)
legend("topright", legend=c("G1 (5 hours); extremely hypoxic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

# ----------------------------------------------------------------------#
# S
attach(sfData_SAerobRoom)
newDose_aer = rep(Dose.Gy., 4)
newF_aer = c(F2725,F2829,F2732,F2743)
plot(newDose_aer[9:16], F2725, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim=c(0.0001, 1), xlim=c(0,22),
     las=1, pch=1, log="y")
points(newDose_aer[9:16], F2829, pch=1)
points(newDose_aer[9:16], F2732, pch=1)
points(newDose_aer[9:16], F2743, pch=1)
detach(sfData_SAerobRoom)

attach(sfData_S_13hHyp)
newDose_hyp = rep(Dose.Gy., 4)
newF_hyp = c(F3797,F3799,F3804,F3848)
points(newDose_hyp[8:14], F3797, pch=16)
points(newDose_hyp[8:14], F3799, pch=16)
points(newDose_hyp[8:14], F3804, pch=16)
points(newDose_hyp[8:14], F3848, pch=16)
detach(sfData_S_13hHyp)

mod_aer = lm(log(newF_aer) ~ newDose_aer + I(newDose_aer^2), na.action = na.exclude)
mod_hyp = lm(log(newF_hyp) ~ newDose_hyp + I(newDose_hyp^2), na.action = na.exclude)
tempvec_aer = predict(mod_aer)[9:16]
tempvec_hyp = predict(mod_hyp)[8:14]

targetmod_aer = nls(log(newF_aer) ~ log(1 - (1 - exp(-v*newDose_aer) )^m ), start = list(v=1, m=1), na.action = na.exclude)
targetmod_hyp = nls(log(newF_hyp) ~ log(1 - (1 - exp(-v*newDose_hyp) )^m ), start = list(v=1, m=1), na.action = na.exclude)
targettempvec_aer = predict(targetmod_aer)[9:16]
targettempvec_hyp = predict(targetmod_hyp)[8:14]

lines(smooth.spline(newDose_aer[9:16], exp(tempvec_aer)), lty=1)
lines(smooth.spline(newDose_hyp[8:14], exp(tempvec_hyp)), lty=1)
lines(smooth.spline(newDose_aer[9:16], exp(targettempvec_aer)), lty=2)
lines(smooth.spline(newDose_hyp[8:14], exp(targettempvec_hyp)), lty=2)
legend("topright", legend=c("S (13 hours)", "Aerobic", "Extremely hypoxic", "Multitarget-singel hit"),
       lty=c(NA,1,1,2), pch=c(NA,1,16,NA), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")
#lines(sort(Dose.Gy.), exp(prd[9:16, 2][order(Dose.Gy.)]), col=2, lty=2)
#lines(sort(Dose.Gy.), exp(prd[9:16, 3][order(Dose.Gy.)]), col=2, lty=2)

summary(mod_aer)
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_aer)$r.squared
summary(mod_aer)$adj.r.squared
summary(mod_hyp)
confint(mod_hyp, conf.level=0.95)
summary(mod_hyp)$r.squared
summary(mod_hyp)$adj.r.squared
summary(targetmod_aer)
summary(targetmod_hyp)

# Diagnostics
resids_aer = resid(mod_aer)
resids_hyp = resid(mod_hyp)
targetresids_aer = resid(targetmod_aer)
targetresids_hyp = resid(targetmod_hyp)

# Quantile normal (Q-Q) plot for LQ
par(mfrow=c(1,2))
qqnorm(resids_aer, las=1, ylab="Standardized residuals e*", pch=1)
qqline(resids_aer)
legend("topleft", legend=c("S (13 hours)", "Aerobic"), lty=c(0,1), pch=c(0,1), bty="o")
qqnorm(resids_hyp, las=1, ylab="Standardized residuals e*", pch=16)
qqline(resids_hyp)
legend("topleft", legend=c("S (13 hours)", "Extremely hypoxic"), lty=c(0,1), pch=c(0,16), bty="o")

# Residual plot
par(mfrow=c(1,1))
matplot(exp(fitted(mod_aer)), resids_aer,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        ylim=c(-1.6, 1.5), xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_aer)), targetresids_aer, pch=2)
abline(0, 0, lty=2)

LOESSresids_aer = loess(resids_aer[!is.na(resids_aer)] ~ fitted(mod_aer)[!is.na(resids_aer)] )
LOESStargetresids_aer = loess(targetresids_aer[!is.na(targetresids_aer)] ~ fitted(targetmod_aer)[!is.na(targetresids_aer)])

lines(sort(exp(fitted(mod_aer))[!is.na(resids_aer)]), 
      predict(LOESSresids_aer)[order(exp(fitted(mod_aer))[!is.na(resids_aer)])], lty=1)
lines(sort(exp(fitted(targetmod_aer))[!is.na(targetresids_aer)]), 
      predict(LOESStargetresids_aer)[order(exp(fitted(targetmod_aer))[!is.na(targetresids_aer)] )], lty=3)
legend("topright", legend=c("S (13 hours); aerobic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

matplot(exp(fitted(mod_hyp)), resids_hyp,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        ylim=c(-1.5, 1), xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_hyp)), targetresids_hyp, pch=2)
abline(0, 0, lty=2)

LOESSresids_hyp = loess(resids_hyp[!is.na(resids_hyp)] ~ fitted(mod_hyp)[!is.na(resids_hyp)] )
LOESStargetresids_hyp = loess(targetresids_hyp[!is.na(targetresids_hyp)] ~ fitted(targetmod_hyp)[!is.na(targetresids_hyp)])

lines(sort(exp(fitted(mod_hyp))[!is.na(resids_hyp)]), 
      predict(LOESSresids_hyp)[order(exp(fitted(mod_hyp))[!is.na(resids_hyp)])], lty=1)
lines(sort(exp(fitted(targetmod_hyp))[!is.na(targetresids_hyp)]), 
      predict(LOESStargetresids_hyp)[order(exp(fitted(targetmod_hyp))[!is.na(targetresids_hyp)] )], lty=3)
legend("topright", legend=c("S (13 hours); extremely hypoxic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

# ----------------------------------------------------------------------#
# G2; 15 hours (hypoxic)
attach(sfData_G2_15hHyp)
newDose = rep(Dose.Gy., 2)
newF = c(F3869,F3870)
plot(Dose.Gy., F3869, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim=c(0.001, 1), xlim=c(0,22),
     las=1, pch=16, log="y")
points(Dose.Gy., F3870, pch=16)
detach(sfData_G2_15hHyp)

mod_hyp = lm(log(newF) ~ newDose + I(newDose^2))
targetmod_hyp = nls(log(newF) ~ log(1 - (1 - exp(-v*newDose) )^m ), start = list(v=1, m=1), na.action = na.exclude)

lines(smooth.spline(newDose[1:7], I(exp(predict(mod_hyp)[1:7]))), lty=1)
lines(smooth.spline(newDose[1:7], I(exp(predict(targetmod_hyp)[1:7]))), lty=2)
legend("topright", legend=c("G2 (15 hours)", "Extremely hypoxic", "Multitarget-singel hit"),
       lty=c(NA,1,2), pch=c(NA,16,NA), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")
#lines(sort(newDose), exp(prd[1:7, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[1:7, 3][order(newDose)]), col=2, lty=2)

summary(mod)
confint(mod, conf.level=0.95) # confidence intervall for model coefficients
summary(mod)$r.squared
summary(mod)$adj.r.squared
summary(targetmod_hyp)

# Diagnostics
resids_hyp = resid(mod_hyp)
targetresids_hyp = resid(targetmod_hyp)

# Residual plot
matplot(exp(fitted(mod_hyp)), resids_hyp,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        ylim=c(-1,0.5), xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_hyp)), targetresids_hyp, pch=2)
abline(0, 0, lty=2)

LOESSresids_hyp = loess(resids_hyp ~ fitted(mod_hyp) )
LOESStargetresids_hyp = loess(targetresids_hyp ~ fitted(targetmod_hyp) )

lines(sort(exp(fitted(mod_hyp))), predict(LOESSresids_hyp)[order(exp(fitted(mod_hyp)))], lty=1)
lines(sort(exp(fitted(targetmod_hyp))), predict(LOESStargetresids_hyp)[order(exp(fitted(targetmod_hyp)))], lty=3)
legend("topright", legend=c("G2; extremely hypoxic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

# ----------------------------------------------------------------------#
# Asynchronous 
attach(sfData_AsynchAerob)
newDose_aer = rep(Dose.Gy., 4)
newF_aer = c(F2726,F2827,F2733,F2744)
plot(newDose_aer[1:7], F2726, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim=c(0.0001, 1), xlim=c(0,25),
     las=1, pch=1, log="y")
points(newDose_aer[1:7], F2827, pch=1)
points(newDose_aer[1:7], F2733, pch=1)
points(newDose_aer[1:7], F2744, pch=1)
detach(sfData_AsynchAerob)

attach(sfData_expAsynchHyp)
newDose_hyp = rep(Dose.Gy., 3)
newF_hyp = c(F3839,F3840, F3842)
points(newDose_hyp[1:7], F3839, pch=16)
points(newDose_hyp[1:7], F3840, pch=16)
points(newDose_hyp[1:7], F3842, pch=16)
detach(sfData_expAsynchHyp)

mod_aer = lm(log(newF_aer) ~ newDose_aer + I(newDose_aer^2))
mod_hyp = lm(log(newF_hyp) ~ newDose_hyp + I(newDose_hyp^2))
targetmod_aer = nls(log(newF_aer) ~ log(1 - (1 - exp(-v*newDose_aer) )^m ), start = list(v=1, m=1), na.action = na.exclude)
targetmod_hyp = nls(log(newF_hyp) ~ log(1 - (1 - exp(-v*newDose_hyp) )^m ), start = list(v=1, m=1), na.action = na.exclude)

lines(smooth.spline(newDose_aer[1:7], I(exp(predict(mod_aer)[1:7]))), lty=1)
lines(smooth.spline(newDose_hyp[1:7], I(exp(predict(mod_hyp)[1:7]))), lty=1)
lines(smooth.spline(newDose_aer[1:7], I(exp(predict(targetmod_aer)[1:7]))), lty=2)
lines(smooth.spline(newDose_hyp[1:7], I(exp(predict(targetmod_hyp)[1:7]))), lty=2)
legend("topright", legend=c("Asynchronous cells", "Aerobic", "Extremely hypoxic", "Multitarget-singel hit"),
       lty=c(NA,1,1,2), pch=c(NA,1,16,NA), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")
#lines(sort(newDose), exp(prd[1:7, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[1:7, 3][order(newDose)]), col=2, lty=2)

summary(mod_aer)
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_aer)$r.squared
summary(mod_aer)$adj.r.squared
summary(mod_hyp)
confint(mod_hyp, conf.level=0.95)
summary(mod_hyp)$r.squared
summary(mod_hyp)$adj.r.squared
summary(targetmod_aer)
summary(targetmod_hyp)

# Diagnostics
resids_aer = resid(mod_aer)
resids_hyp = resid(mod_hyp)
targetresids_aer = resid(targetmod_aer)
targetresids_hyp = resid(targetmod_hyp)

# Quantile normal (Q-Q) plot for LQ
par(mfrow=c(1,2))
qqnorm(resids_aer, las=1, ylab="Standardized residuals e*", pch=1)
qqline(resids_aer)
legend("topleft", legend=c("Asynchronous cells", "Aerobic"), lty=c(0,1), pch=c(0,1), bty="o")
qqnorm(resids_hyp, las=1, ylab="Standardized residuals e*", pch=16)
qqline(resids_hyp)
legend("topleft", legend=c("Asynchronous cells", "Extremely hypoxic"), lty=c(0,1), pch=c(0,16), bty="o")

# Residual plot
par(mfrow=c(1,1))
matplot(exp(fitted(mod_aer)), resids_aer,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        ylim=c(-2,2), xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_aer)), targetresids_aer, pch=2)
abline(0, 0, lty=2)

LOESSresids_aer = loess(resids_aer ~ fitted(mod_aer) )
LOESStargetresids_aer = loess(targetresids_aer ~ fitted(targetmod_aer) )

lines(sort(exp(fitted(mod_aer))), predict(LOESSresids_aer)[order(exp(fitted(mod_aer)))], lty=1)
lines(sort(exp(fitted(targetmod_aer))), predict(LOESStargetresids_aer)[order(exp(fitted(targetmod_aer)))], lty=3)
legend("topright", legend=c("Asynchronous cells; aerobic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")

matplot(exp(fitted(mod_hyp)), resids_hyp,
        ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)",
        ylim=c(-1,1.3), xlim=c(0, 1), las=1, pch=1)
points(exp(fitted(targetmod_hyp)), targetresids_hyp, pch=2)
abline(0, 0, lty=2)

LOESSresids_hyp = loess(resids_hyp ~ fitted(mod_hyp) )
LOESStargetresids_hyp = loess(targetresids_hyp ~ fitted(targetmod_hyp) )

lines(sort(exp(fitted(mod_hyp))), predict(LOESSresids_hyp)[order(exp(fitted(mod_hyp)))], lty=1)
lines(sort(exp(fitted(targetmod_hyp))), predict(LOESStargetresids_hyp)[order(exp(fitted(targetmod_hyp)))], lty=3)
legend("topright", legend=c("Asynchronous cells; extremely hypoxic", "LQ", "Multitarget-singel hit"), 
       lty=c(NA,1,3), pch=c(NA,1,2), bty="o")