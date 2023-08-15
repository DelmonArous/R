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

lines(smooth.spline(newDose_aer[1:5], I(exp(predict(mod_aer)[1:5]))), lty=1)
lines(smooth.spline(newDose_hyp[1:5], I(exp(predict(mod_hyp)[1:5]))), lty=1)
legend("bottomleft", legend=c("Mitosis (0 hours)", "Aerobic", "Extremely hypoxic"), 
       lty=c(0,1,1), pch=c(0,1,16), bty="o")

#prd = predict(mod_aer, level=0.95, interval = "confidence")
#prd = predict(mod_hyp, level=0.95, interval = "confidence")
#lines(sort(newDose), exp(prd[1:5, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[1:5, 3][order(newDose)]), col=2, lty=2)
#paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
#paste("R-Sq =", round(summary(mod)$r.squared, 4)) 

# Residual standard error: standard deviations of the errors/residuals  
# average/typical sized error/residual variations of observations variasjoner av observasjoner rundt linjen, sqrt(MSE)
summary(mod_aer)
coefficients(mod_aer)[2]  # alpha
coefficients(mod_aer)[3]  # beta
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_hyp)
coefficients(mod_hyp)[2]  # alpha
coefficients(mod_hyp)[3]  # beta
confint(mod_hyp, conf.level=0.95)

modTarget = nls(log(newF) ~ log(1 - (1 - exp(-v*newDose) )^m ), start = list(v=1, m=1))
#modTarget = nls(log(newF) ~ -v*newDose + log(m), start = list(v=1, m=1)) # approx for high dose D
lines(smooth.spline(newDose, predict(modTarget)), col=3, lty=1)

# Diagnostics

# Quantile normal (Q-Q) plot for LQ
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Quantile normal (Q-Q) plot for multi-target
sresidsTarget = resid(modTarget)
qqnorm(sresidsTarget, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresidsTarget, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", las=1, col=2, pch=2)
points(fitted(modTarget), sresidsTarget, col=3, pch=3)
abline(0, 0, lty=1)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                     fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
sresidsModTarget = lm(sresidsTarget[sresidsTarget >= -S & sresidsTarget <= S] ~ 
                  fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S] + I(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids>=-S & sresids <= S], predict(sresidsMod)), col=2, lty=2)
lines(smooth.spline(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S], predict(sresidsModTarget)), col=3, lty=3)
legend("topright", legend=c("LQ model", "Multi-target model"), col=c(2,3), lty=c(2,3), pch=c(2,3), bty="o")

#cor(Dose.Gy., f1)       # Pearson's correlation
#cov(Dose.Gy., f1)       # Covariance
#var(f1)                 # Variance
#sd(f1)                  # Standard deviation
#mean(f1, trim=0.10)     # Remove top and bottom 10% of the obervations


detach(sfData_MitoseAerobRoom)

# ----------------------------------------------------------------------#
# G1
attach(sfData_G1AerobRoom)
newDose_aer = rep(Dose.Gy., 4)
newF_aer = c(F2725,F2828,F2731,F2747)
plot(newDose_aer[9:16], F2725, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim=c(0.0005, 1), xlim=c(0,27),
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

lines(smooth.spline(newDose_aer[9:16], exp(tempvec_aer)), lty=1)
lines(smooth.spline(newDose_hyp[1:9], exp(tempvec_hyp)), lty=1)
legend("topright", legend=c("G1 (5 hours)", "Aerobic", "Extremely hypoxic"),
       lty=c(0,1,1), pch=c(0,1,16), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")[9:16]
#lines(sort(newDose), exp(prd[9:16, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[9:16, 3][order(newDose)]), col=2, lty=2)

summary(mod_aer)
coefficients(mod_aer)[2]  # alpha
coefficients(mod_aer)[3]  # beta
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_aer)$r.squared
summary(mod_aer)$adj.r.squared
summary(mod_hyp)
coefficients(mod_hyp)[2]  # alpha
coefficients(mod_hyp)[3]  # beta
confint(mod_hyp, conf.level=0.95)
summary(mod_hyp)$r.squared
summary(mod_hyp)$adj.r.squared

# Diagnostics

# Quantile normal (Q-Q) plot
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=1)
abline(0, 0, lty=2)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids>=-S & sresids <= S], predict(sresidsMod)), col=1, lty=1)

detach(sfData_G1AerobRoom)

# ----------------------------------------------------------------------#
# Mitosis (aerobic at ice temperature)
attach(sfData_MitoseAerobIce)

newDose = rep(Dose.Gy., 3)
newF = c(f1, f2, f3)

plot(Dose.Gy., f1, main="Mitosis (aerobic, ice temperature)", 
     ylab="Survival fraction", xlab="Dose [Gy]", ylim = c(min(newF), max(newF)),
     las=1, col=2, pch=1, log="y")
points(Dose.Gy., f2, col=2, pch=1)
points(Dose.Gy., f3, col=2, pch=1)

mod = lm(log(newF) ~ newDose + I(newDose^2))
lines(smooth.spline(Dose.Gy., I(exp(predict(mod)[1:5]))), col=1, lty=1)
prd = predict(mod, level=0.95, interval = "confidence")
lines(sort(newDose), exp(prd[1:5, 2][order(newDose)]), col=2, lty=2)
lines(sort(newDose), exp(prd[1:5, 3][order(newDose)]), col=2, lty=2)
legend("bottomleft", legend=c("Regression line", "Confidence interval (95%)", 
                              paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
                              paste("R-Sq =", round(summary(mod)$r.squared, 4)) ), 
       col=c(1,2, 0, 0), lty=c(1,2,0, 0), bty="o")

confint(mod, conf.level=0.95) # confidence intervall for model coefficients
coefficients(mod)[2]  # alpha
coefficients(mod)[3]  # beta
summary(mod)$r.squared
summary(mod)$adj.r.squared
summary(mod)

modTarget = nls(log(newF) ~ log(1 - (1 - exp(-v*newDose) )^m ), start = list(v=1, m=1))
#modTarget = nls(log(newF) ~ -v*newDose + log(m), start = list(v=1, m=1)) # approx for high dose D
lines(smooth.spline(newDose, predict(modTarget)), col=3, lty=1)

# Diagnostics

# Quantile normal (Q-Q) plot for LQ
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Quantile normal (Q-Q) plot for multi-target
sresidsTarget = resid(modTarget)
qqnorm(sresidsTarget, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresidsTarget, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=2)
points(fitted(modTarget), sresidsTarget, col=3, pch=3)
abline(0, 0, lty=1)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
sresidsModTarget = lm(sresidsTarget[sresidsTarget >= -S & sresidsTarget <= S] ~ 
                        fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S] + I(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids>=-S & sresids <= S], predict(sresidsMod)), col=2, lty=2)
lines(smooth.spline(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S], predict(sresidsModTarget)), col=3, lty=3)
legend("topright", legend=c("LQ model", "Multi-target model"), col=c(2,3), lty=c(2,3), pch=c(2,3), bty="o")

detach(sfData_MitoseAerobIce)

# ----------------------------------------------------------------------#
# S
attach(sfData_SAerobRoom)
newDose_aer = rep(Dose.Gy., 4)
newF_aer = c(F2725,F2829,F2732,F2743)
plot(newDose_aer[9:16], F2725, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim=c(0.0001, 1), xlim=c(0,25),
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

lines(smooth.spline(newDose_aer[9:16], exp(tempvec_aer)), lty=1)
lines(smooth.spline(newDose_hyp[8:14], exp(tempvec_hyp)), lty=1)
legend("topright", legend=c("S (13 hours)", "Aerobic", "Extremely hypoxic"),
       lty=c(0,1,1), pch=c(0,1,16), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")
#lines(sort(Dose.Gy.), exp(prd[9:16, 2][order(Dose.Gy.)]), col=2, lty=2)
#lines(sort(Dose.Gy.), exp(prd[9:16, 3][order(Dose.Gy.)]), col=2, lty=2)

summary(mod_aer)
coefficients(mod_aer)[2]  # alpha
coefficients(mod_aer)[3]  # beta
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_aer)$r.squared
summary(mod_aer)$adj.r.squared
summary(mod_hyp)
coefficients(mod_hyp)[2]  # alpha
coefficients(mod_hyp)[3]  # beta
confint(mod_hyp, conf.level=0.95)
summary(mod_hyp)$r.squared
summary(mod_hyp)$adj.r.squared

# ----------------------------------------------------------------------#
# Asynchronous 
attach(sfData_AsynchAerob)
newDose_aer = rep(Dose.Gy., 4)
newF_aer = c(F2726,F2827,F2733,F2744)
plot(newDose_aer[1:7], F2726, 
     ylab="Survival fraction", xlab="Dose [Gy]", 
     ylim=c(0.0001, 1), xlim=c(0,30),
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

lines(smooth.spline(newDose_aer[1:7], I(exp(predict(mod_aer)[1:7]))), lty=1)
lines(smooth.spline(newDose_hyp[1:7], I(exp(predict(mod_hyp)[1:7]))), lty=1)
legend("topright", legend=c("Asynchronous cells", "Aerobic", "Extremely hypoxic"),
       lty=c(0,1,1), pch=c(0,1,16), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")
#lines(sort(newDose), exp(prd[1:7, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[1:7, 3][order(newDose)]), col=2, lty=2)

summary(mod_aer)
coefficients(mod_aer)[2]  # alpha
coefficients(mod_aer)[3]  # beta
confint(mod_aer, conf.level=0.95) # confidence intervall for model coefficients
summary(mod_aer)$r.squared
summary(mod_aer)$adj.r.squared
summary(mod_hyp)
coefficients(mod_hyp)[2]  # alpha
coefficients(mod_hyp)[3]  # beta
confint(mod_hyp, conf.level=0.95)
summary(mod_hyp)$r.squared
summary(mod_hyp)$adj.r.squared

modTarget = nls(log(newF) ~ log(1 - (1 - exp(-v*newDose) )^m ), start = list(v=1, m=1))
#modTarget = nls(log(newF) ~ -v*newDose + log(m), start = list(v=1, m=1)) # approx for high dose D
lines(smooth.spline(newDose, predict(modTarget)), col=3, lty=1)

# Diagnostics

# Quantile normal (Q-Q) plot for LQ
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Quantile normal (Q-Q) plot for multi-target
sresidsTarget = resid(modTarget)
qqnorm(sresidsTarget, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresidsTarget, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=2)
points(fitted(modTarget), sresidsTarget, col=3, pch=3)
abline(0, 0, lty=1)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
sresidsModTarget = lm(sresidsTarget[sresidsTarget >= -S & sresidsTarget <= S] ~ 
                        fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S] + I(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids >= -S & sresids <= S], predict(sresidsMod)), col=2, lty=2)
lines(smooth.spline(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S], predict(sresidsModTarget)), col=3, lty=3)
legend("topright", legend=c("LQ model", "Multi-target model"), col=c(2,3), lty=c(2,3), pch=c(2,3), bty="o")

detach(sfData_AsynchAerob)

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

mod = lm(log(newF) ~ newDose + I(newDose^2))

lines(smooth.spline(newDose[1:7], I(exp(predict(mod)[1:7]))), lty=1)
legend("topright", legend=c("G2 (15 hours)", "Extremely hypoxic"),
       lty=c(0,1), pch=c(0,16), bty="o")

#prd = predict(mod, level=0.95, interval = "confidence")
#lines(sort(newDose), exp(prd[1:7, 2][order(newDose)]), col=2, lty=2)
#lines(sort(newDose), exp(prd[1:7, 3][order(newDose)]), col=2, lty=2)

summary(mod)
coefficients(mod)[2]  # alpha
coefficients(mod)[3]  # beta
confint(mod, conf.level=0.95) # confidence intervall for model coefficients
summary(mod)$r.squared
summary(mod)$adj.r.squared


modTarget = nls(log(newF) ~ log(1 - (1 - exp(-v*newDose) )^m ), start = list(v=1, m=1))
#modTarget = nls(log(newF) ~ -v*newDose + log(m), start = list(v=1, m=1)) # approx for high dose D
lines(smooth.spline(newDose, predict(modTarget)), col=3, lty=1)

# Diagnostics

# Quantile normal (Q-Q) plot for LQ
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Quantile normal (Q-Q) plot for multi-target
sresidsTarget = resid(modTarget)
qqnorm(sresidsTarget, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresidsTarget, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=2)
points(fitted(modTarget), sresidsTarget, col=3, pch=3)
abline(0, 0, lty=1)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
sresidsModTarget = lm(sresidsTarget[sresidsTarget >= -S & sresidsTarget <= S] ~ 
                        fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S] + I(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids >= -S & sresids <= S], predict(sresidsMod)), col=2, lty=2)
lines(smooth.spline(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S], predict(sresidsModTarget)), col=3, lty=3)
legend("topright", legend=c("LQ model", "Multi-target model"), col=c(2,3), lty=c(2,3), pch=c(2,3), bty="o")


detach(sfData_G2_15hHyp)

# ----------------------------------------------------------------------#
# G1; 5 hours (hypoxic)
attach(sfData_G1_5hHyp)

newDose = rep(Dose.Gy., 5)
newF = c(F3797,F3799,F3804,F3848,F3852)

plot(Dose.Gy., F3797, main="G1; 5 hours (hypoxic)", 
     ylab="Survival fraction", xlab="Dose [Gy]", ylim=c(0.00061, 1),
     las=1, col=2, pch=1, log="y")
points(Dose.Gy., F3799, col=2, pch=1)
points(Dose.Gy., F3804, col=2, pch=1)
points(Dose.Gy., F3848, col=2, pch=1)
points(Dose.Gy., F3852, col=2, pch=1)

mod = lm(log(newF) ~ newDose + I(newDose^2), na.action = na.exclude)
tempvec = c(predict(mod)[1:4], predict(mod)[23:27])
lines(smooth.spline(Dose.Gy., exp(tempvec)), col=1, lty=1)
prd = predict(mod, level=0.95, interval = "confidence")
lines(Dose.Gy., exp(tempvec), col=2, lty=2)
lines(Dose.Gy., exp(tempvec), col=2, lty=2)
legend("bottomleft", legend=c("Regression line", "Confidence interval (95%)", 
                              paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
                              paste("R-Sq =", round(summary(mod)$r.squared, 4)) ), 
       col=c(1,2, 0, 0), lty=c(1,2,0, 0), bty="o")


matplot(Dose.Gy., log(cbind(F3797,F3799,F3804,F3848,F3852)),
        main="G1; 5 hours (hypoxic)", 
        ylab="Survival fraction", xlab="Dose [Gy]", 
        las=1, col=2, pch=1)


mod = lm(log(newF) ~ newDose + I(newDose^2), na.action = na.exclude)
tempvec = c(newF[1:4], newF[23:27])
lines(sort(newDose), fitted(mod)[order(newDose)], col=1, lty=1)
prd = predict(mod, level=0.95, interval = "confidence")
lines(sort(newDose), prd[, 2][order(newDose)], col=2, lty=2)
lines(sort(newDose), prd[, 3][order(newDose)], col=2, lty=2)
legend("bottomleft", legend=c("Regression line", "Confidence interval (95%)", 
                              paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
                              paste("R-Sq =", round(summary(mod)$r.squared, 4)) ), 
       col=c(1,2, 0, 0), lty=c(1,2,0, 0), bty="o")

confint(mod, conf.level=0.95) # confidence intervall for model coefficients
coefficients(mod)[2]  # alpha
coefficients(mod)[3]  # beta
summary(mod)$r.squared
summary(mod)$adj.r.squared
summary(mod)

# Diagnostics

# Quantile normal (Q-Q) plot
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=1)
abline(0, 0, lty=2)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids>=-S & sresids <= S], predict(sresidsMod)), col=1, lty=1)

detach(sfData_G1_5hHyp)

# ----------------------------------------------------------------------#
# S; 13 hours (hypoxic)
attach(sfData_S_13hHyp)

newDose = rep(Dose.Gy., 4)
newF = c(F3797,F3799,F3804,F3848)

plot(Dose.Gy., F3797, main="S; 13 hours (hypoxic)", 
     ylab="Survival fraction", xlab="Dose [Gy]", ylim=c(0.0005, 1),
     las=1, col=2, pch=1, log="y")
points(Dose.Gy., F3799, col=2, pch=1)
points(Dose.Gy., F3804, col=2, pch=1)
points(Dose.Gy., F3848, col=2, pch=1)

mod = lm(log(newF) ~ newDose + I(newDose^2), na.action = na.exclude)
tempvec = predict(mod)[8:14]
lines(smooth.spline(Dose.Gy., exp(tempvec)), col=1, lty=1)
prd = predict(mod, level=0.95, interval = "confidence")
lines(sort(Dose.Gy.), exp(prd[8:14, 2][order(Dose.Gy.)]), col=2, lty=2)
lines(sort(Dose.Gy.), exp(prd[8:14, 3][order(Dose.Gy.)]), col=2, lty=2)
legend("bottomleft", legend=c("Regression line", "Confidence interval (95%)", 
                              paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
                              paste("R-Sq =", round(summary(mod)$r.squared, 4)) ), 
       col=c(1,2, 0, 0), lty=c(1,2,0, 0), bty="o")

confint(mod, conf.level=0.95) # confidence intervall for model coefficients
coefficients(mod)[2]  # alpha
coefficients(mod)[3]  # beta
summary(mod)$r.squared
summary(mod)$adj.r.squared
summary(mod)

# Diagnostics

# Quantile normal (Q-Q) plot
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=1)
abline(0, 0, lty=2)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids>=-S & sresids <= S], predict(sresidsMod)), col=1, lty=1)


detach(sfData_S_13hHyp)
# ----------------------------------------------------------------------#
# Mitosis; 0 hours (hypoxic)
attach(sfData_Mitose_0hHyp)

newDose = rep(Dose.Gy., 2)
newF = c(F3823,F3824)

plot(Dose.Gy., F3823, main="Mitosis; 0 hours (hypoxic)", 
     ylab="Survival fraction", xlab="Dose [Gy]", ylim = c(min(newF), max(newF)),
     las=1, col=2, pch=1, log="y")
points(Dose.Gy., F3824, col=2, pch=1)

mod = lm(log(newF) ~ newDose + I(newDose^2))
lines(smooth.spline(Dose.Gy., I(exp(predict(mod)[1:5]))), col=1, lty=1)
prd = predict(mod, level=0.95, interval = "confidence")
lines(sort(Dose.Gy.), exp(prd[1:5, 2][order(Dose.Gy.)]), col=2, lty=2)
lines(sort(Dose.Gy.), exp(prd[1:5, 3][order(Dose.Gy.)]), col=2, lty=2)
legend("bottomleft", legend=c("Regression line", "Confidence interval (95%)", 
                              paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
                              paste("R-Sq =", round(summary(mod)$r.squared, 4)) ), 
       col=c(1,2, 0, 0), lty=c(1,2,0, 0), bty="o")

confint(mod, conf.level=0.95) # confidence intervall for model coefficients
coefficients(mod)[2]  # alpha
coefficients(mod)[3]  # beta
summary(mod)$r.squared
summary(mod)$adj.r.squared
summary(mod)

modTarget = nls(log(newF) ~ log(1 - (1 - exp(-v*newDose) )^m ), start = list(v=1, m=1))
#modTarget = nls(log(newF) ~ -v*newDose + log(m), start = list(v=1, m=1)) # approx for high dose D
lines(smooth.spline(newDose, predict(modTarget)), col=3, lty=1)

# Diagnostics

# Quantile normal (Q-Q) plot for LQ
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Quantile normal (Q-Q) plot for multi-target
sresidsTarget = resid(modTarget)
qqnorm(sresidsTarget, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresidsTarget, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=2)
points(fitted(modTarget), sresidsTarget, col=3, pch=3)
abline(0, 0, lty=1)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
sresidsModTarget = lm(sresidsTarget[sresidsTarget >= -S & sresidsTarget <= S] ~ 
                        fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S] + I(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids >= -S & sresids <= S], predict(sresidsMod)), col=2, lty=2)
lines(smooth.spline(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S], predict(sresidsModTarget)), col=3, lty=3)
legend("topright", legend=c("LQ model", "Multi-target model"), col=c(2,3), lty=c(2,3), pch=c(2,3), bty="o")

detach(sfData_Mitose_0hHyp)

# ----------------------------------------------------------------------#
# Exponential growth; asynchronous (hypoxic)
attach(sfData_expAsynchHyp)

newDose = rep(Dose.Gy., 3)
newF = c(F3839,F3840, F3842)

plot(Dose.Gy., F3839, main="Exponential growth; asynchronous (hypoxic)", 
     ylab="Survival fraction", xlab="Dose [Gy]", ylim = c(min(newF), max(newF)),
     las=1, col=2, pch=1, log="y")
points(Dose.Gy., F3840, col=2, pch=1)
points(Dose.Gy., F3842, col=2, pch=1)

mod = lm(log(newF) ~ poly(newDose, 2, raw=TRUE))
mod = lm(log(newF) ~ newDose + I(newDose^2))
lines(smooth.spline(Dose.Gy., I(exp(predict(mod)[1:7]))), col=1, lty=1)
prd = predict(mod, level=0.95, interval = "confidence")
lines(sort(Dose.Gy.), exp(prd[1:7, 2][order(Dose.Gy.)]), col=2, lty=2)
lines(sort(Dose.Gy.), exp(prd[1:7, 3][order(Dose.Gy.)]), col=2, lty=2)
legend("bottomleft", legend=c("Regression line", "Confidence interval (95%)", 
                              paste("R-Sq(adj) =", round(summary(mod)$adj.r.squared, 4)), 
                              paste("R-Sq =", round(summary(mod)$r.squared, 4)) ), 
       col=c(1,2, 0, 0), lty=c(1,2,0, 0), bty="o")

confint(mod, conf.level=0.95) # confidence intervall for model coefficients
coefficients(mod)[2]  # alpha
coefficients(mod)[3]  # beta
summary(mod)$r.squared
summary(mod)$adj.r.squared
summary(mod)

modTarget = nls(log(newF) ~ log(1 - (1 - exp(-v*newDose) )^m ), start = list(v=1, m=1))
#modTarget = nls(log(newF) ~ -v*newDose + log(m), start = list(v=1, m=1)) # approx for high dose D
lines(smooth.spline(newDose, predict(modTarget)), col=3, lty=1)

# Diagnostics

# Quantile normal (Q-Q) plot for LQ
sresids = resid(mod)
qqnorm(sresids, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresids, col=1)

# Quantile normal (Q-Q) plot for multi-target
sresidsTarget = resid(modTarget)
qqnorm(sresidsTarget, las=1, ylab="Standardized residuals e*", col=2)
qqline(sresidsTarget, col=1)

# Residual plot
matplot(fitted(mod), sresids, ylab="Standardized residuals e*", xlab="Estimated survival fraction (fitted values)", col=2, pch=2)
points(fitted(modTarget), sresidsTarget, col=3, pch=3)
abline(0, 0, lty=1)
S = sqrt(deviance(mod)/df.residual(mod))
sresidsMod = lm(sresids[sresids >= -S & sresids <= S] ~ 
                  fitted(mod)[sresids>=-S & sresids <= S] + I(fitted(mod)[sresids>=-S & sresids <= S]^2) )
sresidsModTarget = lm(sresidsTarget[sresidsTarget >= -S & sresidsTarget <= S] ~ 
                        fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S] + I(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S]^2) )
lines(smooth.spline(fitted(mod)[sresids >= -S & sresids <= S], predict(sresidsMod)), col=2, lty=2)
lines(smooth.spline(fitted(modTarget)[sresidsTarget >=-S & sresidsTarget <= S], predict(sresidsModTarget)), col=3, lty=3)
legend("topright", legend=c("LQ model", "Multi-target model"), col=c(2,3), lty=c(2,3), pch=c(2,3), bty="o")

detach(sfData_expAsynchHyp)

# Assumptions made when fitting in a linear regression model
# 1. The Y-values (or the errors, "e") are independent
# 2. The Y-values can be expressed as a linear function of the X variable
# 3. Variation of observations around the regression line (the residual SE)
#     is constant (homoscedasticity)
# 4. For given value of X, Y values (or the error) are normally distributed
# Assumption 2-4 can be checked by examining the residuals or error