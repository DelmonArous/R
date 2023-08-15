# Clean up the memory before start
rm(list = ls(all = TRUE))

# Read the olympic data into R by the command:
cirrhosis = read.table( 
  "https://www.uio.no/studier/emner/matnat/math/STK4900/data/cirrhosis.txt", 
  sep = "\t", header = T, na.strings = ".")

# Attach the R-library for survival analysis:
library(survival)


### Problem 3a) ###

par(mfrow=c(2,2))

# Compute Kaplan-Meier estimates for the two treatment groups (prednisone, 
# placebo) (without confidence intervals)
survpred.treat = survfit(Surv(time,status) ~ treat, data = cirrhosis, 
                   conf.type = "plain")
print(survpred.treat)
plot(survpred.treat, lty = 1:2, main = "Survival of the treatment covariate", 
     xlab = "Time (days)", ylab = "Survival")
legend("topright", c("Prednisone","Placebo"), lty = 1:2)

# Compute Kaplan-Meier estimates for the two gender groups (females, males) 
# (without confidence intervals)
survpred.sex = survfit(Surv(time,status) ~ sex, data = cirrhosis, 
                   conf.type = "plain")
print(survpred.sex)
plot(survpred.sex, lty = 1:2, main = "Survival of the sex covariate", 
     xlab = "Time (days)", ylab = "Survival")
legend("topright", c("Female","Male"), lty = 1:2)

# Compute Kaplan-Meier estimates for the three ascites groups (none, slight, 
# marked ascites at start of treatment) (without confidence intervals)
survpred.asc = survfit(Surv(time,status) ~ asc, data = cirrhosis, 
                       conf.type = "plain")
print(survpred.asc)
plot(survpred.asc, lty = 1:3, main = "Survival of the ascites covariate", 
     xlab = "Time (days)", ylab = "Survival")
legend("topright", c("None ascites","Slight ascites","Marked ascites"), lty = 1:3)

# Compute Kaplan-Meier estimates for the three age groups (<50 years, 50-65 
# years, >65 years) (without confidence intervals)
survpred.agegr = survfit(Surv(time,status) ~ agegr, data = cirrhosis, 
                       conf.type = "plain")
print(survpred.agegr)
plot(survpred.agegr, lty = 1:3, main = "Survival of the grouped age covariate", 
     xlab = "Time (days)", ylab = "Survival")
legend("topright", c("<50 years","50-65 years",">65 years"), lty = 1:3)

# Print out the median survival time for each covariate
print(survpred.treat)
print(survpred.sex)
print(survpred.asc)
print(survpred.agegr)


### Problem 3b) ###

# Log-rank test for difference between the two treatment groups (prednisone, 
# placebo):
survdiff(Surv(time,status) ~ treat, data = cirrhosis)

# Log-rank test for difference between the two gender groups (females, males):
survdiff(Surv(time,status) ~ sex, data = cirrhosis)

# Log-rank test for difference between the three ascites groups (none, slight, 
# marked ascites at start of treatment):
survdiff(Surv(time,status) ~ asc, data = cirrhosis)

# Log-rank test for difference between the three age groups (<50 years, 50-65 
# years, >65 years):
survdiff(Surv(time,status) ~ agegr, data = cirrhosis)


### Problem 3c) ###

# Cox regresion with treatment, sex, ascites severity and age as covariates:
fit.full = coxph(Surv(time,status) ~ factor(treat) + factor(sex) + factor(asc) 
                 + I(age - 17), data = cirrhosis)
summary(fit.full)
