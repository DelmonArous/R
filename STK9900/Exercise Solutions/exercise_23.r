
# STK 4900
# Exercise 23

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# The data are from a old study of children with leukemia.
# The children were treated so that they had no symptoms (one say that they are in remission).
# They were then given an active treatment ("Drug") to keep the symptoms from coming back or placebo ("Control").
# Time until the symptoms came back (called relapse) were measured.
# For some children the symptoms had not come back at the end of the study. These give rise to censored observations.

# Read data.
leukemia.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/data/gehan.txt", header=TRUE)
# Take a look at the data.
leukemia.data
# "time" is time in weeks to relapse or censoring
# "cens" is 1 for relapse and 0 for censoring
# "treat" is 1 for "Control" and 2 for "Drug"

# Load the R-library for survival analysis.
library(survival)


# QUESTION 1
# Compute Kaplan-Meier estimates for the two groups (without confidence intervals)
surv.obj.1 = survfit(Surv(time,cens)~treat, conf.type="none", data=leukemia.data)
summary(surv.obj.1)
# Make sure you understand what the output tells!
# See Slide 12, Lecture 9.

# Plot the Kaplan-Meier estimates.
plot(surv.obj.1, lty=1:2)
abline(h=0.5, col="red", lty=2)
# Interpret the plots.
# Read from the plots (approximately) what is the median time to relapse for the two groups
# Check the results by giving the command print(surv.obj.1)


# QUESTION 2
# Compute Kaplan-Meier estimates for the two groups with confidence intervals
# (the default confidence interval in R is not a good choice, so we make an explicit choice of the type of confidence interval)
surv.obj.2 = survfit(Surv(time,cens)~treat, conf.type="plain", data=leukemia.data)
summary(surv.obj.2)

plot(surv.obj.2, conf.int=TRUE, lty=1:2)
abline(h=0.5, col="red", lty=2)
# Interpret the output and the plot.


# QUESTION 3
# Log-rank test for difference between the groups:
survdiff(Surv(time,cens)~treat, data=leukemia.data)
# What does the output tell you?


# QUESTION 4
# Cox regresion with treatment group as covariate:
cox.obj = coxph(Surv(time,cens)~factor(treat), data=leukemia.data)
summary(cox.obj)
# Interpret the results!
