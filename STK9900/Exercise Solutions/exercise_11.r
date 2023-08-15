
# STK 4900
# Exercise 11

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# Check and change working directory.
#getwd()
#setwd("C:/Users/username/Documents")

# In this exercise we will use data from the Heart and Estrogen/Progestin Study (HERS), a clinical trial of
# hormone therapy for prevention of recurrent heart attacks and death among post-menopausal women with existing coronary heart disease.

# The aim of this exercise is to study how the change in low-density lipoprotein (LDL) cholesterol over the first year of
# the HERS study depends on the baseline value of LDL (i.e. the value at entry to the study) and use of hormone therapy (HT).

# Read data
HERS.data = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/hers.txt", header=TRUE, sep="\t", na.strings=".")

# Before we start doing our analysis we define the change in LDL, denoted LDLch.
HERS.data[,"LDLch"] = HERS.data[,"LDL1"] - HERS.data[,"LDL"]
#We also defined the centered LDL at baseline (by subtracting the mean value 145 mg/dL), denoted cLDL
HERS.data[,"cLDL"] = HERS.data[,"LDL"] - 145


# a)
# Fit a linear model with the change in LDL as the response and hormone therapy (HT) and baseline LDL (not centered) as covariates:
lm.obj.1 = lm(LDLch~HT+LDL, data=HERS.data)
summary(lm.obj.1)
# Give an interpretation of the estimates in the fitted model.



# b)
# We then fit a model with HT and centered LDL at baseline as covariates:
lm.obj.2 = lm(LDLch~HT+cLDL, data=HERS.data)
summary(lm.obj.2)
# Compare the estimates for this model with those in question a).
# What is the interpretation of the estimate for the intercept?

# The package 'scales' permits the use of the function 'alpha' which controls the transparency of the colours in R.
# It is relly useful trick for crowded scatterplots!
par(mfrow=c(1,2))
# HT==0
plot(LDLch~cLDL, data=HERS.data[HERS.data$HT==0, ], main="Control", cex=.6, col=scales::alpha("darkgrey", .5), xlim=c(-100, 200), ylim=c(-200, 150))
abline(a=lm.obj.2$coefficients[1], b=lm.obj.2$coefficients[3], col="green3", lwd=2)
# HT==1
plot(LDLch~cLDL, data=HERS.data[HERS.data$HT==1, ], main="Treatment", cex=.6, col=scales::alpha("darkgrey", .5), xlim=c(-100, 200), ylim=c(-200, 150))
abline(a=lm.obj.2$coefficients[1]+lm.obj.2$coefficients[2], b=lm.obj.2$coefficients[3], col="darkorange", lwd=2)
par(mfrow=c(1,1))



# c)
# We then fit a model with interaction.
lm.obj.3 = lm(LDLch~HT+cLDL+I(HT*cLDL), data=HERS.data)
summary(lm.obj.3)
# Is there a significant interaction?
# Given an interpretation of the estimates in the model (cf. slide 19 of Lecture 4)

# In particular, what is the effect of baseline LDL for those with no hormone therapy?
# And what is the effect of baseline LDL for those on hormone therapy?

# Let's compare the coefficients

lm.obj.2$coefficients
lm.obj.3$coefficients

par(mfrow=c(1,2))
# HT==0
plot(LDLch~cLDL, data=HERS.data[HERS.data$HT==0, ], main="Control", cex=.6, col=scales::alpha("darkgrey", .5), xlim=c(-100, 200), ylim=c(-200, 150))
abline(a=lm.obj.3$coefficients[1], b=lm.obj.3$coefficients[3], col="green3", lwd=2)
# HT==1
plot(LDLch~cLDL, data=HERS.data[HERS.data$HT==1, ], main="Treatment", cex=.6, col=scales::alpha("darkgrey", .5), xlim=c(-100, 200), ylim=c(-200, 150))
abline(a=lm.obj.3$coefficients[1]+lm.obj.3$coefficients[2], b=lm.obj.3$coefficients[3]+lm.obj.3$coefficients[4], col="darkorange", lwd=2)
par(mfrow=c(1,1))



#### Bonus: confidence vs. prediction intervals ####

# We start by fittin a simple model
lmSimple = lm(LDLch~cLDL, data=HERS.data)
summary(lmSimple)

# Now we prepare data for "prediction"
newdata = data.frame(cLDL=seq(min(HERS.data$cLDL, na.rm=T)-1, max(HERS.data$cLDL, na.rm=T)+1, length.out=1000))
confInts = predict(lmSimple, newdata=newdata, interval="confidence")
predInts = predict(lmSimple, newdata=newdata, interval="prediction")

# Let's plot both intervals and see what happens.
plot(HERS.data$cLDL,
     HERS.data$LDLch,
     cex=.8,
     pch=16,
     col=scales::alpha("darkgrey", .5),
     xlab="cLDL",
     ylab="LDLch")
abline(a=lmSimple$coefficients[1], b=lmSimple$coefficients[2], col="red", lwd=2)

lines(newdata$cLDL, confInts[,"lwr"], lty=3, col="blue", lwd=1.5)
lines(newdata$cLDL, confInts[,"upr"], lty=3, col="blue", lwd=1.5)

lines(newdata$cLDL, predInts[,"lwr"], lty=3, col="darkorange", lwd=1.5)
lines(newdata$cLDL, predInts[,"upr"], lty=3, col="darkorange", lwd=1.5)

