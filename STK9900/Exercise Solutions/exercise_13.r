
# STK 4900
# Exercise 13

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# a)
# Read data.
salinity.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/data/exer5.dat", header=FALSE)
colnames(salinity.data)=c("salt","saltprev","trend","discharge")
# Take a look at the data. 
salinity.data

# Get an overview of the data:
summary(salinity.data)

# Create a scatter plot matrix.
plot(salinity.data)


# Make sure that you understand what the summary measures tell you!
# What do you see from the scatter plots?

# Fit multiple linear regression with all three covariates.
lm.obj = lm(salt~saltprev+trend+discharge, data=salinity.data)
summary(lm.obj)
# How important are each of the covariates? How does this agree with the scatter plots?


# b)
# Extract fitted values and residuals, and give a summary of the residuals:
y.hat = lm.obj$fit
residual.val = lm.obj$res

# Check
salinity.data[,"salt"] - (y.hat + residual.val)
# Why are these 0?


# c)
# We will make various plots of the residuals

# (i) Checking normality
par(mfrow=c(1,2))
# Histogram and Q-Q plot (make one plot at a time)
hist(residual.val)
# Q-Q plot
qqnorm(residual.val)
qqline(residual.val)

par(mfrow=c(1,1))


# What do the plots tell you?


# (ii) Checking homoscedasticy
# Plot of residuals versus fitted values:
plot(y.hat, residual.val, xlab="Fitted values", ylab="Residuals")


# What does the plot tell you?


# (iii) Checking linearity
# We will make CPR plots for each of the predictors
# We use the "car" library (package)
library(car)
crPlots(lm.obj, terms=~.)


# What do the plots tell you? Are there indications of deviation from linearity?
# (Note that the green line in the plots may be influenced quite a lot by a few observations,
# so one should be careful not to overinterpret the plots)

# (iv) Studying influence
# We will make boxplots of the standardized dfbetas to check if some
# observations have a large influence of the estimates of the fitted model.

boxplot(dfbetas(lm.obj)[,-1])


# You will find that one observation is very influential for the estimated effect of "discharge".
# A simple way of identifying this observations is to just list
# the observations together with the standardized defbetas
# for "discharge" (which is the 4th parameter when we also count the intercept)
cbind(salinity.data, dfbetas(lm.obj)[,4])
# Which observation is the influential one?

# A more "advanced" way to identify the influential observation, is to use
# the identify-command (cf. the slides from the lectures).
# If time allows you may try the following commands for identifying the observation
# that is influential for the estimate of "discharge"
#db = dfbeta(lm.obj)
#j = 4
#boxplot(db[,j])
#identify(rep(1,dim(db)[1]),db[,j], labels=1:dim(db)[1])

# Or we can sort the dataframes with values in decreasing order!


# d)
# Try fitting a model with a second degree term for discharge
lm.obj.2 = lm(salt~saltprev+trend+discharge+I(discharge^2), data=salinity.data)
summary(lm.obj.2)
# Check the model assumptions for this model. Is it improved comare to the model without a second degree term?

crPlots(lm.obj.2)

# Try yourself other models (e.g. with interactions and/or more second order terms)
# Which model would you suggest to use for predicting salinity?





