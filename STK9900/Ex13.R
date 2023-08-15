# Read the data into a dataframe, give names to the variables, and inspect 
# the data:
salinity = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/exer5.dat")
names(salinity) = c("salt", "saltprev", "trend", "discharge")
salinity

# Get an overview of the data:
summary(salinity)
plot(salinity)
# Make sure that you understand what the summary measures tell you!
# What do you see from the scatter plots?

# Do linear regression with all three covariates and inspect the results:
lmfull = lm(salt ~ saltprev + trend + discharge, data = salinity)
summary(lmfull)
# How important are each of the covariates? How does this agree with the 
# scatter plots?

# Extract fitted values and residuals, and give a summary of the residuals:
saltfit = lmfull$fit
saltres = lmfull$res
summary(saltres)

# (i) Checking normality
# Histogram and Q-Q plot (make one plot at a time)
boxplot(saltres)
hist(saltres)
qqnorm(saltres); qqline(saltres)    # alternative command: plot(lmfull,2)
# What do the plots tell you?

# (ii) Checking homoscedasticy
# Plot of residuals versus fitted values:
plot(saltfit, saltres, xlab = "Fitted values", ylab = "Residuals")
# What does the plot tell you?
# Alternative plots:
plot(lmfull,1)
plot(lmfull,3)

# (iii) Checking linearity
# We will make CPR plots for each of the predictors
# We then use the "car" library (package)
library(car)
crPlots(lmfull, terms=~.)
# What do the plots tell you? Are there indications of deviation from linearity?
# (Note that the green line in the plots may be influenced quite a lot by a 
# few observations, so one should be careful not to overinterpret the plots )  

# (iv) Studying influence
# We will make boxplots of the standardized dfbetas to check if some 
# observations have a large influence of the estimates of the fitted model:
boxplot(dfbetas(lmfull)[,-1])

# You will find that one observation is very influential for the estimated 
# effect of "discharge"
# A simple way of identifying this observations is to just list the 
# observations together with the standardized defbetas for "discharge" 
# (which is the 4th parameter when we also count the intercept)
cbind(salinity,dfbetas(lmfull)[,4])
# Which observation is the influential one?

# A more "advanced" way to identify the influential observation, is to use the 
# identify-command (cf. the slides from the lectures).

# If time allows you may try the following commands for identifying the 
# observation that is influential for the estimate of "discharge"
db = dfbeta(lmfull)
j = 4
boxplot(db[,j])
identify(rep(1, dim(db)[1]), db[,j], labels = 1:dim(db)[1])

# Try fitting a model with a second degree term for discharge
lmfull2 = lm(salt ~ saltprev + trend + discharge + I(discharge^2), 
             data = salinity)
summary(lmfull2)

