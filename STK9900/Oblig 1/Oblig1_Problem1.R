# Clean up the memory before start
rm(list = ls(all = TRUE))

# Read the data into R by the command:
no2 = read.table(
  "https://www.uio.no/studier/emner/matnat/math/STK4900/v21/obliger/no2.txt", 
  header = T)

### Problem 1a) ###

# Compute mean and median of the variables log.no2 and log.cars: 
summary(no2[c("log.no2", "log.cars")])

# Make a boxplot of the variables log.no2 and log.cars: 
boxplot(no2[c("log.no2", "log.cars")], xlab = "Variable", ylab = "Value") 

# Make a scatter plot of the log(no2) versus log(cars):
plot(no2$log.cars, no2$log.no2, pch = 19)

# Problem 1b)

# Fit a linear regression model to the data, using log concentration of NO2 as 
# outcome and amount of traffic measured by log(number of cars per hour) as 
# covariate:
fit1 = lm(log.no2 ~ log.cars, data = no2, x = T, y = T)
summary(fit1)
abline(fit1)

### Problem 1c) ###

# Extract fitted values and residuals:
no2fit1 = fit1$fit
no2res1 = fit1$res

# (i) Checking linearity
# Make CPR plots for each of the predictors (require the "car" package)
library(car)
crPlots(fit1, terms=~.)

# (ii) Checking homoscedasticity
# Plot of residuals versus fitted values:
plot(fit1,1)
plot(fit1,3)

# (iii) Checking normality
# Histogram and Q-Q plot (make one plot at a time)
boxplot(no2res1)
hist(no2res1)
plot(fit1,2)

### Problem 1d) ###

# Compute the correlation between the four covariates and make scatterplots 
# of the covariates:
cor(no2[,2:5])
plot(no2[,2:5])

# Fit models according to this ordering where the temperature is measured in 
# Kelvin (in order to log transform the temperature):
fit2 = lm(log.no2 ~ log.cars + I(temp+273.15) + wind.speed + hour.of.day, 
          data = no2, x = T, y = T)
fit3 = lm(log.no2 ~ log.cars + log(I(temp+273.15)) + wind.speed + hour.of.day, 
          data = no2, x = T, y = T)
fit4 = lm(log.no2 ~ log.cars + I(temp+273.15) + log(wind.speed) + hour.of.day, 
          data = no2, x = T, y = T)
fit5 = lm(log.no2 ~ log.cars + I(temp+273.15) + wind.speed + log(hour.of.day), 
          data = no2, x = T, y = T)

summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5)

# Compute adjusted R-squared cross-validated R-squared for the models:
summary(fit1)$adj.r.squared
summary(fit2)$adj.r.squared
summary(fit3)$adj.r.squared
summary(fit4)$adj.r.squared
summary(fit5)$adj.r.squared

cv.R2(fit1)
cv.R2(fit2)
cv.R2(fit3)
cv.R2(fit4)
cv.R2(fit5)

### Problem 1e) ###

# Extract fitted values and residuals:
no2fit4 = fit4$fit
no2res4 = fit4$res

# (i) Checking linearity
# Make CPR plots for each of the predictors (require the "car" package)
crPlots(fit4, terms=~.)

# (ii) Checking homoscedasticity
# Plot of residuals versus fitted values:
plot(fit4,1)
plot(fit4,3)

# (iii) Checking normality
# Histogram and Q-Q plot (make one plot at a time)
boxplot(no2res4)
hist(no2res4)
plot(fit4,2)

cv.R2=function(lmobj,y=lmobj$y,x=lmobj$x) {
  
  a = t(x)%*%x
  d = diag(1/eigen(a)$values)
  e = eigen(a)$vector
  inv.a = e %*% d %*% t(e)
  v = diag(x%*%inv.a%*%t(x))
  SSkryss = sum((lmobj$res/(1-v))^2)
  SStot = sum((y-mean(y))^2)
  R2kryss = 1-SSkryss/SStot
  R2kryss
  
}
