
# STK 4900
# Week 1, Exercise 9

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# Read the data into a dataframe and give names to the variables.
insurance.data = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/exer3_2.dat", header=FALSE)
names(insurance.data)=c("income","risk.avs","amount")

# Take a look at the data.
insurance.data
summary(insurance.data)
# Make sure that you understand what the summary measures tell you!

# Create scatter plots between 'amount' and other variables.

# 1 by 2 grid of plots.
par(mfrow=c(1,2))
# Scatter plot: 'amount' vs 'income'
plot(x=insurance.data[,"income"],
     y=insurance.data[,"amount"],
     xlab="Annual income",
     ylab="Amount of life insurance",
     main = "")
# Scatter plot: 'amount' vs 'risk.avs'
plot(x=insurance.data[,"risk.avs"],
     y=insurance.data[,"amount"],
     xlab="Risk aversion score",
     ylab="Amount of life insurance",
     main = "")
par(mfrow=c(1,1))

# Alternative: scatter plot matrix
plot(insurance.data)


# What do the plots tell you?



# Compute the correlation between the variables.
cor(insurance.data)
# Do the correlations agree with what you saw from the plots?

# Carry out univariate regression analysis of 'amount' versus each of the other two variables.
lm.obj.1 = lm(amount~income, data=insurance.data)
summary(lm.obj.1)
lm.obj.2 = lm(amount~risk.avs, data=insurance.data)
summary(lm.obj.2)

# Which of the two variables, income and risk aversion, is most important for explaining the variation in the amount of life insurance carried?
# Does any of the variables (alone) have a significant effect?

# Carry out a regression analysis including both income and risk aversion:
lm.obj.3 = lm(amount~income+risk.avs, data=insurance.data)
summary(lm.obj.3)

# What does this model tell you? Does it look better than the best of the two models with only one covariate?
# Try yourself models with second order terms for income and/or risk aversion.
lm.obj.1.poly = lm(amount~income+I(income^2), data=insurance.data)
summary(lm.obj.1.poly)

lm.obj.2.poly = lm(amount~risk.avs+I(risk.avs^2), data=insurance.data)
summary(lm.obj.2.poly)

x.grid = seq(min(insurance.data$risk.avs)-1, max(insurance.data$risk.avs)+1, length.out=1000)
y.predicted = predict.lm(lm.obj.2.poly, newdata=data.frame(risk.avs=x.grid))

par(mfrow=c(1,2))
# Scatter plot: 'amount' vs 'income'
plot(x=insurance.data[,"income"],
     y=insurance.data[,"amount"],
     xlab="Annual income",
     ylab="Amount of life insurance",
     main = "")
abline(lm.obj.1, col="red")
# Scatter plot: 'amount' vs 'risk.avs'
plot(x=insurance.data[,"risk.avs"],
     y=insurance.data[,"amount"],
     xlab="Risk aversion score",
     ylab="Amount of life insurance",
     main = "")
lines(x=x.grid,
      y=y.predicted,
      col="blue")
par(mfrow=c(1,1))
