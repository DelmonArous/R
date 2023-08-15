rm(list = ls(all = TRUE))

# Read the data into a dataframe, give names to the variables, 
# and inspect the data:
insurance = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/exer3_2.dat")
names(insurance) = c("income", "riskave", "amount")
insurance

# Attach the dataframe:
attach(insurance)

# Compute summary measures for the variables:
summary(insurance)
    
# Make plots (side by side) of amount versus each of the other two variables:
par(mfrow = c(1,2))
plot(income, amount)
plot(riskave, amount)
par(mfrow = c(1,1))

# Compute the correlation between the variables:
cor(insurance)
# Do the correlations agree with what you saw from the plots?

# Do univariate regression analyses of amount versus each of the other 
# two variables:
fit1 = lm(amount ~ income)
fit2 = lm(amount ~ riskave)
summary(fit1)
summary(fit2)
# Which of the two variables, income and risk aversion, is most important for 
# explaining the variation in the amount of life insurance carried?
# Does any of the variables (alone) have a significant effect?

# Do a regression analysis including both income and risk aversion:
fit3 = lm(amount ~ income + riskave)
summary(fit3)
# What does this model tell you? Does it look better than the best of the 
# two models with only one covariate?

