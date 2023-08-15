rm(list = ls(all = TRUE))

# Read the data into R by the command:
pef     = c(494, 395, 516, 434, 476, 413, 442, 433)
minipef = c(512, 430, 520, 428, 500, 364, 380, 445)

# Make scatter plot
plot(pef, minipef, pch = 19)  

# Compute the Pearson correlation coefficient
cov(pef, minipef)

cov(pef, minipef)/(sd(pef)*sd(minipef))
cor(pef, minipef) # check

# 95% confidence interval for the (theoretical) correlation coefficient
cor.test(pef, minipef)

# Fit a linear regression model with "minipef" as the outcome and 
# "pef" as the predictor
fit     = lm(minipef ~ pef)
output  = summary(fit)
abline(fit)

#  Multiply the estimated slope from question d with the ratio of the 
# empirical standard deviations of "pef" and "minipef"
beta_1  = output$coefficients[2,1]
s_x     = output$coefficients[2,2]
s_y     = output$coefficients[1,2]

rhat = beta_1*(s_x/s_y)


