
# STK 4900
# Exercise 18

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# In this exercise we will use data from the Western Collaborative Group Study (WCGS), a large epidemiological
# study designed to study risk factors for coronary heart disease (CHD).
# More than 3000 men were recruited to the study, and a number of (potential) risk factors were recorded at entry.
# The men were then followed for about ten years and it was recorded if they developed CHD or not.
# The WCGS data are used in many of the examples in Section 3.4 and Chapter 6 of the course text book byVittinghoff et al.

# Read data.
WCGS.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/data/wcgs.txt", header=TRUE, sep="\t", na.strings=".")
# In this exercise we will restrict ourselves to study the effect of smoking and age on the risk for CHD.

# a)
# We first fit a logistic regression model with smoke as covariate.
glm.obj.1 = glm(chd69~smoke, data=WCGS.data, family=binomial)
summary(glm.obj.1)
# Is there a significant effect of smoking on the risk of developing CHD? (Cf. slide 11)

# Compute a 95% confidence interval for the regression coefficient for smoking (cf slide 9)
# Also estimate the odds ratio for smoking and determine a 95% confidence interval for the odds ratio (cf slide 9).

# b)
# In order to estimate the odds ratio with 95% confidence limits from a fitted logistic model, you may use the following function (cf slide 10).
exp.coef.func = function(glm.obj) {
  
  alpha = 0.05
  coef.mat = summary(glm.obj)$coef
  
  lower = exp(coef.mat[,1] - qnorm(p=1-alpha/2)*coef.mat[,2])
  upper = exp(coef.mat[,1] + qnorm(p=1-alpha/2)*coef.mat[,2])
  result = cbind(estimate=exp(coef.mat[,1]), lower, upper)
  
  return(result)
}
# Use this function to estimate the odds ratio for smoking with 95% confidence limits.
# Check that you get the same results as in question a).
exp.coef.func(glm.obj.1)

# We can also use the R's built-in solution.
exp(confint.default(glm.obj.1))


# c)
# We then use logistic regression to study the effect of age (at entry to the study) for the risk of developing CHD.
glm.obj.2 = glm(chd69~age, data=WCGS.data, family=binomial)
summary(glm.obj.2)
# Use the "exp.coef.func()" to estimate the odds ratio for a one-year increase in age with 95% confidence interval.
exp.coef.func(glm.obj.2)
# Alternative
exp(confint.default(glm.obj.2))

# d)
# When interpreting the effect of age, it may be reasonable to give the odds ratio corresponding to
# a ten-year increase in age (rather than a one-year increase as we did in question c).
# The easiest way to achieve this is to fit the logistic model using age/10 as covariate.
glm.obj.3 = glm(chd69~I(age/10), data=WCGS.data, family=binomial)
summary(glm.obj.3)

# Compare the estimates and standard errors of this model with the one in question c). What do you see?
# You may then use the "exp.coef.func()" to estimate the odds ratio for a ten-year increase in age with 95% confidence interval. (Why?)
exp(confint.default(glm.obj.3))

