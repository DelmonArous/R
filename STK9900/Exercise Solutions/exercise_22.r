
# STK 4900
# Exercise 22

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# Read data.
cancer.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/data/cancer.txt", header=FALSE)
names(cancer.data)=c("age","cig","pyr","cancer")
# Take a look at the data.
head(cancer.data)
tail(cancer.data)

# Questions 1 & 2)
# We first consider the model  E(Y) = n*exp{b0+b1*s+b2*a}, where
#     Y = number of cancer cases (= cancer),
#     n = number of person years (= pyr),
#     s = number of cigarettes smoked per day (= cig)
#     a = age in years (= age)
# We may write the model on the form  E(Y)= exp{1*log(n)+b0+b1*s+b2*a}.

# Note that log(n) appears as a sort of "covariate" where we know that the regression coefficient takes the value 1.
# This is called an OFFSET.

# We fit the model and look at the result.
glm.obj.1 = glm(cancer~offset(log(pyr))+age+cig, data=cancer.data, family=poisson)
summary(glm.obj.1)
# Make sure that you understand what the output tells you!
# Are there significant effects of age and the number of cigarettes smoked?

# It is common to report the results of a Poisson regression by means of the rate ratio RR = exp(beta) with confidence limits.
# To this end we may use the function "exp.coef.func()" from Exercise 18.
# Use the function to compute rate ratios for age and number of cigarettes:
exp.coef.func(glm.obj.1)

# Give an interpretation of what the table tells you about the effect of age and the number of cigarettes smoked


# QUESTION 3)
# We then look at a model with second order terms and interaction:
glm.obj.2 = glm(cancer~offset(log(pyr))+age+I(age^2)+cig+I(cig^2)+age:cig, data=cancer.data, family=poisson)
# Reduce the model by (step-wise) eliminating non-significant covariates.
# (Use Wald tests from the summary-command and/or deviances from the anova-command.)
# Discuss the interpretation of your "final model".

# ADDITIONAL QUESTION:
# Age and the number of cigarettes smoked are reported in intervals.
# We may alternatively consider these covariates as categorical.
# Such a model is fitted by the command:

glm.obj.3 = glm(cancer~offset(log(pyr))+factor(age)+factor(cig), data=cancer.data, family=poisson)
summary(glm.obj.3)
# Give an interpretation of this model.
# Discuss how the model may be used to assess the fit of your "final model" from question 3.

anova(glm.obj.3, glm.obj.2, test="Chisq")
