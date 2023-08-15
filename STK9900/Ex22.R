# Clean up the memory before start
rm(list=ls(all=TRUE))

# Read the data into a dataframe, give names to the variables, 
# and inspect the data:
cancerdata = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/cancer.txt")
names(cancerdata) = c("age", "cig", "pyr", "cancer")
cancerdata

# Make sure that the data are the same as given in the exercise.

# Questions 1 & 2)
# We first consider the model  E(Y) = n*exp{b0+b1*s+b2*a}, where
#       Y = number of cancer cases (= cancer),
#       n = number of person years (= pyr),
#       s = number of cigarettes smoked per day (= cig)
#       a = age in years (= age)
# We may write the model on the form  E(Y) = exp{1*log(n) + b0 + b1*s + b2*a}.
# Note that  log(n)  appears as a sort of "covariate" where we know that the 
# regression coefficient takes the value 1. This is called an OFFSET.


# We fit the model and look at the result:
cancerfit.1 = glm(cancer ~ offset(log(pyr)) + age + cig, data = cancerdata, 
                  family = poisson)
summary(cancerfit.1)

# Make sure that you understand what the output tells you.
# Are there significant effects of age and the number of cigarettes smoked?


# It is common to report the results of a Poisson regression by means of the 
# rate ratio RR = exp(beta) with confidence limits.  
# To this end we may use the function expcoef from exercise 18.
# Use the function to compute rate ratios for age and number of cigarettes:
expcoef(cancerfit.1)

# Give an interpretation of what the table tells you about the effect of age 
# and the number of cigarettes smoked


# Question 3)
# We then look at a model with second order terms and interaction:
cancerfit.2 = glm(cancer ~ offset(log(pyr)) + age+I(age^2) + cig + I(cig^2) + 
                    age:cig, data = cancerdata, family = poisson)
summary(cancerfit.2)
# Reduce the model by (step-wise) eliminating non-significant covariates.
# (Use Wald tests from the summary-command and/or deviances from the 
# anova-command). Discuss the interpretation of your "final model".

anova(cancerfit.1, cancerfit.2, test = "Chisq")

# ADDITIONAL QUESTION:
# Age and the number of cigarettes smoked are reported in intervals.
# We may alternatively consider these covariates as categorical.
# Such a model is fitted by the command:
cancerfit.a = glm(cancer ~ offset(log(pyr)) + factor(age) + factor(cig), 
                  data = cancerdata, family = poisson)
summary(cancerfit.a)
# Give an interpretation of this model.
# Discuss how the model may be used to assess the fit of your "final model" 
# from question 3.

anova(cancerfit.a, cancerfit.2, test = "Chisq")


expcoef = function(glmobj)
{
  regtab  = summary(glmobj)$coef
  expcoef = exp(regtab[,1])
  lower   = expcoef*exp(-1.96*regtab[,2])
  upper   = expcoef*exp(1.96*regtab[,2])
  cbind(expcoef, lower, upper)
}
