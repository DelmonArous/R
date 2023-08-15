rm(list = ls(all = TRUE))

# In this exercise we will use data from the Heart and Estrogen/Progestin 
# Study (HERS), a clinical trial of hormone therapy for prevention of 
# recurrent heart attacks and death among post-menopausal women with 
# existing coronary heart disease.

# The aim of this exercise is to study how the change in low-density 
# lipoprotein (LDL) cholesterol over the first year of the HERS study 
# depends on the baseline value of LDL (i.e. the value at entry to the study) 
# and use of hormone therapy (HT)

# You may read the HERS data into R by the command:
hers = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/v17/hers.txt", 
  sep = "\t", header = T, na.strings = ".")

# Before we start doing our analysis we define the change in LDL, denoted LDLch:
hers$LDLch = hers$LDL1 - hers$LDL

# We also defined the centered LDL at baseline (by subtracting the mean 
# value 145 mg/dL), denoted cLDL
hers$cLDL = hers$LDL - 145 

# Fit a linear model with the change in LDL as the response and hormone therapy 
# (HT) and baseline LDL (not centered) as covariates:
fit.a = lm(LDLch ~ HT + LDL, data = hers)
summary(fit.a)
# Give an interpretation of the estimates in the fitted model.

# We then fit a model with HT and centered LDL at baseline as covariates:
fit.b = lm(LDLch ~ HT + cLDL, data = hers)
summary(fit.b)
# Compare the estimates for this model with those in question a.
# What is the interpretation of the estimate for the intercept?

# All the estimates are the same, except for the intercept (both beta0 and 
# se(beta0)). The interpretation of beta0 is when all the covariates are equal
# to zero. This is done to have more interpretable results.

# We then fit a model with interaction:
fit.c = lm(LDLch ~ HT + cLDL + HT:cLDL, data = hers)
summary(fit.c)
# Is there a significant interaction?
# Given an interpretation of the estimates in the model (cf. slide 19 of 
# Lecture 4). In particular, what is the effect of baseline LDL for those 
# with no hormone therapy? And what is the effect of baseline LDL for those 
# on hormone therapy?

# There is no interaction between the HT covariate and the 
# HT*cLDL as the HT covariate estimate is exactly the same as before. There is
# no significant interaction (correlation) between the cLDL and HT*cLDL 
# covariates as the estimated value of the cLDL covariate is approximately 
# the same - orthogonality. There is a better effect on the baseline LDL if 
# HT is used compared to no HT, as this will reduce LDL values.



