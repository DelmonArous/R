# Clean up the memory before start
rm(list=ls(all=TRUE))

# In this exercise we will use data from the Western Collaborative Group Study 
# (WCGS), a large epidemiological study designed to study risk factors for 
# coronary heart disease (CHD). More than 3000 men were recruited to the study, 
# and a number of (potential) risk factors were recorded at entry. The men were 
# then followed for about ten years and it was recorded if they developed CHD 
# or not. The WCGS data are used in many of the examples in Section 3.4 and 
# Chapter 6 of the course text book byVittinghoff et al.

# Read the WCGS data into R by the command:
wcgs = read.table( 
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/wcgs.txt", 
  sep = "\t", header = T, na.strings = ".")


# Question a)
# Our starting point in this exercise is the model for CHD risk using the 
# predictors age (per 10 years), cholesterol (per 50 mg/dL), systolic blood 
# pressure (per 50 mmHg), body mass index (per 10 kg/m2), smoking (yes, no), 
# and behavioral pattern (with the four groups 1=A1, 2=A2, 3=B3, 4=B4).

# You may fit this model by the commands on slide 19:
wcgs$behcat = factor(wcgs$behpat)
wcgs.beh    = glm(chd69 ~ age_10 + chol_50 + sbp_50 + bmi_10 + smoke + behcat, 
                  data = wcgs, family = binomial, subset= (chol<600))
summary(wcgs.beh)

# Perform the commands and check that you get the results reported on slide 19.
# Comment on the results. In particular discuss the effects of the behavioral 
# patterns


# Question b)
# In question a) we consider a model with four behavioral groups. One may wonder 
# if it is sufficient to consider only two behavioral groups (A-behavior and 
# B-behavior). To this end we fit a model with the binary covariate dibpat 
# (which is coded as 0 for B3 and B4, and as 1 for A1 and A2):
wcgs.beh2 = glm(chd69 ~ age_10 + chol_50 + sbp_50 + bmi_10 + smoke + dibpat, 
                data = wcgs, family = binomial, subset = (chol<600))

# To test the null hypothesis that the effects behavioral patterns A1 and A2 
# are the same and the effects behavioral patterns B1 and B2 are the same, 
# we compare the (residual) deviances for the model in a) and the model 
# considered here:
anova(wcgs.beh2, wcgs.beh, test = "Chisq")

# Perform the command and explain what the output tells you (cf slide 26). 
# What do you conclude from the test?


# Question c)
# We then fit a model without behavioral pattern and compare the three models 
# in an analysis of deviance table:
wcgs.resc = glm(chd69 ~ age_10 + chol_50 + bmi_10 + sbp_50 + smoke, 
                data = wcgs, family = binomial, subset=  (chol<600))

# Perform the command and explain what you may conclude from the analysis of 
# deviance table?
anova(wcgs.resc, wcgs.beh2, wcgs.beh, test = "Chisq")


# Question d)
# Which of the models in question a, b, and c would you prefer?
# Use the expcoef-command from Exercise 18 to obtain the estimated odds ratios 
# with 95% confidence intervals for your preferred model. Discuss what these 
# odds ratios tell you about the risk factors for CHD.
