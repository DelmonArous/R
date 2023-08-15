
# STK 4900
# Exercise 17

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# In this exercise we will use data from the Western Collaborative Group Study (WCGS), a large
# epidemiological study designed to study risk factors for coronary heart disease (CHD).
# More than 3000 men were recruited to the study, and a number of (potential) risk factors were recorded at entry.
# The men were then followed for about ten years and it was recorded if they developed CHD or not.
# The WCGS data are used in many of the examples in Section 3.4 and Chapter 6 of the course text book byVittinghoff et al.

# Read data.
WCGS.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/data/wcgs.txt", header=TRUE, sep="\t", na.strings=".")

# In this exercise we study the effect of smoking and age on the risk for CHD. We will only look at the estimates (standard errors, confidence intervals, tests etc will be considered in Exercise 18).

# a)
# We first look at the effect of smoking (which is a binary covariate: 0=non-smoker, 1=smoker).
# We start out by making a 2x2 table for the development of CHD or not for smokers and nonsmokers:

WCGS.table = table(WCGS.data[,"smoke"], WCGS.data[,"chd69"])
# From the table you will see that there were 1502 smokers and that 159 of these developed CHD.
# Of the 1652 non-smokers, 98 developed CHD

# We then estimate the proportion who developed CHD among smokers and non-smokers.
# We also estimate the relative risk and the odds ratio (cf. slides 17 and 19):

p1 = WCGS.table[2,2]/sum(WCGS.table[2,])
p0 = WCGS.table[1,2]/sum(WCGS.table[1,])

# Relative risk (Slide 17, Lecture 6)
RR = p1/p0
# Odds ratio (Slide 19, Lecture 6)
OR = (p1/(1-p1))/(p0/(1-p0))

cbind(p1,p0,RR,OR)

# Perform the commands and discuss what the results tell you.


# b)
# We then use a logistic regression model with smoke as covariate:
glm.obj.1 = glm(chd69~smoke, data=WCGS.data, family=binomial)
summary(glm.obj.1)
# From a fitted logistic model we may estimate the odds ratio corresponding to a one unit's increase in a covariate (cf. slide 25, 26). Use this result to compute the odds ratio for CHD for a smoker compared to a non-smoker. Compare with the result from question a.


# c)
# We then use logistic regression to study the effect of age (at entry to the study) for the risk of developing CHD.
glm.obj.2 = glm(chd69~age, data=WCGS.data, family=binomial)
summary(glm.obj.2)

# From a fitted logistic model we may estimate the odds ratio corresponding to a given increase in a covariate (cf. slide 25, 26). Compute the odds ratio for CHD for a one-year and a ten-year increase in age.
# Check that you get the same results as on pages 162-163 in the text book byVittinghoff et al.
