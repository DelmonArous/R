
# STK 4900
# Week 1, Exercise 10

# Clean up the memory before we start.
rm(list=ls(all=TRUE))


# The Heart and Estrogen/Progestin Study (HERS) is a clinical trial of hormone therapy for prevention of recurrent
# heart attacks and death among post-menopausal women with existing coronary heart disease.
# The HERS data are used in many of the examples in Chapters 3 and 4 of the course text book byVittinghoff et al.
# In this exercise we will study how different variables may influence the glucose level in the blood for the non-diabetic women in the cohort,
# in particular we are interested to see if exercise may help to reduce the glucose level (cf. Section 4.1 in Vittinghoff et al.).

# Read the data.
HERS.data = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/hers.txt", header=TRUE, sep="\t", na.strings=".")
# Keep only women without diabetes.
HERS.data.1 = HERS.data[(HERS.data[,"diabetes"] == 0), ]
rownames(HERS.data.1) = NULL

# We will start out by investigating (in questions a-c) how the glucose levels are for women who exercise at least
# three times a week (coded as exercise=1) and women who exercise less than three times a week (coded as exercise=0).


# a)
# Make a summary and boxplot of the glucose levels according to the level of exercise:
# Create summary
summary(HERS.data.1[(HERS.data[,"exercise"] == 0), "glucose"])
summary(HERS.data.1[(HERS.data[,"exercise"] == 1), "glucose"])
table(HERS.data[,"exercise"])

# Create a boxplot.
# Draw boxplot.
boxplot(glucose~exercise, data=HERS.data.1, names=c("No exercise","Exercise"))

# Discuss what the summaries and boxplot tell you.

# Another way to explore the data might be using histograms to inspect distributions.
par(mfrow=c(2,1))
hist(HERS.data.1[HERS.data.1$exercise==0, "glucose"], breaks=40, xlim=c(40, 140), main="No exercise", xlab="glucose")
hist(HERS.data.1[HERS.data.1$exercise==1, "glucose"], breaks=40, xlim=c(40, 140), main="Exercise", xlab="glucose")
par(mfrow=c(1,1))



# b)
# Test whether there is a difference in glucose level and make a confidence interval:
t.test(glucose~exercise, var.equal=TRUE, data=HERS.data.1)
# What may you conclude for the test and the confidence interval?



# c)
# Perform a simple linear regression with glucose level as outcome and exercise as predictor:
lm.obj = lm(glucose~exercise, data=HERS.data.1)
summary(lm.obj)
# Discuss how the result of the simple linear regression relates to those in question b).



# d)
# The women who exercise at least three times a week and the women who exercise less than
# three times a week may differ in many ways. For example they may be younger and have
# a lower BMI (body mass index). We will therefore perform a multiple linear regression analysis
# where we adjust for these to variables.

# Perform a multiple linear regression with glucose level as outcome and exercise, age, and BMI as predictors.
lm.obj.2 = lm(glucose~exercise+age+BMI, data=HERS.data.1)
summary(lm.obj.2)
# Discuss the result of this analysis in relation with the result in question c). Also discuss how age and BMI influence the glucose level.



