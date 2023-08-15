rm(list = ls(all = TRUE))

# You may read the HERS data into R and extract the women without diabetes 
# by the commands: 
hers = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/hers.txt", 
  sep = "\t", header = T, na.strings = ".")
hers$exercise = factor(hers$exercise)
hers.no = hers[hers$diabetes == 0, ]

# We will start out by investigating (in questions a-c) how the glucose 
# levels are for women who exercise at least three times a week 
# (coded as exercise=1) and women who exercise less than three 
# times a week (coded as exercise=0).

# Make a summary and boxplot of the glucose levels according to the level 
# of exercise:
hers.no$glucose
summary(hers.no$glucose[hers.no$exercise == 0])
summary(hers.no$glucose[hers.no$exercise == 1])
boxplot(hers.no$glucose ~ hers.no$exercise)
# Discuss what the summaries and boxplot tell you

# Test if there is a difference in glucose level and make a confidence 
# interval:
t.test(glucose ~ exercise, var.equal = T, data = hers.no)
# What may you conclude for the test and the confidence interval?

# The p-value is significant, telling that there is a significant difference 
# between the means of the two exercise groups. The 95% CI does not contain 
# zero, so it is also an indicator that there is a significance in the
# difference between the means. This is because of the large sample size, 
# although the means are very similar

# Perform a simple linear regression with glucose level as outcome and 
# exercise as predictor:
fit.c = lm(glucose ~ exercise, data = hers.no)
summary(fit.c)
# Discuss how the result of the simple linear regression relates to 
# those in question b.

# Here we get a significance of the 'exercise' covariate precicely because of 
# the large sample size - even if the means of the two exercise groups are
# very similar, the large sample size will pronounce even small effects.


# The women who exercise at least three times a week and the women who 
# exercise less than three times a week may differ in many ways. 
# For example they may be younger and have a lower BMI (body mass index). 
# We will therefore perform a multiple linear regression analysis where we 
# adjust for these to variables.

# Perform a simple linear regression with glucose level as outcome and 
# exercise, age, and BMI as predictors:
fit.d = lm(glucose ~ exercise + age + BMI, data = hers.no)
summary(fit.d)
# Discuss the result of this analysis in relation with the result in 
# question c. Also discuss how age and BMI influence the glucose level.

