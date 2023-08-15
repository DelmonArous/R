rm(list = ls(all = TRUE))

# Read the data into R by the command:
hers.sample = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/v17/hers.sample.txt", 
  header = T)

# Make a plot of the systolic blood pressure versus age:
plot(hers.sample$age, hers.sample$sbp, pch = 19)

# Fit a linear regression model to the data, using systolic blood pressure as 
# the outcome and age as the predictor:
hers.fit.b = lm(sbp ~ age, data = hers.sample)
summary(hers.fit.b)
abline(hers.fit.b)

# In order to get a more meaningful interpretation of the intercept, it may 
# be useful to "center" the predictor by subtraction of its mean (or a 
# quantity close to its mean).  Fit a linear regression using "age - 67" as 
# predictor:
hers.fit.c = lm(sbp ~ I(age-67), data = hers.sample)
summary(hers.fit.c)
abline(hers.fit.c)

# Sometimes it may be useful to use another unit for measuring a predictor 
# than the one given in the data file. To illustrate this point, we will here 
# fit a linear regression model where age is measured in units of then years:
hers.fit.d = lm(sbp ~ I(age/10), data = hers.sample)
summary(hers.fit.d)
abline(hers.fit.d)