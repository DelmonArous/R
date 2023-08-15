
# STK 4900
# Week 1, Exercise 6

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# a)
# Read the data into R (cf. above) and inspect the data. Make a plot of the systolic blood pressure versus age.
HERS.data = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v17/hers.sample.txt", header=T)

# First observations in the data.
head(HERS.data)

# Last observations in the data.
tail(HERS.data)

# Create a scatter plot
plot(x=HERS.data[,"age"],
     y=HERS.data[,"sbp"],
     xlab="age",
     ylab="sbp")

# Discuss what the plot tells you.



# b) Fit a linear regression model to the data, using systolic blood pressure as the outcome and age as the predictor.
# Fit a linear model.
lm.obj.1 = lm(sbp~age, data=HERS.data)
# Show the result.
summary(lm.obj.1)

# Extract estimated parameters
beta.0.hat = lm.obj.1$coefficients[1]
names(beta.0.hat) = NULL
beta.1.hat = lm.obj.1$coefficients[2]
names(beta.1.hat) = NULL

print(beta.0.hat)
print(beta.1.hat)

# Create a scatter plot with fitted regression line
plot(x=HERS.data[,"age"],
     y=HERS.data[,"sbp"],
     xlab="age",
     ylab="sbp")
# Plot the fitted linear model.
abline(a=beta.0.hat, b=beta.1.hat, col="blue")


# Give an interpretation of the slope of the regression line (i.e. the effect of age on sbp).
# Does age have a significant effect on systolic blood pressure?
# Can you give a meaningful interpretation to the intercept?

# c) In order to get a more meaningful interpretation of the intercept, it may be useful
# to "center" the predictor by subtraction of its mean (or a quantity close to its mean).
# Fit a linear regression using "age - 67" as predictor:
# Compare the least squares estimates with the ones you got in question b).
# How can you now interpret the intercept?

# Fit a linear regression with "age - 67" as predictor.
lm.obj.2 = lm(sbp~I(age-67), data=HERS.data)
# Show the result.
summary(lm.obj.2)

# Summary of the result.
result.table = data.frame(lm.1 = lm.obj.1$coefficients,
                          lm.2 = lm.obj.2$coefficients)
rownames(result.table) = c("beta.0.hat","beta.1.hat")
result.table

# Pay attention when you produce a plot *after* having transformed data.
plot(sbp~I(age-67), data=HERS.data)
abline(lm.obj.2, col="red")
abline(v=0, lty=3)

# d) Sometimes it may be useful to use another unit for measuring a predictor than the one given in the data file.
# To illustrate this point, we will here fit a linear regression model where age is measured in units of then years.

# Fit a linear regression with age/10 as predictor.
lm.obj.3 = lm(sbp~I(age/10), data=HERS.data)
# Show the result.
summary(lm.obj.3)

# Summary of the result.
result.table = cbind(result.table, lm.3 = lm.obj.3$coefficients)
result.table

# In this case, the scatterplot will be
plot(sbp~I(age/10), data=HERS.data)
abline(lm.obj.3, col="darkorange")

# Compare the least squares estimates with the ones you got in question b).
# How can you now interpret the slope?
