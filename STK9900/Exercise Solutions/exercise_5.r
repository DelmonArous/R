
# STK 4900
# Week 1, Exercise 5

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# a) Read the data into R and make a plot of the data.

# Enter the data.
pef.data = data.frame(pef = c(494,395,516,434,476,413,442,433),
                      minipef = c(512,430,520,428,500,364,380,445))

# Show the entered data.
str(pef.data)
summary(pef.data)

# Create a scatter plot
plot(x=pef.data[,"pef"],
     y=pef.data[,"minipef"],
     xlab="pef",
     ylab="minipef")

# Or...
plot(pef.data[,"minipef"]~pef.data[,"pef"],
     xlab="pef",
     ylab="minipef",
     col="red",
     pch=2)


# b) Compute the Pearson correlation coefficient. (See slide 29 from the lectures for R help.)

# Compute correlation (use the built-in function.)
cor.val = cor(pef.data[,"pef"], pef.data[,"minipef"])
cor.val
# Compute correlation (manually)
cor.two = cov(pef.data[,"pef"], pef.data[,"minipef"])/(sd(pef.data[,"pef"])*sd(pef.data[,"minipef"]))
cor.two


# c) Find 95% confidence interval for the (theoretical) correlation coefficient.
# (See slide 31 from the lectures for R help.)

# Test for correlation.
cor.test.obj = cor.test(pef.data[,"pef"], pef.data[,"minipef"])
cor.test.obj

# See what is in the object.
names(cor.test.obj)
# Extract confidence interval.
cor.test.conf.int = cor.test.obj$conf.int
cor.test.conf.int


# d) Fit a linear regression model with "minipef" as the outcome and "pef" as the predictor.
# Give an interpretation of the estimated slope, i.e. the least square estimate for "pef".
# (See slide 36 from the lectures for R help.)

# Fit a linear model.
lm.obj = lm(minipef~pef, data=pef.data)
# Print the result.
summary(lm.obj)

# Extract estimated parameters
beta.0.hat = lm.obj$coefficients[1]
names(beta.0.hat) = NULL
beta.1.hat = lm.obj$coefficients[2]
names(beta.1.hat) = NULL


# Create a scatter plot.
plot(x=pef.data[,"pef"],
     y=pef.data[,"minipef"],
     xlab="pef",
     ylab="minipef")
# Plot the fitted linear model.
abline(a=beta.0.hat, b=beta.1.hat, col="red")

# Or...
plot(x=pef.data[,"pef"],
     y=pef.data[,"minipef"],
     xlab="pef",
     ylab="minipef")

abline(lm.obj, col="red")


# e) Multiply the estimated slope from question d) with the ratio of the empirical standard deviations of "pef" and "minipef".
# Compare the result with the Pearson correlation from question b). What do you see?

# Compute the quantity.
beta.1.hat*(sd(pef.data[,"pef"])/sd(pef.data[,"minipef"]))
