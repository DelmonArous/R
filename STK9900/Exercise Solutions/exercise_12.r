
# STK 4900
# Exercise 12

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# Read the data.
gun.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/v17/gun.dat", col.names=c("method","phys","team","rounds"))
# Take a look at the data.
gun.data


# a)
# Compute correlations.
cor(gun.data)
# How are the correlations between the covariates ("method","phys","team")?
# Can you explain the reason for this?
# How are the correlations between the covariates and "rounds"?

# Hint: create a scatter plot matrix.
# Open a pdf device.

plot(gun.data)

table(gun.data$method)
table(gun.data$phys)
table(gun.data$team)



# b), c), d)
# Define the covariates as factors (categorical covariates).
gun.data[,"method"] = factor(gun.data[,"method"])
gun.data[,"phys"] = factor(gun.data[,"phys"])
gun.data[,"team"] = factor(gun.data[,"team"])

# Now the correlation command gives an error - with categorical variables is not possible
# to compute the "usual" correlation, we must use Kendall Tau for example.
cor(gun.data)

# Fit a model with main effects and interactions and write the anova table:
lm.obj = lm(rounds~method*phys*team, data=gun.data)
summary(lm.obj)

# Notice the aren't any "1"s in the model, those are the baselines!


anova(lm.obj)


# What does the anova table tell you? (See slide 33, Lecture 5)
# Which interactions and main effects are significant?




#### Part b ####

# For the data in R-exercise 12 on the gun data we could be interested in predicting the rounds given specified levels of method, physique and team.
# and find confidence interval for estimated expected values as well as prediction intervals for new observations given the levels of these factors.


# a) Set up a data frame for values where you would like to make predictions, e.g.
gun.test=data.frame(method=factor(c(1,2,1,2)),
                    phys=factor(c(1,1,2,3)),
                    team=factor(c(1,2,3,1)))

# Then find fitted/predicted values for your favourite model gfitfav from R-exercise 12 by
predict(lm.obj, newdata=gun.test)


# b) Then obtain confidence intervals for the expected values at this levels of the factors by
predict(lm.obj, newdata=gun.test, interval="confidence")


# c) Next find the corresponding prediction intervals by
predict(lm.obj, newdata=gun.test, interval="prediction")


# Compare and discuss the difference between the confidence and prediction intervals.

