# Read the data into a dataframe, give names to the variables, and inspect the data:
gun = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/v17/gun.dat", 
  col.names = c("method", "phys", "team", "rounds"))
gun
# Check that the data correspond to those given in the exercise.

# Compute correlations:
cor(gun)
# How are the correlations between the covariates ("method","phys","team")?
# Can you explain the reason for this?
# How are the correlations between the covariates and "rounds"?

# There are no correlation between the covariates as the experiment is planned,
# and have a balanced design -> orthogonality between the covariates. There is
# only a strong correlation between the response and the 'method' covariate.

# Define the covariates as factors (categorical covariates):
gun$method  = factor(gun$method)
gun$phys    = factor(gun$phys)
gun$team    = factor(gun$team)

# Fit a model with main effects and interactions and write the anova table:
gfit = lm(rounds ~ method*phys*team, data = gun) 
anova(gfit)
# What does the anova table tell you? Which interactions and main effects 
# are significant?

# Look at the estimates:
summary(gfit)
# Give an interpretation of the (most important) estimates.

# We could be interested in predicting the rounds given specified levels of 
# method, physique and team, and find confidence interval for estimated 
# expected values as well as prediction intervals for new observations 
# given the levels of these factors.

# Set up a data frame for values where you would like to make predictions, e.g.
testdata = data.frame(method = factor(c(1,2,1,2)), 
                      phys = factor(c(1,1,2,3)), 
                      team = factor(c(1,2,3,1)))

# Then find fitted/predicted values for your favorite model gfitfav 
# from R-exercise 12 by
predict(gfitfav , newdata = testdata)

# Then obtain confidence intervals for the expected values at this levels of 
# the factors by
predict(gfitfav , newdata = testdata, interval = "confidence")

# Next find the corresponding prediction intervals by
predict(gfitfav , newdata = testdata, interval = "prediction")

# Compare and discuss the difference between the confidence and prediction intervals.

