# Clean up the memory before start
rm(list = ls(all = TRUE))

# Read the crabs data into R by the command:
crabs = read.table( 
  "https://www.uio.no/studier/emner/matnat/math/STK4900/data/crabs.txt", 
  header = T, na.strings = ".")

### Problem 1b) ###

# Fit a logistic regression model with width as covariate:
fit.width = glm(y ~ width, data = crabs, family = binomial)
summary(fit.width)

### Problem 1c) ###

# Fit a logistic regression model with weight as covariate:
fit.weight = glm(y ~ weight, data = crabs, family = binomial)
summary(fit.weight)

# Fit a logistic regression model with color as covariate:
fit.color = glm(y ~ factor(color), data = crabs, family = binomial)
summary(fit.color)

# Fit a logistic regression model with spine as covariate:
fit.spine = glm(y ~ factor(spine), data = crabs, family = binomial)
summary(fit.spine)

### Problem 1d) ###

# Fit a logistic regression model with width, weight, color and spine as 
# covariate:
fit.all = glm(y ~ width + weight + factor(color) + factor(spine), data = crabs, 
              family = binomial)
summary(fit.all)

# Fit a logistic regression model by iteratively appending width, weight, color 
# and spine as covariate:
fit.model1 = glm(y ~ width, data = crabs, family = binomial)
fit.model2 = glm(y ~ width + weight, data = crabs, family = binomial)
fit.model3 = glm(y ~ width + weight + factor(color), data = crabs, 
                 family = binomial)
fit.model4 = glm(y ~ width + weight + factor(color) + factor(spine), 
                 data = crabs, family = binomial)

# To test the null hypothesis that q of the beta parameters are equal to zero, 
# we pairwise compare the (residual) deviances for model 1, model 2, model 3 
# and model 4:
anova(fit.model1, fit.model2, fit.model3, fit.model4, test = "Chisq")

# Make a summary of model 2; the logistic regression model with width and weight
# as numerical covariates:
summary(fit.model2)

# Compute the correlation between the four covariates and make scatterplots 
# of the covariates:
cor(crabs[,2:5])
plot(crabs[,2:5])


### Problem 1e) ###

# Fit a logistic regression model by sequentially appending all the covariates 
# including the interactions between each covariate:
fit.model5 = glm(y ~ width + weight + factor(color) + factor(spine) 
                 + width:weight + width:factor(color) + width:factor(spine) 
                 + weight:factor(color) + weight:factor(spine) 
                 + factor(color):factor(spine), data = crabs, family = binomial)
summary(fit.model5)

# Make an Analysis of Deviance Table
anova(fit.model5, test = "Chisq")
