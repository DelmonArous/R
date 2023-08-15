# Clean up the memory before start
rm(list = ls(all = TRUE))

# Read the olympic data into R by the command:
olympic = read.table( 
  "https://www.uio.no/studier/emner/matnat/math/STK4900/data/olympic.txt", 
  sep = "\t", header = T, na.strings = ".")


### Problem 2b) ###

# We consider the Poisson regression model 
#       E(Y) = exp{1*log(w) + b0 + b1*x1 + b2*x2 + b3*x3 }, where
#       Y  = number of medals won by the nation under the 2000 Olympic games (= Total2000),
#       x1 = number of medals won by the nation in the 1996 Olymic games (= Total1996),
#       x2 = logarithm of the nation's population size per 1000 (= Log.population),
#       x3 = the per capita Gross Domestic Product (GDP) of the nation (= GDP.per.cap)
#       W  = number of athletes representing the nation (= Log.athletes)
# Note that log(w) appears as a sort of "covariate" where we know that the 
# regression coefficient takes the value 1. This is called an OFFSET.


# We fit the Poisson regression model with x1, x2 and x3 as covariates and 
# look at the summary:
fit.model1 = glm(Total2000 ~ offset(Log.athletes) + Total1996 + Log.population 
                 + GDP.per.cap, data = olympic, family = poisson)
summary(fit.model1)

# Make an Analysis of Deviance Table
anova(fit.model1, test = "Chisq")

# We fit the Poisson regression model with x2 and x3 as covariates and 
# look at the summary:
fit.model2 = glm(Total2000 ~ offset(Log.athletes) + Log.population + 
                   GDP.per.cap, data = olympic, family = poisson)
summary(fit.model2)

# Make an Analysis of Deviance Table
anova(fit.model2, test = "Chisq")

# We fit the Poisson regression model with x2 as covariate and look at the 
# summary:
fit.model3 = glm(Total2000 ~ offset(Log.athletes) + Log.population, 
                 data = olympic, family = poisson)
summary(fit.model3)

# Compute the correlation between the four covariates and make scatterplots 
# of the covariates:
cor(olympic[,2:6])
plot(olympic[,2:6])
