rm(list = ls(all = TRUE))

# Read the data into a dataframe and give names to the variables:
cafe = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/exer3_1.dat")
names(cafe) = c("no", "sale")

# Take a look at the data (make sure they correspond to those given in 
# the exercise):
cafe

# Attach the dataframe (making the variables available):
attach(cafe)

# Make a plot of sale as a function of the number or dispensers:
plot(no,sale)
# Inspect the plot. How is the relation between the number of dispensers 
# and the coffee sale?

# Fit a straight line and draw it on the plot:
linfit = lm(sale ~ no)
summary(linfit) 
abline(linfit)
# How well does the straight line describe the relation between the number 
# of dispensers and the coffee sale?

# Fit a second order polynomial:
sqfit = lm(sale ~ no + I(no^2))
summary(sqfit)

# Compute and draw the fitted second order polynomial:
x = seq(0,7,0.1)
koef = lm(sale ~ no + I(no^2))$coef
lines(x, koef[1] + koef[2]*x + koef[3]*x^2, lty = 2)
# Do the straight line or the second order polynomial provide the best 
# description of the relation between the number of dispensers and the 
# coffee sale?

# Compare adjusted R-squared
summary(linfit)$adj.r.squared
summary(sqfit)$adj.r.squared

# All coefficients are statistically significant. The intercepts aren't 
# physically impossible since they are not negative. And the estimated 
# coefficient for 'no' is higher for the quadratic model as the quadratic 
# coefficient is negative. This implies that there is a correlation between
# the 'no' and 'no^2' terms. If we consider the same constant value for the 
# 'no' covariate, then there is a negative (significant) effect on the sales 
# from the 'no^2' covariate that adjusts for the linear term, more strongly for
# larger number of dispencers
