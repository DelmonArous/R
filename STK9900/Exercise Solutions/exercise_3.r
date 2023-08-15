# R-help to Exercise 3


# Read the data into R and look at the data:

speed=scan("http://www.uio.no/studier/emner/matnat/math/STK4900/v17/exer2.dat")
speed



# QUESTION a)

# Plot the data in various ways:

hist(speed) # histogram
hist(speed, breaks=20) # a little bit more precise...

plot(ecdf(speed), verticals=T, do.points=F) # empirical distribution function

boxplot(speed)             # box plot



# What do the different plots tell you?
# Are there indications for "outliers" in the data?





# QUESTION b)



# Compute the (empirical) mean and the median, which are the two most common measures of location
# (alternatively you may use the command summary(speed) to compute the mean and the median)

mean(speed)                 #mean
median(speed)              #median

summary(speed)

plot(ecdf(speed), verticals=T, do.points=F)            #empirical distribution function

# We can highlight important quantiles, like the median, this way!
abline(h=.5, lty=3, col="red")
abline(v=median(speed), lty=3, col="red")

# What do the two measures of location tell you?



# QUESTION c)

# Compute the (empirical) standard deviation and the interquartile range, which are two common measures of spread:

sd(speed)                     # standard deviation
IQR(speed)                  # interquartile range



# What do the two measures of spread tell you?





# QUESTION d)


# Compute t-based confidence interval using all data:

t.test(speed)



# Compute t-based confidence interval without the two "outliers":

t.test(speed[speed>0])



# What do the two intervals tell you? Which one do you find most reasonable?


