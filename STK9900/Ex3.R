#speed   = c(28, -44, 29, 30, 24, 28, 37, 32, 36, 27, 26, 28, 
#            29, 26, 27, 22, 23, 20, 25, 25, 36, 23, 31, 32,
#            24, 27, 33, 16, 24, 29, 36, 21, 28, 26, 27, 27, 
#            32, 25, 28, 24, 40, 21, 31, 32, 28, 26, 30, 27,
#            26, 24, 32, 29, 34, -2, 25, 19, 36, 29, 30, 22,
#            28, 33, 39, 25, 16, 23)


# Read the data into R and look at the data:
speed   = scan("http://www.uio.no/studier/emner/matnat/math/STK4900/v17/exer2.dat")
n       = length(speed)

### Ex. 3a ###

# Make a histogram
hist(speed)

# Plot the empirical distribution function
plot(ecdf(speed), verticals = T, do.points = F)

# Make a boxplot
boxplot(speed)

### Ex. 3b ###

# Compute mean and median 
mu  = mean(speed)
med = median(speed)
summary(speed)

### Ex. 3c ###

# Compute standard deviation and the interquartile range (i.e. the difference 
# between the third and the first quartile). These are both measures of spread
SD  = sd(speed)
iqr = IQR(speed)
quantile(speed, 0.75) - quantile(speed, 0.25) # lik svar som iqr
summary(speed)[5] - summary(speed)[2]         # ogsa lik svar som iqr

### Ex. 3d ###

# Compute t-based confidence interval using all data:
t.test(speed)

# Compute t-based confidence interval without the two "outliers":
t.test(speed[speed > 0])

# R has readymade commands for t-tests with corresponding confidence intervals.
# Use the command "t.test" to compute the confidence interval (this gives a one-sided confidence interval)::
t.test(speed, alternative = "greater", mu = 33.02)
t.test(speed[speed > 0], alternative = "greater", mu = 33.02)

# Compute lower and upper limit of the 95% confidence interval:
# (Check that you get the same result as in the lectures (cf slide 18))
L = mu - qt(0.975, n-1)*(SD/sqrt(n))     # lower limit
U = mu + qt(0.975, n-1)*(SD/sqrt(n))     # upper limit

# Compute t-statistic:
tstat = (mu - 33.02)/(SD/sqrt(n))       # t-statistic
tstat

# Compute P-value:
1 - pt(tstat, n-1)
