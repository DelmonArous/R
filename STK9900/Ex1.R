# Start by reading the data
age = c(249, 254, 243, 268, 253, 269, 287, 241, 273, 306, 303, 280, 260, 256, 278, 344, 304, 283, 310)
n   = length(age)

# Compute mean, median and standard deviation: 
mu  = mean(age)
med = median(age)
SD  = sd(age)

# Make a histogram (cf. slide 4)
hist(age)

# Plot the empirical distribution function (cf. slide 4)
plot(ecdf(age))                                 # Basic plot
plot(ecdf(age), verticals = T, do.points = F)   # Nicer looking plot

# Compute min, first quartile, median, third quartile, and max (cf. slide 5)
quantile(age)

# Make a boxplot (cf. slide 5)
boxplot(age)

# Compute the 97.5% percentile of the t-distribution with 18 degrees of freedom:
qt(0.975, n-1)

# Compute lower and upper limit of the 95% confidence interval:
# (Check that you get the same result as in the lectures (cf slide 18))
L = mean(age) - qt(0.975, n-1)*(sd(age)/sqrt(n))     # lower limit
U = mean(age) + qt(0.975, n-1)*(sd(age)/sqrt(n))     # upper limit

# Compute t-statistic:
tstat = (mu-265)/(SD/sqrt(n))       # t-statistic
tstat

# Compute P-value:
# (Check that you get the same result as in the lectures (cf slide 22))
1-pt(tstat, n-1)

# R has readymade commands for t-tests with corresponding confidence intervals.
# Use the command "t.test" to compute the confidence interval (this gives a two-sided test):
t.test(age, mu = 265)

# Use the command "t.test" to compute a one-sided test (this gives a one-sided confidence interval):
# (Check that you get the same results as above)
t.test(age, alternative = "greater", mu = 265)






