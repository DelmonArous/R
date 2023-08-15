# Exercise 1:  one-sample t-test and confidence interval


# You may copy the commands below from the web-browser into the command-window of R (or into a R-script)
# A line that starts with # is a comment, and R will disregard such lines.
#At the lectures we looked an example with the age of mineral samples (cf. slide 3 from the lectures)
#We will in this exercise see how the computations for this examples may be done in R.
#Start by reading the data into R. This may be done by the command:

age=c(249, 254, 243, 268, 253, 269, 287, 241, 273, 306, 303, 280, 260, 256, 278, 344, 304, 283, 310)



#Compute mean, median and standard deviation:

mean(age)

median(age)

sd(age)

# Check that you get the same result as in the lectures (cf slide 3)



# Make a histogram (cf. slide 4)

hist(age)

# the argument "breaks" helps with the precision of the histogram

hist(age, breaks=10)



# Plot the empirical cumulative distribution function (cf. slide 4)

plot(ecdf(age))                                        # Basic plot

plot(ecdf(age),verticals=T, do.points=F)          # Nicer looking plot




# Compute min, first quartile, median, third quartile, and max (cf. slide 5)

quantile(age)

# You can also specify which quantiles

quantile(age, seq(0, 1, .1))




# Make a boxplot (cf. slide 5)

boxplot(age)



#We will then consider confidence intervals and hypothesis testing.
# We will illustrate direct calculations of the quantities involved as well as the use of a special R-command.



# Compute the 97.5% percentile of the t-distribution with 18 degrees of freedom:

qt(0.975,18)



# Compute lower and upper limit of the 95% confidence interval:

mean(age) - qt(0.975,18)*(sd(age)/sqrt(19))      # lower limit
mean(age) + qt(0.975,18)*(sd(age)/sqrt(19))     # upper limit



# Check that you get the same result as in the lectures (cf slide 18)



# Compute t-statistic:

tstat=(mean(age)-265)/(sd(age)/sqrt(19))       #t-statistic
tstat



# Compute P-value:

1-pt(tstat,18)


# Check that you get the same result as in the lectures (cf slide 22)
# p-value is the probability to observe a test statistic at least as extreme as the one we have, if H0 is true.



# R has readymade commands for t-tests with corresponding confidence intervals.
# Use the command "t.test" to compute the confidence interval (this gives a two-sided test):

t.test(age,mu=265)



# Use the command "t.test" to compute a one-sided test (this gives a one-sided confidence interval):

t.test(age,alternative="greater",mu=265)



# Check that you get the same results as above.




#### Exercise 1b ####

# Exercise 1b: Bootstrapping
# In this exercise we will explain bootstrapping in some more detail than the in Lecture 1.
# We will use the same data set as in Exercise 1, age of mineral samples.
# First generate one bootstrap sample by the command

# It is best to set a seed when producing pseudo-random numbers
set.seed(123)
bootsamp=sample(age,replace=T)

# Compare this sampled data with the original data by writing

sort(bootsamp)

sort(age)

# Compute the mean and median of the bootstrap sample and compare with the corresponding
# values in the original data.

mean(bootsamp)
mean(age)

median(bootsamp)
median(age)

# Draw another bootstrap sample and repeat the comparison.
# One will typically draw a large number of bootstrap samples, say 1000, calculate the
# statistic for which we want to find a confidence interval and use the 2.5 and 97.5 percentiles
# as the confidence limits. This can be done, for the mean, using the commands

bootagemean<-numeric(0)
for (i in 1:1000) bootagemean[i]<-mean(sample(age,replace=T))
sort(bootagemean)[c(25,975)]

# Compare with the results using the t-interval. Also, inspect the ten first bootstrap estimates
# of the mean by writing

bootagemean[1:10]

c(mean(age) - qt(0.975,18)*(sd(age)/sqrt(19)), mean(age) + qt(0.975,18)*(sd(age)/sqrt(19)))
sort(bootagemean)[c(25,975)]

# Increase the number of bootstrap samples to 10000 and 50000. Comment.
# Make bootstrap confidence intervals for the median.

bootagemedian <- sapply(1:1000, function(i) median(sample(age, replace=T)))
sort(bootagemedian)[c(25,975)]

