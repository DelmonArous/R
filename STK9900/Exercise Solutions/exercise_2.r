#Exercise 2:  two-sample t-test and confidence interval



# At the lectures we looked an example on bone mineral density (cf. slide 25 from the lectures)
# We will in this exercise see how the computations for this example may be done in R.

# Start by reading the data into R. This may be done by the command:

cont=c(0.228, 0.207, 0.234, 0.220, 0.217, 0.228, 0.209, 0.221, 0.204, 0.220, 0.203, 0.219, 0.218, 0.245, 0.210)
treat=c(0.250, 0.237, 0.217, 0.206, 0.247, 0.228, 0.245, 0.232, 0.267, 0.261, 0.221, 0.219, 0.232, 0.209, 0.255)



# Find the means and standard deviations, and check that you get the same results as in the lectures (slide 25)

summary(cont)
summary(treat)



# Use the command "t.test" to compute the confidence interval, t-statistic and P-value:

t.test(treat, cont, var.equal=TRUE)


# Make sure that you understand the output from the command and check that you get the same results as in the lectures (slides 27 and 28)



# Optional: Use the formulas given on slide 26 to compute the pooled standard deviation, the standard error of the effect of treatment, and the 95% confidence interval.

sp <- sqrt(.5*sd(cont)^2+.5*sd(treat)^2)
se <- sp*sqrt(1/15+1/15)
c(mean(treat)-mean(cont)-qt(.975, df=28)*se,mean(treat)-mean(cont)+qt(.975, df=28)*se)


# Optional: Use the formulas given on slide 28 to compute the t-statistic and the P-value.

tStat <- (mean(treat)-mean(cont))/se
tStat

