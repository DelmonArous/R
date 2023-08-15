#Start by reading the data into R. This may be done by the command:
cont  = c(0.228, 0.207, 0.234, 0.220, 0.217, 0.228, 0.209, 0.221, 0.204, 0.220, 0.203, 0.219, 0.218, 0.245, 0.210)
treat = c(0.250, 0.237, 0.217, 0.206, 0.247, 0.228, 0.245, 0.232, 0.267, 0.261, 0.221, 0.219, 0.232, 0.209, 0.255)
n1    = length(cont)
n2    = length(treat)

# Find the means and standard deviations, and check that you get the same 
# results as in the lectures (slide 25)
mu1   = mean(cont)
mu2   = mean(treat)
SD1   = sd(cont)
SD2   = sd(treat)

# Use the command "t.test" to compute the confidence interval, t-statistic and P-value:
t.test(treat, cont, var.equal = T)

# Make sure that you understand the output from the command and check that 
# you get the same results as in the lectures (slides 27 and 28)





# Optional: Use the formulas given on slide 26 to compute the pooled standard 
# deviation, the standard error of the effect of treatment, and the 95% 
# confidence interval.


# Optional: Use the formulas given on slide 28 to compute the t-statistic 
# and the P-value.