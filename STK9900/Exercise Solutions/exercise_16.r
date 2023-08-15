
# STK 4900
# Exercise 16

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# In this exercise we will use the results from the opinion polls from Norstat
# from February and March 2017 to investigate the change in the support for some political parties.
# In February 2017 Norstat asked n1=935 individuals which party they would support
# if there had been election to the parliament tomorrow.
# Of these y1=230 would have voted Høyre.
# One month later, in March,  n2=927 persons were interviewed and y2=208 of these would have voted Høyre.
# See http://www.pollofpolls.no/?cmd=Maling&gallupid=3096

# a)
# We start out by estimating the change in the support for Høyre with a 95 % confidence interval (cf. slide 6)
# Try to program such an interval yourself in R.
# Use Slide 6 of Lecture 6.
n1 = 935
y1 = 230
p1 = y1/n1
se1 = sqrt(p1*(1-p1)/n1)
n2 = 927
y2 = 208
p2 = y2/n2
se2 = sqrt(p2*(1-p2)/n2)
se = sqrt(se1^2+se2^2)
change = p1 - p2
margin = qnorm(p=0.975)*se
lower = change - margin
upper = change + margin

cbind(change,margin,lower,upper)


# b)
# We then test the null hypothesis that the support for Høyre has not changed from February to March (cf. slide 8)
# Use Slide 8 of Lecture 6.
p = (y1+y2)/(n1+n2)
se0 = sqrt(p*(1-p)/n1+p*(1-p)/n2)
z = (p1-p2)/se0
p.val = 2*(1-pnorm(abs(z)))

cbind(z,p.val)

# Perform these commands and comment on the results.
# Is the null hypothesis rejected or not? How does this relate to the confidence interval computed earlier?


# c)
# R has a command for comparing two proportions
n1 = 935
n2 = 927
y1 = 230
y2 = 208
hoyre = matrix(c(y1,y2,n1-y1,n2-y2), nrow=2) # give the data for Høyre in a 2x2 table (cf. slide 10)
prop.test(hoyre, correct=F)

# Perform these commands and check that the results agree with those obtained earlier.
# The prop.test-command give a chi squared statistic, not a z-value as we computed earlier. What is the relation between the two?


# d)
# We will then take a look at the results for Senterpartiet (Sp). In February 80 of the 935  persons who were interviewed would have
# voted Senterpartiet; while in March 101 of the 927 interviewed would have voted Senterpartiet.

# Estimating the change in the support for Senterpartiet with a 95 % confidence interval.
# Also test the null hypothesis that the support for Senterpartiet has not changed from February to March.

# What are your conclusions?


# Method 1:
# Manual computation

# Confidence interval
n1 = 935
y1 = 80
p1 = y1/n1
se1 = sqrt(p1*(1-p1)/n1)
n2 = 927
y2 = 101
p2 = y2/n2
se2 = sqrt(p2*(1-p2)/n2)
se = sqrt(se1^2+se2^2)
change = p1 - p2
margin = qnorm(p=0.975)*se
lower = change - margin
upper = change + margin

cbind(change,margin,lower,upper)

# Hypothesis testing
p = (y1+y2)/(n1+n2)
se0 = sqrt(p*(1-p)/n1+p*(1-p)/n2)
z = (p1-p2)/se0
p.val = 2*(1-pnorm(abs(z)))

cbind(z,p.val)


# Method 2:
# Built-in solution in R
n1 = 935
n2 = 927
y1 = 80
y2 = 101
hoyre = matrix(c(y1,y2,n1-y1,n2-y2), nrow=2) # give the data for Høyre in a 2x2 table (cf. slide 10)
prop.test(hoyre, correct=F)
