  
# STK 4900
# Exercise 15

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# a)
# Of the n=935 persons who were interviewed by Norstat, y=sum(y_i)=309 would have voted Ap
# The following calculations reproduce the result from the lectures (cf. slide 4)

n = 935
y = 309
p = y/n
se = sqrt(p*(1-p)/n)

margin= 1.96*se
lower = p - margin
upper = p + margin

cbind(p, margin, lower, upper)
# Do the calculations and check that you get the result from the lectures.

# b)
# In the opinion poll, 122 of the persons interviewed would have voted Fremskrittspartiet (FrP) and 80 would have voted Senterpartiet (Sp).
# Repeat the calculations above for Fremskrittspartiet and Senterpartiet.
# How is the "margin of error" for these parties compared to the "margin of error" for Ap (cf. slide 4)?

n = 122
y = 80
p = y/n
se = sqrt(p*(1-p)/n)

margin= 1.96*se
lower = p - margin
upper = p + margin

cbind(p, margin, lower, upper)

