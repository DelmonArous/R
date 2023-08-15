# At the lectures we looked an example concerning the opinion poll from 
# February 2017 (cf. slide 2 from the lectures). We will in this exercise 
# consider this example further.

rm(list = ls(all = TRUE))

# Question a)
# Of the n=935 persons who were interviewed by Norstat, y=sum(yi)=309 would 
# have voted Ap. The following calculations reproduce the result from the 
# lectures (cf. slide 4)
n   = 935
y   = 309

# Sample proportion:
p   = y/n

# Standard error of the sample proportion:
se  = sqrt(p*(1-p)/n)

# 95% confidence interval: 
margin = 1.96*se
lower  = p - margin
upper  = p + margin

# Display results and check that you get the result from the lectures:
cbind(p, margin, lower, upper)

# Question b)

# In the opinion poll, 122 of the persons interviewed would have voted 
# Fremskrittspartiet (FrP) and 80 would have voted Senterpartiet (Sp).

# Repeat the calculations above for Fremskrittspartiet and Senterpartiet.
# How is the "margin of error" for these parties compared to the 
# "margin of error" for Ap (cf. slide 4)?

# Fremskrittspartiet (FrP)
y   = 122
p   = y/n
se  = sqrt(p*(1-p)/n)

margin = 1.96*se
lower  = p - margin
upper  = p + margin

cbind(p, margin, lower, upper)

# Senterpartiet (Sp)
y   = 80
p   = y/n
se  = sqrt(p*(1-p)/n)

margin = 1.96*se
lower  = p - margin
upper  = p + margin

cbind(p, margin, lower, upper)
