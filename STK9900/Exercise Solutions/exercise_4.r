
# STK 4900
# Week 1, Exercise 4

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# a) Read the data from the web
solvents.data = as.data.frame(read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/solvents.txt", header=T))

# let's peek at the data
str(solvents.data)

# 'type' should be a categorical variable
# Check
is.factor(solvents.data[,"type"])
# Convert into a categorical variable
solvents.data[,"type"] = factor(solvents.data[,"type"])
# Check again
is.factor(solvents.data[,"type"])

# In fact...
summary(solvents.data)


# Look at the first 6 rows of data
head(solvents.data)
# Last 6 rows of data
tail(solvents.data)
# Number of observations per type
table(solvents.data[,"type"])


# b) Make a boxplot of the data.
boxplot(rate~type, 
        data=solvents.data, 
        names=c("1. Aromatics", "2. Chloroalkanes", "3. Esters"))


# c) Perform a one-way ANOVA.
aov.obj = aov(rate~type, 
              data=solvents.data)

aov.obj

summary(aov.obj)
