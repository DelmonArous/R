rm(list = ls(all = TRUE))

# Read the data into R by the command:
solvents = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/solvents.txt", 
  header = T)

# Make a boxplot of the data:
boxplot(rate~type, data = solvents)

# Defines diet to be a categorical variable
solvents$type = factor(solvents$type)

# One-way ANOVA
aov.solvents = aov(rate~type, data = solvents)
summary(aov.solvents)
