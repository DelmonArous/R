# Clean up the memory before start
rm(list = ls(all = TRUE))

# Read the data into R by the command:
blood = read.table(
  "https://www.uio.no/studier/emner/matnat/math/STK4900/v21/obliger/blood.txt", 
  header = T)

# Age variable is considered as categorical (factors):
blood$age = factor(blood$age)

### Problem 2a) ###

# Make a summary and boxplot of the blood pressure according to the level 
# of age group:
summary(blood$Bloodpr[blood$age == 1])
summary(blood$Bloodpr[blood$age == 2])
summary(blood$Bloodpr[blood$age == 3])
boxplot(blood$Bloodpr ~ blood$age)

### Problem 2b) ###

# One-way ANOVA
bloodANOVA = aov(Bloodpr ~ age, data = blood)
summary(bloodANOVA)

### Problem 2c) ###

#  Regression model with age group as categorical predictor variable
lm.blood = lm(Bloodpr ~ age, data = blood)
summary(lm.blood)
