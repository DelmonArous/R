# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# QUESTION a)
# Read the data and variable names into a data frame, and take a look at the data
aserum = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/serum.txt", 
  header=T)
aserum

# Check that the data are the same as given in the table in the exercise
# Make sure that you understand how the data are coded.


# QUESTION b)
# Plot serum response against time for the different persons (identified by 
# colors) and treatments (identified by solid or dotted line) using the matplot 
# command for multiple lines in a plot.

hours.mat   = matrix(c(1,2,3,6), nrow = 4, ncol = 10)
drug.A      = matrix(aserum$serum[aserum$drug == 1], nrow = 4)
drug.B      = matrix(aserum$serum[aserum$drug == 2], nrow = 4)
serum.mat   = cbind(drug.A, drug.B)
matplot(hours.mat, serum.mat, type = "l", lty = c(1,1,1,1,1,2,2,2,2,2), 
        col = c(1,2,3,4,5,1,2,3,4,5), xlab = "Hours", ylab = "Serum", lwd = 2)

# Think about what the plot tells you!


# QUESTION c)
# We first compute the AUC for each person and each drug:
auc = matrix(0, nrow = 5, ncol = 2)
for (i in 1:5) # for-loop over subjects
  for (j in 1:2) # for-loop over drugs
  {
    # Show which iteration we are in. 
    # (This is to help understanding how a for-loop works.)
    cat("i = ", i, ", j = ", j, "\n", sep = "") 
    # Extract serum levels for a specific subject under specific drug. 
    s         = aserum$serum[aserum$subject == i & aserum$drug == j]
    print(s)
    # Compute AUC
    auc[i,j]  = s[1] + s[2] + 2*s[3] + 1.5*s[4]
    # Clean up before we go to the next iteration.
    remove(s)
  }
print(auc)
# Perform the computations and make sure that you get the AUCs

# Perform a paired t-test
t.test(auc[,1], auc[,2], paired = T)
# What you may conclude from this hypothesis testing?


# QUESTION d)
# In order to fit a random effects model we will use the nlme-package.
# In order to fit the modell, you then give the commands:
library(nlme)        # load the nlme-package
fit.random = lme(serum ~ factor(drug) + factor(time), random = ~1|subject, 
                 data = aserum, method = "ML")
summary(fit.random)


# QUESTION e)
# In order to test if drug has an effect, you may fit a model without drug 
# and use the anova-command to compare the model with the one from question d)

# Test the effect of drug
lme.obj.2 = lme(serum ~ factor(time), random=~1|subject, data = aserum, 
                method = "ML")
summary(lme.obj.2)
anova(fit.random, lme.obj.2)

# Test the effect of time
lme.obj.3 = lme(serum ~ factor(drug), random=~1|subject, data = aserum, 
                method = "ML")
summary(lme.obj.3)
anova(fit.random, lme.obj.3)


# QUESTION f)
# Use Slide 11, 12, Lecture 10
var.subj  = as.numeric(VarCorr(lme.obj.2)[1,2])^2
var.sigma = as.numeric(VarCorr(lme.obj.2)[2,2])^2

var.subj/(var.subj+var.sigma)