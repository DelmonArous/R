
# STK 4900
# Exercise 25

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# QUESTION a)
# Read the data and variable names into a data frame, and take a look at the data
serum.data = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/serum.txt", header=TRUE)
serum.data

# Plot serum response against time for the different persons (identified by colors) and treatments
# (identified by solid or dotted line) using the matplot command for multiple lines in a plot.
hours.mat = matrix(c(1,2,3,6), nrow=4, ncol=10)
drug.A = matrix(serum.data$serum[serum.data$drug==1], nrow=4)
drug.B = matrix(serum.data$serum[serum.data$drug==2], nrow=4)
serum.mat = cbind(drug.A, drug.B)

matplot(hours.mat, serum.mat, type="l", lty=rep(x=c(1,2), each=5), col=rep(x=1:5, times=2), xlab="Hours", ylab="Serum", lwd=2)
# Think about what the plot tells you!


# QUESTION c)
AUC.matrix = as.data.frame(matrix(0,nrow=5,ncol=2))
colnames(AUC.matrix) = c("drug.A","drug.B")
for (i in 1:5) {# for-loop over subjects
  for (j in 1:2) {# for-loop over drugs
    
    # Show which iteration we are in. (This is to help understanding how a for-loop works.)
    cat("i = ", i, ", j = ", j, "\n", sep="")
    # Extract serum levels for a specific subject under specific drug. 
    serum.in.hours = serum.data[(serum.data[,"subject"]==i & serum.data[,"drug"]==j),"serum"]
    # Compute AUC.
    AUC.matrix[i,j] = serum.in.hours[1] + serum.in.hours[2] + 2*serum.in.hours[3] + 1.5*serum.in.hours[4]
    # Clean up before we go to the next iteration.
    remove(serum.in.hours)
  } 
}
print(AUC.matrix)
# Perform the computations and make sure that you get the AUCs

# Perform a paired t-test
t.test(AUC.matrix[,"drug.A"], AUC.matrix[,"drug.B"], paired=TRUE)
# What you may conclude from this hypothesis testing?


# QUESTION d)
# In order to fit a random effects model we will use the nlme-package.
# If this has not been installed, you will have to do so. At the end of the introduction to R, it is described how you may install new R-packages.

# In order to fit the modell, you then give the commands:
library(nlme)      # load the nlme-package
lme.obj.1 = lme(serum~factor(drug)+factor(time), random=~1|subject, data=serum.data, method="ML")
summary(lme.obj.1)


# QUESTION e)
# In order to test if drug has an effect, you may fit a model without drug and use
# the anova-command to compare the model with the one from question d.

# Test the effect of drug
lme.obj.2 = lme(serum~factor(time), random=~1|subject, data=serum.data, method="ML")
summary(lme.obj.2)
anova(lme.obj.1, lme.obj.2)

# Test the effect of time
lme.obj.3 = lme(serum~factor(drug), random=~1|subject, data=serum.data, method="ML")
summary(lme.obj.3)
anova(lme.obj.1, lme.obj.3)

# QUESTION f)
# Use Slide 11, 12, Lecture 10
var.subj = as.numeric(VarCorr(lme.obj.2)[1,2])^2
var.sigma = as.numeric(VarCorr(lme.obj.2)[2,2])^2

var.subj/(var.subj+var.sigma)
