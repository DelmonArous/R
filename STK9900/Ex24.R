# Clean up the memory before start
rm(list=ls(all=TRUE))

# Read the data into a dataframe:
melanom = read.table(
  "http://www.uio.no/studier/emner/matnat/math/STK4900/data/melanoma.txt", 
  header = T)

# Attach the R-library for survival analysis (if you have not already done so):
library(survival)


# QUESTION a)
# Commands:
fit.sex = survfit(Surv(lifetime, status == 1) ~ sex, data = melanom)
plot(fit.sex, lty = 1:2, mark.time = F)
survdiff(Surv(lifetime,status == 1) ~ sex, data = melanom)


# QUESTION b) 

# Take a look at the distribution of grthick (grouped tumor thickness)
table(melanom[,"grthick"])

# Compute Kaplan-Meier estimates for the two groups (males, females) (without confidence intervals).
surv.obj = survfit(Surv(lifetime, status == 1) ~ grthick, data = melanom)
# Plot the Kaplan-Meier estimates.
plot(surv.obj, lty = 1:3, mark.time = FALSE)
# Log-rank test for difference between the groups:
survdiff(Surv(lifetime,status == 1) ~ grthick, data = melanom)


# Question c)

# Take a look at the distribution of ulcer (ulceration)
table(melanom[,"ulcer"])

# Compute Kaplan-Meier estimates for the two groups (males, females) (without confidence intervals).
surv.obj = survfit(Surv(lifetime, status == 1) ~ ulcer, data = melanom)
# Plot the Kaplan-Meier estimates.
plot(surv.obj, lty = 1:2, mark.time = FALSE)
# Log-rank test for difference between the groups:
survdiff(Surv(lifetime, status == 1) ~ ulcer, data = melanom)


# QUESTIONS d) and e)
# Commands for sex
coxfit.sex = coxph(Surv(lifetime, status == 1) ~ factor(sex), data = melanom)
summary(coxfit.sex)

# The other covariates in question d are done by similar commands, and so is question e

