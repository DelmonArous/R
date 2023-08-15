
# STK 4900
# Exercise 20

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# Read data.
insects.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/data/insects.txt", header=TRUE)
# Take a look at the data. 
insects.data

# a)
# Compute the proportion dead at the various doses and plot the proportions versus logdose.
insects.data[,"proportion"] = insects.data[,"DEAD"]/insects.data[,"NUMBER"]

# Draw a scatter plot between "LOGDOSE" and "proportion"
# Open a pdf device.
pdf("./Exercise_20_a.pdf", width = 6, height = 6)
# Scatter plot
plot(x=insects.data[,"LOGDOSE"], y=insects.data[,"proportion"], ylim=c(0,1), pch=16)
# Turn off the device.
dev.off()

# Try to make the plot look nice by giving appropriate labels for the x-axis and the y-axis.

# b)
# Fit a logistic regression model with logdose as covariate and look at the result:
# (Note that since we have grouped data, not binary, the response has to be specified as

# cbind(y,n-y) where n=number of individuals in a group and y=number of "successes" in the group.)
glm.obj = glm(cbind(DEAD,NUMBER-DEAD)~LOGDOSE, data=insects.data, family=binomial)
summary(glm.obj)


# c)
# grid over x axis to visualize predicted values from the model.
x.grid = seq(from=0.3, to=1.2, by=0.01)
# Compute the probabilities obtained from the fitted model.
y.hat = predict(object=glm.obj, newdata=data.frame(LOGDOSE=x.grid), type="response")

# Draw a scatter plot between "LOGDOSE" and "proportion"
# Open a pdf device.
pdf("./Exercise_20_c.pdf", width = 6, height = 6)
# Scatter plot
plot(x=insects.data[,"LOGDOSE"], y=insects.data[,"proportion"], ylim=c(0,1), pch=16)
# Visualize the fitted model.
points(x=x.grid, y=y.hat, type="l", col="blue")
# Turn off the device.
dev.off()


# d)
# Use the formula from slide 24, Lecture 6.
# Extract coefficients from the model.
beta.0.hat = glm.obj$coefficients[1]
names(beta.0.hat) = NULL
beta.1.hat = glm.obj$coefficients[2]
names(beta.1.hat) = NULL
# Target probability
p = 0.5
# Convert p into LD50.
LD50 = (log(p/(1-p)) - beta.0.hat)/beta.1.hat
print(LD50)
# Check whether the solution is correct.
predict(object=glm.obj, newdata=data.frame(LOGDOSE=LD50), type="response")
