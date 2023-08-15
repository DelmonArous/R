
# STK 4900
# Week 1, Exercise 7

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# Check and change working directory.

# In this exercise we will perform some simulations that may help to get intuition on the concept
# of correlation and the performance of the Pearson correlation coefficient.
# In order to keep it simple, we will assume throughout the exercise that the expected values
# are 0 and the standard deviations are 1.

# Before you start doing the computations, you have to load the "MASS" library into R by the command.
library(MASS)


# a)
# Generate 25 observations (x,y) from the bivariate normal distribution with correlation 0.30
# (See slide 24 from the lectures for explanation of the bivariat normal.
# Actual commands for generating bivariat normal data are given below.)
# Compute the Pearson correlation coefficient and plot the observations.

# Number of data points to simulate.
n = 25
# Correlation value that we want to use for simulation.
rho = 0.30

# The true mean vector.
mu = matrix(data=c(0,0), nrow=2, ncol=1)
# The true covariance matrix.
S = matrix(data=c(1,rho,rho,1), nrow=2, ncol=2)
# The non-diagonal elements of "S" should be covariance.
# However, we used Pearson correlation coefficient (rho) in the non-diagonal elements of "S".
# Why is this ok in this case?

# Simulate data
# Use set.seed() for replicability.
set.seed(1)
sim.data = mvrnorm(n=n, mu=mu, Sigma=S)


# Display the simulated data
head(sim.data)
tail(sim.data)
summary(sim.data)

# Compute empirical correlation.
sample.cor = cor(sim.data[,1], sim.data[,2])
sample.cor

# Repeat the commands a number of times. Note how the Pearson correlation coefficient and the plot vary.

# Create a scatter plot
# Open a pdf device.
#pdf("./Exercise_7_a.pdf", width = 6, height = 6)
# Create a scatter plot.
plot(x=sim.data[,1],
     y=sim.data[,2],
     xlab="x",
     ylab="y",
     main = paste("Correlation = ", round(x=sample.cor, digits=2)))
# Turn off the device.
#dev.off()


# More extensive plot of the result.

# Create a 2 by 1 grid of plots.
par(mfrow=c(2,1))

# Histogram of x
hist(sim.data[,1], xlab="x", main="")
# Display sample mean.
abline(v=mean(sim.data[,1]), col="blue", lty=2)
# Display true mean.
abline(v=mean(mu[1]), col="blue", lty=1)
# Display legend.
legend("topright", 
       legend=c("sample mean", "true mean"), 
       col=c("blue", "blue"), 
       lty=c(2,1), 
       cex=0.8)

# Histogram of y
hist(sim.data[,2], xlab="y", main="")
# Display sample mean.
abline(v=mean(sim.data[,2]), col="blue", lty=2)
# Display true mean.
abline(v=mean(mu[2]), col="blue", lty=1)
# Display legend.
legend("topleft", legend=c("sample mean", "true mean"), col=c("blue", "blue"), lty=c(2,1), cex=0.8)

par(mfrow=c(1,1))


# Create a scatter plot.
plot(x=sim.data[,1],
     y=sim.data[,2],
     xlab="x",
     ylab="y",
     main=paste("Correlation = ", round(x=sample.cor, digits=2)))



# b) Repeat a) for correlation 0.60 and correlation 0.90.
# Note how the plots look like when the correlation is 0.60 and 0.90.

# Since we have to repeat the simulation multiple times, we rewrite the simulation as a function.


# This function simulates 'n' data points from bivariate normal
# distribution (with mean = 0  variance = 0 and rho, which is to be specified by the user).
# Then, it creates a plot.
# This function takes following argumets as its inputs:
#   n: Number of data points to simulate
#   rho: Correlation value that we want to use for simulation
#   seed: Seed for replicability.

simulate.and.plot.func = function(n, rho=0, seed) {
  
  # The true mean vector.
  mu.vec = matrix(data=c(0,0), nrow=2, ncol=1)
  # The true covariance matrix.
  cov.mat = matrix(data=c(1,rho,rho,1), nrow=2, ncol=2)
  
  # Simulate data
  # Use set.seed() for replicability.
  set.seed(seed)
  sim.data = mvrnorm(n=n, mu=mu.vec, Sigma=cov.mat)
  
  # Compute empirical correlation.
  sample.cor = cor(sim.data[,1], sim.data[,2])
  
  # Create a scatter plot.
  plot(x=sim.data[,1],
       y=sim.data[,2],
       xlab="x",
       ylab="y",
       main=paste("Sample cor = ", round(x=sample.cor, digits=2), "\n", "n = ", n, ", True cor = ", rho),
       xlim=c(-3, 3),
       ylim=c(-3, 3))
  
} # End of function.

# Now, we can simply use the function we wrote.
# Create a 1 by 3 grid of plots.
par(mfrow=c(1,3))
simulate.and.plot.func(n=25, rho=.3, seed=1)
# Repeat a) for rho = 0.6
simulate.and.plot.func(n=25, rho=.6, seed=1)
# Repeat a) for rho = 0.9
simulate.and.plot.func(n=25, rho=.9, seed=1)
par(mfrow=c(1,1))



# c)
# Repeat a) and b) for n=100 and n=500
# Note how the variation in the Pearson correlation coefficient depends on the sample size.

# Open a pdf device.
#pdf("./Exercise_7_c.pdf", width = 12, height = 8)
# Create a 1 by 2 grid of plots.
par(mfrow=c(2,3))

# Simulation with n = 100, rho = 0.3
simulate.and.plot.func(n=100, rho=0.3, seed=1)
# Simulation with n = 100, rho = 0.6
simulate.and.plot.func(n=100, rho=0.6, seed=1)
# Simulation with n = 100, rho = 0.9
simulate.and.plot.func(n=100, rho=0.9, seed=1)

# Simulation with n = 500, rho = 0.3
simulate.and.plot.func(n=500, rho=0.3, seed=1)
# Simulation with n = 500, rho = 0.6
simulate.and.plot.func(n=500, rho=0.6, seed=1)
# Simulation with n = 500, rho = 0.9
simulate.and.plot.func(n=500, rho=0.9, seed=1)

par(mfrow=c(1,1))


