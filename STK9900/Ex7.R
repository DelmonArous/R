rm(list = ls(all = TRUE))

# Load the "MASS" library into R by the command:
library(MASS)

# Generate 25 observations (x,y) from the bivariate normal distribution with 
# correlation 0.30 (see slide 24 from the lectures for explanation of the 
# bivariat normal. Actual commands for generating bivariat normal data are 
# given below) 

# Compute the Pearson correlation coefficient and plot the observations:
n_vec   = c(25, 100, 400)
rho_vec = c(0.30, 0.60, 0.90)

for (rho in rho_vec) {
  for (n in n_vec) {
    
    m   = matrix(c(0,0), nrow = 2)
    S   = matrix(c(1,rho,rho,1), nrow = 2)
    obs = mvrnorm(n, m, S)
    x   = obs[ , 1]
    y   = obs[ , 2]
    cor(x, y)
    plot(x, y,  main = paste("rho=", rho, ", n=", n))

  }
}


