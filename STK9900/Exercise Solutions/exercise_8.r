
# STK 4900
# Week 1, Exercise 8

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# Read the data into a dataframe and give names to the variables.
coffee.data = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/exer3_1.dat", header=FALSE)
names(coffee.data)=c("n.dispenser","sale")

# Take a look at the data.
coffee.data


# Create a scatter plot.
plot(x=coffee.data[,"n.dispenser"],
     y=coffee.data[,"sale"],
     xlab="Number of dispensers",
     ylab="Coffee sales")

# Inspect the plot. How is the relation between the number of dispensers and the coffee sale?

# Fit a linear regression.
lm.obj = lm(sale~n.dispenser, data=coffee.data)
summary(lm.obj)

# Fit a second order polynomial.
poly.lm.obj = lm(sale~n.dispenser+I(n.dispenser^2), data=coffee.data)
summary(poly.lm.obj)


# Visualize two models.

# Create a scatter plot.
plot(x=coffee.data[,"n.dispenser"],
     y=coffee.data[,"sale"],
     xlab="Number of dispensers",
     ylab="Coffee sales",
     main="")
# Plot the fitted linear model.
abline(lm.obj, col="blue")
# Notice that abline gives a warning now.
abline(poly.lm.obj, col="red")

# Plot the fitted polynomial model.
# Create a grid for x-value
x.grid = seq(from=min(coffee.data[,"n.dispenser"])-1, to=max(coffee.data[,"n.dispenser"])+1, length.out=1000)
# Predict from the fitted model.
y.hat = predict(object=poly.lm.obj, newdata=data.frame(n.dispenser=x.grid))
# Help: how to predict from a fitted linear model.
?predict.lm
# Saving the coefficients of the polynomial model
koef = lm(sale~n.dispenser+I(n.dispenser^2), data=coffee.data)$coef


# Create a scatter plot.
plot(x=coffee.data[,"n.dispenser"],
     y=coffee.data[,"sale"],
     xlab="Number of dispensers",
     ylab="Coffee sales",
     main="")
# Plot the fitted linear model.
abline(lm.obj, col="blue")
# Plot the polynomila fit
lines(x.grid,
      koef[1]+koef[2]*x.grid+koef[3]*x.grid^2,
      col="green3")


# list of colors in R:
# https://www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/

# Display legend.
legend("topleft", legend=c("linear model", "2. poly"), col=c("blue", "green3"), lty=c(1,1), cex=0.8)


# How well does the straight line describe the relation between the number of dispensers and the coffee sale?
# Does the linear model or the second order polynomial provide the best description of the relation between the number of dispensers and the coffee sale?
