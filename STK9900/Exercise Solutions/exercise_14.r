
# STK 4900
# Exercise 14

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# a)
# Read data.
rock.data = read.table(file="http://www.uio.no/studier/emner/matnat/math/STK4900/data/rock.dat", col.names=c("area","peri","shape","perm"), header=FALSE)
# Take a look at the data. 
rock.data

# Descriptive statistics, histogram and boxplot of "perm" (make one plot at a time):
summary(rock.data)

par(mfrow=c(1,2))
hist(rock.data[,"perm"])
boxplot(rock.data[,"perm"])
par(mfrow=c(1,1))

# How is the distribution of "perm"? Symetric or skewed?


# Try a log-transform:
rock.data[,"perm.log"] = log(rock.data[,"perm"])
summary(rock.data[,"perm.log"])


par(mfrow=c(1,2))
hist(rock.data[,"perm.log"])
boxplot(rock.data[,"perm.log"])
par(mfrow=c(1,1))


# How is the distribution of "perm.log"? Symetric or skewed?
# (We choose to use "perm.log" in question c).)



# b)
# Compute the correlation between the three covariates and make scatterplots of the covariates.
cor(rock.data[,1:3])


# Create a scatter plot
# Scatter plot matrix
plot(rock.data[,1:3])


# Is there a strong correlation between some of the covariates?
# How could this affect the regression analysis?


# c)
# Create scatter plots between 'perm.log' and other variables.


# 1 by 3 grid of plots.
par(mfrow=c(1,3))
# Scatter plot: 'perm.log' vs 'area'
plot(x=rock.data[,"area"],
     y=rock.data[,"perm.log"],
     main = "")
# Scatter plot: 'perm.log' vs 'peri'
plot(x=rock.data[,"peri"],
     y=rock.data[,"perm.log"],
     main = "")
# Scatter plot: 'perm.log' vs 'shape'
plot(x=rock.data[,"shape"],
     y=rock.data[,"perm.log"],
     main = "")
par(mfrow=c(1,1))


# Alternative: correlation matrix
cor(rock.data[c(1,2,3,5),c(1,2,3,5)])

# Which of the covariates seem to be of importance?


# Linear regression of "perm.log" using all three covariates.
# (the options "x=TRUE" and "y=TRUE" are included because of the computations in question d)
lm.obj.1 = lm(perm.log~area+peri+shape, data=rock.data, x=TRUE, y=TRUE)
summary(lm.obj.1)
names(lm.obj.1)

# Which of the covariates seem to be of importance?
# Is this in agreement with what you could see from the scatterplots?

lm.area <- lm(perm.log~area, data=rock.data)
summary(lm.area)

lm.peri <- lm(perm.log~peri, data=rock.data)
summary(lm.peri)

lm.shape <- lm(perm.log~poly(shape, 2), data=rock.data)
summary(lm.shape)



# d)
# According to their t-values the order of the covariates is: peri > area > shape
# Fit simpler regression models (according to the given order).
lm.obj.2 = lm(perm.log~area+peri, data=rock.data, x=TRUE, y=TRUE)
lm.obj.3 = lm(perm.log~peri, data=rock.data, x=TRUE, y=TRUE)


# A function that computes cross-validated R-squared.
cv.R2.func = function(lm.obj) {
  
  # Safety check
  if (("x" %in% names(lm.obj)) == FALSE) {stop("The given 'lm.obj' doesn't contain attribute 'x'.")}
  if (("y" %in% names(lm.obj)) == FALSE) {stop("The given 'lm.obj' doesn't contain attribute 'y'.")}
    
  # Extract 'x' and 'y' from the given 'lm.obj'.
  x = lm.obj$x
  y = lm.obj$y
  
  # Compute cross-validated R-squared.
  a = t(x)%*%x
  d = diag(1/eigen(a)$values)
  e = eigen(a)$vector
  inv.a = e%*%d%*%t(e)
  v = diag(x%*%inv.a%*%t(x))
  SSkryss = sum((lm.obj$res/(1-v))^2)
  SStot = sum((y-mean(y))^2)
  cross.validated.R2 = 1 - SSkryss/SStot
  
  # Return the cross validated R-squared.
  return(cross.validated.R2)
}

# Compute cross-validated R2:
cv.R2.func(lm.obj.1)
cv.R2.func(lm.obj.2)
cv.R2.func(lm.obj.3)
# Which of the three models give the best prediction?
# Compare the cross-validated R2 with R2 from question c). What do you see?


# e)
# Fit a model with second order terms and interactions:
lm.obj.4 = lm(perm.log~area+peri+shape+I(area^2)+I(peri^2)+I(shape^2)+area:peri+area:shape+peri:shape, data=rock.data, x=T, y=T)
summary(lm.obj.4)

# Which variables seem to be most important?
# According to their t-values the variables are ordered as: area:peri > area^2 > peri > area > peri^2 > area:shape > peri:shape > shape^2 > shape

# Fit models according to this ordering:
lm.obj.5 = lm(perm.log~area+peri+I(area^2)+I(peri^2)+I(shape^2)+area:peri+area:shape+peri:shape, data=rock.data, x=T, y=T)
lm.obj.6 = lm(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri+area:shape+peri:shape, data=rock.data, x=T, y=T)
lm.obj.7 = lm(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri+area:shape, data=rock.data, x=T, y=T)
lm.obj.8 = lm(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri, data=rock.data, x=T, y=T)
lm.obj.9 = lm(perm.log~area+peri+I(area^2)+area:peri, data=rock.data, x=T, y=T)
lm.obj.10 = lm(perm.log~peri+I(area^2)+area:peri, data=rock.data, x=T, y=T)
lm.obj.11 = lm(perm.log~I(area^2)+area:peri, data=rock.data, x=T, y=T)
lm.obj.12 = lm(perm.log~area:peri, data=rock.data, x=T, y=T)

# Compute cross-validated R2 for the models:
cv.R2.func(lm.obj.5)
cv.R2.func(lm.obj.6)
cv.R2.func(lm.obj.7)
cv.R2.func(lm.obj.8)
cv.R2.func(lm.obj.9)
cv.R2.func(lm.obj.10)
cv.R2.func(lm.obj.11)
cv.R2.func(lm.obj.12)

# Which model gives the best prediction?


# f)
# Use the same "recipe" as in question e).
# Which model gives the best prediction? Is it the same one as in question e?

# K-Fold CV function
kFoldCV <- function(lm.formula, K=5, seed=123) {
  n <- nrow(rock.data)
  
  kFolds <- rep(1:K, length.out=n)
  kFolds <- sample(kFolds, size=n, replace=F)
  MSEs <- c()
  
  for(k in 1:K) {
    Ids <- which(kFolds==k)
    
    train.data <- rock.data[-Ids, ]
    test.data <- rock.data[Ids, ]
    
    lm.cv <- lm(lm.formula, data=train.data)
    yhat <- predict(lm.cv, newdata=test.data)
    
    MSEs[k] <- mean((test.data$perm.log-yhat)^2)
  }
  
  return(mean(MSEs))
}


# Let's use it on the models of point e)
kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+I(shape^2)+area:peri+area:shape+peri:shape))
kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri+area:shape+peri:shape))
kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri+area:shape))
kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri))
kFoldCV(formula(perm.log~area+peri+I(area^2)+area:peri))
kFoldCV(formula(perm.log~peri+I(area^2)+area:peri))
kFoldCV(formula(perm.log~I(area^2)+area:peri))
kFoldCV(formula(perm.log~area:peri))


# Now we repeat CV 10 times and average over result
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+I(shape^2)+area:peri+area:shape+peri:shape), seed=i)))
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri+area:shape+peri:shape), seed=i)))
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri+area:shape), seed=i)))
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~area+peri+I(area^2)+I(peri^2)+area:peri), seed=i)))
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~area+peri+I(area^2)+area:peri), seed=i)))
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~peri+I(area^2)+area:peri), seed=i)))
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~I(area^2)+area:peri), seed=i)))
mean(sapply(1:10, function(i) kFoldCV(formula(perm.log~area:peri), seed=i)))



# g)
# Make various plots of the residuals.
# See previous exercise for help with R-commands.





