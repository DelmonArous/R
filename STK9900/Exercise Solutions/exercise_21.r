
# STK 4900
# Exercise 21

# Clean up the memory before we start.
rm(list=ls(all=TRUE))

# The following function may be used to calculate quantities mentioned in 
# Slide 8 - 14, Lecture 8
poisson.expected.table.func = function(observed.frequency) {
  
  # Number of total observations.
  n = sum(observed.frequency)
  # Frequency.
  y = 0:(length(observed.frequency)-1)
  
  # y.bar, as defined in Slide 8, Lecture 8.
  mean.y = sum(y*observed.frequency)/n
  # s^2, as defined in Slide 8, Lecture 8.
  var.y = ((n-1)^(-1))*sum(observed.frequency*(y-mean.y)^2)
  # Coefficient of dispersion (CD), as defined in Slide 8, Lecture 8.
  CD = var.y/mean.y
  # Degrees of freedom, as defined in Slide 9, Lecture 8.
  df = length(observed.frequency) - 2
  
  # Poisson denisty, as defined in Slide 3, Lecture 8.
  prob.y = dpois(x=0:df, lambda=mean.y)
  
  # Note that "1 - ppois(q=df, lambda=mean.y)" and "1 - sum(prob.y)" should give
  # exact the same result. Why?
  
  # Aggregate the frequencies that are equal to or higher than
  # the highest observed number of frequency, as mentioned in Slide 10, Slide 11, Lecture 8.
  prob.y = c(prob.y, 1 - ppois(q=df, lambda=mean.y))
  # Alternative:
  #prob.y = c(prob.y, 1 - sum(prob.y))
  
  # Compute the expected number of frequency (E_j) by using the definition on Slide 10, Lecture 8.
  # (Recognize Poisson density from Slide 3, Lecture 10.)
  expected.frequency = n*prob.y
  
  # Summarize the result into a matrix.
  result.matrix = cbind(
    # Frequency.
    y,
    # Observed number of frequency.
    observed.frequency,
    # Expected number of frequency.
    expected.frequency,
    # Difference between observed number of frequency and expected number of frequency.
    observed.frequency-expected.frequency,
    # Relative difference between observed number of frequency and expected number of frequency.
    # This is the term inside the summation part of Pearson chi-squared statistic, in Slide 10, Lecture 8.
    (observed.frequency-expected.frequency)^2/expected.frequency
  )
  result.matrix = as.data.frame(result.matrix)
  # Give proper names to the columns.
  colnames(result.matrix) = c("y","Ob.freq","Ex.freq","O-E","(O-E)^2/E")
  # Round off some decimals
  result.matrix = round(x=result.matrix, digits=2)
  
  # Compute Pearson chi-squared statistic, as defined in Slide 10, Lecture 8.
  Chi.sq = sum(result.matrix[,"(O-E)^2/E"])
  # Compute the p-value that belongs to the Pearson chi-squared statistic.
  p.val = 1 - pchisq(q=Chi.sq, df=df)
  # what does "pchisq()" do?
  # Hint: Use "?pchisq"
  # Why should we use "1 - pchisq(q=Chi.sq, df=df)" and not "pchisq(q=Chi.sq, df=df)"?
  
  # Print some results.
  cat("mean.y = ", mean.y, "\n", sep="")
  cat("var.y = ", var.y, "\n", sep="")
  cat("CD = ", CD, "\n", sep="")
  cat("Chi.sq = ", Chi.sq, "\n", sep="")
  cat("p-value = ", p.val, "\n", sep="")
  cat("\n")
  
  # Return the result matrix.
  return(result.matrix)
}


# Use the function:
yeast = c(75,103,121,54,30,13,4)
poisson.expected.table.func(observed.frequency=yeast)

# Try to understand what the function does and how the function produces the output given.
# Make sure you understand the output, and use the output to answer the questions in the exercise.
# Repeat the analysis for the two other data sets in the exercise.

moss = c(100,9,6,8,1,0,2)
poisson.expected.table.func(observed.frequency=moss)

azuki = c(61,50,1)
poisson.expected.table.func(observed.frequency=azuki)
