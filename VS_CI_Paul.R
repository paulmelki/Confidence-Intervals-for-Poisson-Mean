# ----------------------- MATHEMTATICAL STATISTICS II ------------------------ #
# ------------ Project 2: Confidence Interval for a Poisson Mean ------------- #
# ---------------------------------------------------------------------------- #
# ------- Yanxin LI ---- Paul MELKI ---- Mira RAHAL ---- Zefeng ZHANG -------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Created on 08/03/2020 ------------------------------------------------------ #
# Iteration 3.0.1 ---------- Last updated on 16/03/2020 -------- by Paul MELKI #
# R Version 3.6.2 on Windows 10 ---------------------------------------------- #
# ---------------------------------------------------------------------------- #


# ##### ---------------------------------------------------------------------- #
# 1. 
# PRELIMINARIES  
# ---------------------------------------------------------------------------- #

# 1.1. Libraries' Imports
library(COUNT)          # ver. 1.3.4
library(fitdistrplus)   # ver. 1.0-1.4
library(ggplot2)        # ver. 3.2.1
library(Cairo)

# 1.2. Data Import
# Our study will be conducted on the 'mdvis' dataset, found in the 'COUNT'
# package. 

# Include the dataset
data("mdvis")

# Let's take a closer look at the dataset
help("mdvis")
# We can see that the data frame contains 2227 observations and 13 variables.
# It appears that all variables are discrete variables, with the exception of 
# the variable "loginc" being a continuous one.
# The main variable studied in this dataset is the first column, "numvisit", 
# specifying the number of patient visits to a doctor's office during a 3 month
# period. It is, of course, a discrete variable and thus we can fit a Poisson
# distribution to it. 



# ##### ---------------------------------------------------------------------- #
# 2. 
# FITTING POISSON DISTRIBUTION TO THE DATA  
# ---------------------------------------------------------------------------- #

# 2.1. Export the first column of the data to a new variable
numvisit <- as.numeric(mdvis[, 1])

# 2.2. Look at the sample average of our obtained vector
numvisit.mean <- mean(numvisit)

# 2.3. Let's take a look at the distibution of the sample
plotdist(numvisit, "pois", para = list(lambda = numvisit.mean))

# 2.4. We now try to fit the Poisson distribution to the data using 2 methods:
# Maximum Likelihood Estimation (MLE) and Method of Moments Estimation (MME).
numvisit.MLE.fit <- fitdist(data = numvisit,
                            distr = "pois",
                            method = "mle")
numvisit.MME.fit <- fitdist(data = numvisit,
                            distr = "pois",
                            method = "mme")

# 2.5. Let's take a look at the mean estimates for each of the methods
numvisit.MLE.fit$estimate
numvisit.MME.fit$estimate
numvisit.mean
# We notice that the two estimates are extremely close to each other, and to the
# actual sample mean. We choose to work with the MLE estimate of the mean.
estimatedMean <- numvisit.MLE.fit$estimate



# ##### ---------------------------------------------------------------------- #
# 3. 
# GENERATION OF 500 SAMPLES FROM FITTED POISSON DISTRIBUTION 
# ---------------------------------------------------------------------------- #

# 3.1. Generate 500 samples from a Poisson distribution, having a mean equal to
# our MLE estimated mean.
listOfSamples <- list()      # list to hold our generated samples
numberOfSamples <- 500       # number of samples to generate 
set.seed(100)                # in order to allow for code replication
# Generate the samples
for (i in 1:numberOfSamples) {
  sample.i <- rpois(nrow(mdvis), numvisit.MLE.fit$estimate)
  listOfSamples[[i]] <- sample.i
}



# ##### ---------------------------------------------------------------------- #
# 4. 
# VARIANCE-STABILIZING CONFIDENCE INTERVAL
# ---------------------------------------------------------------------------- #

# For this section, the name of all results will be of the form "vs.[Result]"
# where "vs" refers to "variance-stablizing"

# 4.1. We create a data frame in which we save the lower and upper bounds, and 
# the length of each confidence interval, calculated from each of the 500 
# samples. We also save a binary variable "containsMean" which takes the value
# 1 if the true mean is in the interval and 0 if not.
vs.confidenceIntervals <- data.frame(
  "LB" = rep(0, numberOfSamples),
  "UB" = rep(0, numberOfSamples),
  "length" = rep(0, numberOfSamples),
  "containsMean" = rep(0, numberOfSamples)
)

# 4.2. Get Z_{alpha/2} and save the number of observations in the dataset
z <- qnorm(0.025)
n <- nrow(mdvis)

# 4.3. For each of the generated samples, calculate the variance-stablizing 
# confidence interval, as specified in (Barker, 2002, p. 87).
for (i in 1:numberOfSamples) {
  
  # Compute the mean of the given sample
  xBar <- mean(listOfSamples[[i]])
 
  # Compute the Upper Bound and Lower Bound
  lowerBound <- xBar + (z^2) / (4*n) + z * sqrt(xBar / n)
  upperBound <- xBar + (z^2) / (4*n) - z * sqrt(xBar / n)
  
  # Save the Upper Bound and Lower Bound in the created data frame
  vs.confidenceIntervals$LB[i] <- lowerBound
  vs.confidenceIntervals$UB[i] <- upperBound
  
  # Save the length of the calculated interval
  vs.confidenceIntervals$length[i] <- upperBound - lowerBound
  
  # Check whether the CI contains the mean or not
  vs.confidenceIntervals$containsMean[i] <- ifelse(
    numvisit.MLE.fit$estimate >= lowerBound &
    numvisit.MLE.fit$estimate <= upperBound,
    1, 0
  )
}

# 4.4. Compute the Empirical Coverage Probability defined as the empirical mean
# of the number of times the true mean was found in the variance-stablizing CI.
# ECP -> Empirical Coverage Probability
vs.ECP <- mean(vs.confidenceIntervals$containsMean)

# 4.5. Compute the Empirical Expected Length of the variance-stablizing CI which
# is the empirical mean of the lengths of all obtained CIs.
# EEL -> Empirical Expected Length
vs.EEL <- mean(vs.confidenceIntervals$length)



# ##### ---------------------------------------------------------------------- #
# 5. 
# BOOTSTRAP INTERVALS
# ---------------------------------------------------------------------------- #

# For each of the 500 previously generated samples from Poisson distribution, 
# generate B bootstrap samples.
# Review last slide in the slides set "Confidence Intervals" by Prof. Thomas

# 5.1. This will be a "list of lists" containing, for each of the 500 samples we
# have, B bootstrap samples
bootstrapSamples <- list()

# 5.2. This will be a "list of vectors" containing, for each of the 500 samples we 
# have, B means for each of the B bootstrap samples generated.
bootstrapMeans <- list()

# 5.3. We create a data frame in which we save the lower and upper bounds, and 
# the length of each confidence interval, calculated from each of the 500 
# samples. We also save a binary variable "containsMean" which takes the value
# 1 if the true mean is in the interval and 0 if not.
boot.confidenceIntervals <- data.frame(
  "LB" = rep(0, numberOfSamples),
  "UB" = rep(0, numberOfSamples),
  "length" = rep(0, numberOfSamples),
  "containsMean" = rep(0, numberOfSamples)
)

# 5.4. Set number of bootstrap samples to generate
B <- 100

# 5.5. For each of the 500 generated samples, generate B bootstrap samples and 
# calculate Bootstrap CIs as specified in slides by Prof. Thomas
for (i in 1:numberOfSamples) {
  
  # Initialise the list of bootstrap samples
  bootstrapSamples[[i]] <- list()
  # Initialise the vector of bootstrap means
  bootstrapMeans[[i]] <- rep(0, B)
  
  for (j in 1:B) {
    # Generate random indices of elements to include in bootstrap samples
    indices <- sample(1:n, n, replace = TRUE)
    # Generate bootstrap sample and save
    bootstrapSamples[[i]][[j]] <- mdvis$numvisit[indices]
    # Compute the mean of the generated bootstrap sample
    bootstrapMeans[[i]][j] <- mean(bootstrapSamples[[i]][[j]])
    # Compute the pivotal quantity, as specified in the slides
    bootstrapMeans[[i]][j] <- bootstrapMeans[[i]][j] - numvisit.MLE.fit$estimate
  }
  
  # Compute lower and upper bounds of bootstrap confidence interval for each 
  # of the samples as specified in the slides
  lowerBound <- numvisit.mean - quantile(bootstrapMeans[[i]], 0.975)
  upperBound <- numvisit.mean - quantile(bootstrapMeans[[i]], 0.025)
  
  # Save the Upper Bound and Lower Bound in the created data frame
  boot.confidenceIntervals$LB <- lowerBound
  boot.confidenceIntervals$UB <- upperBound
  
  # Save the length of the calculated interval
  boot.confidenceIntervals$length <- upperBound - lowerBound
  
  # Check whether the CI contains the mean or not
  boot.confidenceIntervals$containsMean <- ifelse(
    numvisit.MLE.fit$estimate >= lowerBound &
      numvisit.MLE.fit$estimate <= upperBound,
    1, 0
  )
  
}

# Let's take a look at the obtained results
View(boot.confidenceIntervals)
# Nice!

# 5.6. Compute the Empirical Coverage Probability defined as the empirical mean
# of the number of times the true mean was found in the Bootstrap CI.
# ECP -> Empirical Coverage Probability
boot.ECP <- mean(boot.confidenceIntervals$containsMean)

# 4.5. Compute the Empirical Expected Length of the variance-stablizing CI which
# is the empirical mean of the lengths of all obtained CIs.
# EEL -> Empirical Expected Length
boot.EEL <- mean(boot.confidenceIntervals$length)


# ---------------------------------------------------------------------------- #


# ????? THIS SECTION IS QUESTIONABLE 
# Looking at population coverage probability
# This method needs to be approximated numerically. Look at answer by Prof. 
# Thomas on Moodle
popcovs = 0
for (i in 0:100) {
  
  lower_bound <- i + (z^2) / (4*n) + z * sqrt(x_bar / n)
  upper_bound <- i + (z^2) / (4*n) - z * sqrt(x_bar / n)
  
  indicator <- ifelse(
    numvisit.MLE.fit$estimate >= lower_bound &
      numvisit.MLE.fit$estimate <= upper_bound,
    1, 0
  )
  
  popcovs <- popcovs + indicator * exp(
    -numvisit.MLE.fit$estimate * n 
  ) * i^(n * numvisit.MLE.fit$estimate)
  
}
