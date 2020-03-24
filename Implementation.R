# ----------------------- MATHEMTATICAL STATISTICS II ------------------------ #
# ------------ Project 2: Confidence Interval for a Poisson Mean ------------- #
# ---------------------------------------------------------------------------- #
# ------- Yanxin LI ---- Paul MELKI ---- Mira RAHAL ---- Zefeng ZHANG -------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Created on 08/03/2020 ------------------------------------------------------ #
# Iteration 4.2.1 ---------- Last updated on 24/03/2020 -------- by Paul MELKI #
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

# 1.3. Set working directory, in order to save plots, later on.
setwd("P:\\College Material\\Semester 2\\Mathematical Statistics 2\\Project\\Code")



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

# 3.1. Save the number of observations in the dataset
n <- nrow(mdvis)

# 3.2. Generate 500 samples from a Poisson distribution, having a mean equal to
# our MLE estimated mean.
listOfSamples <- list()      # list to hold our generated samples
numberOfSamples <- 500       # number of samples to generate 
set.seed(100)                # in order to allow for code replication
# Generate the samples
for (i in 1:numberOfSamples) {
  sample.i <- rpois(n, numvisit.MLE.fit$estimate)
  listOfSamples[[i]] <- sample.i
}



# ##### ---------------------------------------------------------------------- #
# 4. 
# WALD CONFIDENCE INTERVAL
# ---------------------------------------------------------------------------- #

# For this section, the names of all results will be of the form "wald.[Result]"

# 4.1. We create a data frame in which we save the lower and upper bounds, and 
# the length of each confidence interval, calculated from each of the 500 
# samples. We also save a binary variable "containsMean" which takes the value
# 1 if the true mean is in the interval and 0 if not.

wald.confidenceIntervals <- data.frame(
  "LB" = rep(0, numberOfSamples),
  "UB" = rep(0, numberOfSamples),
  "length" = rep(0, numberOfSamples),
  "containsMean" = rep(0, numberOfSamples)
)

# 4.2. Get the required z_{\alpha_2} value
z <- abs(qnorm(0.025))

# 4.3. For each of the generated samples, calculate the Wald confidence interval
# as specified in (Barker, 2002, p. 86).
for (i in 1:numberOfSamples) {
  
  # Compute the mean of the sample
  xBar <- mean(listOfSamples[[i]])
  
  # Compute the upper and lower bounds
  upperBound <- xBar + z * (xBar/n)^(0.5)
  lowerBound <- xBar - z * (xBar/n)^(0.5)
  
  # Save the upper and lower bounds in the created dataframe
  wald.confidenceIntervals$UB[i] <- upperBound
  wald.confidenceIntervals$LB[i] <- lowerBound
  
  # Compute the length of the confidence interval
  wald.confidenceIntervals$length[i] <- upperBound - lowerBound
  
  # Check whether the CI contains mean or not
  wald.confidenceIntervals$containsMean[i] <- ifelse(
    numvisit.MLE.fit$estimate >= lowerBound &
      numvisit.MLE.fit$estimate <= upperBound,
    1, 0
  )
}

# 4.4. Compute the Empirical Coverage Probability defined as the empirical mean
# of the number of times the true mean was found in the Wald CI.
# ECP -> Empirical Coverage Probability
wald.ECP <- mean(wald.confidenceIntervals$containsMean)

# 4.5. Compute the Empirical Expected Length of the Wald CI which is the
# empirical mean of the lengths of all obtained CIs.
# EEL -> Empirical Expected Length
wald.EEL <- mean(wald.confidenceIntervals$length)

# 4.6. Compute the numerical approximation of the population coverage 
# probability (APCP) as defined in (Barker, 2002, p. 85)
wald.APCP <- 0
# We choose the value of 6000 in order to keep the computations tractable...
for (i in 1:6000) {
  
  # Compute the upper and lower bounds
  lowerBound <- i - z * (i^0.5)
  upperBound <- i + z * (i^0.5)
  
  # Create create an indicator of the event n*\theta between upper & lower bound
  indic <- ifelse(n * numvisit.MLE.fit$estimate >= lowerBound &
                  n * numvisit.MLE.fit$estimate <= upperBound, 1, 0)
  
  wald.APCP <- wald.APCP + indic * dpois(i, n * numvisit.MLE.fit$estimate)
}

# Let's take a look at the obtained approximation of the population coverage
# probability
wald.APCP # Nice!



# ##### ---------------------------------------------------------------------- #
# 5. 
# SCORE CONFIDENCE INTERVAL
# ---------------------------------------------------------------------------- #

# In this section, the names of all results will be of the form "score.[Result]"

# 5.1. We create a data frame in which we save the lower and upper bounds, and 
# the length of each confidence interval, calculated from each of the 500 
# samples. We also save a binary variable "containsMean" which takes the value
# 1 if the true mean is in the interval and 0 if not.

score.confidenceIntervals <- data.frame(
  "LB" = rep(0, numberOfSamples),
  "UB" = rep(0, numberOfSamples),
  "length" = rep(0, numberOfSamples),
  "containsMean" = rep(0, numberOfSamples)
)

# 5.2. For each of the generated samples, calculate the Score confidence 
# interval as specified in (Barker, 2002, p. 86).

for (i in 1:numberOfSamples) {
  
  # Compute the mean of the sample
  xBar <- mean(listOfSamples[[i]])
  print(xBar)
  
  # Compute the upper and lower bounds
  upperBound <- xBar + z^2/(2 * n)- z * sqrt(4 * xBar + (z^2)/n)/(sqrt(4 * n))
  lowerBound <- xBar + z^2/(2 * n)+ z * sqrt(4 * xBar + (z^2)/n)/(sqrt(4 * n))
  print(paste(upperBound, lowerBound))
  
  # Save the upper and lower bounds in the created dataframe
  score.confidenceIntervals$UB[i] <- upperBound
  score.confidenceIntervals$LB[i] <- lowerBound
  
  # Compute the length of the confidence interval
  score.confidenceIntervals$length[i] <- upperBound - lowerBound
  
  # Check whether the CI contains mean or not
  score.confidenceIntervals$containsMean[i] <- ifelse(
    numvisit.MLE.fit$estimate >= lowerBound &
      numvisit.MLE.fit$estimate <= upperBound,
    1, 0
  )
}

# 5.3. Compute the Empirical Coverage Probability defined as the empirical mean
# of the number of times the true mean was found in the Score CI.
# ECP -> Empirical Coverage Probability
score.ECP <- mean(score.confidenceIntervals$containsMean)

# 5.4. Compute the Empirical Expected Length of the Score CI which is the
# empirical mean of the lengths of all obtained CIs.
# EEL -> Empirical Expected Length
score.EEL <- mean(score.confidenceIntervals$length)

# 5.5. Compute the numerical approximation of the population coverage 
# probability (APCP) as defined in (Barker, 2002, p. 85)
lowerBound <- NULL
upperBound <- NULL
score.APCP <- 0
# We choose the value of 6000 in order to keep the computations tractable...
for (i in 1:6000) {
  
  # Compute the upper and lower bounds
  upperBound <- i/n + z^2/(2 * n) + z * sqrt(4 * i/n + (z^2)/n)/(sqrt(4 * n))
  lowerBound <- i/n + z^2/(2 * n) - z * sqrt(4 * i/n + (z^2)/n)/(sqrt(4 * n))

  # Create create an indicator of the event n*\theta between upper & lower bound
  indic <- ifelse(numvisit.MLE.fit$estimate >= lowerBound &
                  numvisit.MLE.fit$estimate <= upperBound, 1, 0)
  
  score.APCP <- score.APCP + indic * dpois(i, n * numvisit.MLE.fit$estimate)
}

# Let's take a look at the obtained approximation of the population coverage
# probability
score.APCP # Nice!



# ##### ---------------------------------------------------------------------- #
# 6. 
# VARIANCE-STABILIZING CONFIDENCE INTERVAL
# ---------------------------------------------------------------------------- #

# For this section, the name of all results will be of the form "vs.[Result]"
# where "vs" refers to "variance-stablizing"

# 6.1. We create a data frame in which we save the lower and upper bounds, and 
# the length of each confidence interval, calculated from each of the 500 
# samples. We also save a binary variable "containsMean" which takes the value
# 1 if the true mean is in the interval and 0 if not.
vs.confidenceIntervals <- data.frame(
  "LB" = rep(0, numberOfSamples),
  "UB" = rep(0, numberOfSamples),
  "length" = rep(0, numberOfSamples),
  "containsMean" = rep(0, numberOfSamples)
)

# 6.2. Get Z_{alpha/2}
z <- qnorm(0.025)

# 6.3. For each of the generated samples, calculate the variance-stablizing 
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

# 6.4. Compute the Empirical Coverage Probability defined as the empirical mean
# of the number of times the true mean was found in the variance-stablizing CI.
# ECP -> Empirical Coverage Probability
vs.ECP <- mean(vs.confidenceIntervals$containsMean)

# 6.5. Compute the Empirical Expected Length of the variance-stablizing CI which
# is the empirical mean of the lengths of all obtained CIs.
# EEL -> Empirical Expected Length
vs.EEL <- mean(vs.confidenceIntervals$length)

# 6.6. Compute the numerical approximation of the population coverage 
# probability (APCP) as defined in (Barker, 2002, p. 85)
vs.APCP <- 0
# We choose the value of 6000 in order to keep the computations tractable...
for (i in 1:6000) {
  
  # Compute the upper and lower bounds
  upperBound <- i/n + (z^2) / (4*n) + z * sqrt(i) / n
  lowerBound <- i/n + (z^2) / (4*n) - z * sqrt(i) / n

  # Create create an indicator of the event n*\theta between upper & lower bound
  indic <- ifelse(numvisit.MLE.fit$estimate >= lowerBound &
                  numvisit.MLE.fit$estimate <= upperBound, 1, 0)
  # Update value
  vs.APCP <- vs.APCP + indic * dpois(i, n * numvisit.MLE.fit$estimate)
}

# Let's take a look at the obtained approximation of the population coverage
# probability
vs.APCP # Nice!



# ##### ---------------------------------------------------------------------- #
# 7. 
# EXACT METHOD CONFIDENCE INTERVAL
# ---------------------------------------------------------------------------- #

# For this section, the name of all results will be of the form "exact.[Result]"

# Since we know that the value of the estimated lambda paremeter is 2.589133, 
# we start with a list of values of lambda between 2 and 3, in order to find the
# lower and upper bounds using the Exact Method. 

# 7.1. Initialise the list of lambdas and save its length to be used later on
listOfLambdas <- seq(from = 2.0, to = 3.0, by = 0.001)
numberOfLambdas <- length(listOfLambdas)

# 7.2. We create a data frame in which we save the lower and upper bounds, and 
# the length of each confidence interval, calculated from each of the 500 
# samples. We also save a binary variable "containsMean" which takes the value
# 1 if the true mean is in the interval and 0 if not.
exact.confidenceIntervals <- data.frame(
  "LB" = rep(0, numberOfSamples),
  "UB" = rep(0, numberOfSamples),
  "length" = rep(0, numberOfSamples),
  "containsMean" = rep(0, numberOfSamples)
)

# 7.3. For each of the generated samples, calculate the "exact" confidence 
# interval, as specified in (Barker, 2002, p. 88).
for (i in 1:numberOfSamples) {
  
  # In order to calculate the upper and lower bounds of the Exact CI, we need 
  # Poisson CDFs for each of the lambdas in our listOfLambdas. We create the 
  # vectors that will carry these values
  upperBoundCDF <- lowerBoundCDF <- rep(0, numberOfLambdas)
  
  # We defined the T_observed, as defined in (Barker, 2002, p.86), as the sum
  # of observations in the sample
  tObs <- sum(listOfSamples[[i]])
  
  # Calculate the lower and upper bound Poisson CDFs
  for (j in 1:numberOfLambdas) {
    upperBoundCDF[j] <- ppois(tObs, lambda = listOfLambdas[j] * n)
    lowerBoundCDF[j] <- ppois(tObs, lambda = listOfLambdas[j] * n, lower = F)
  }
  
  # Taking alpha = 5%, we choose all those values in the upperBoundCDF and
  # lowerBoundCDF that are less than alpha / 2
  upperBoundCDF.requiredIndices <- which(upperBoundCDF <= 0.025)
  lowerBoundCDF.requiredIndices <- which(lowerBoundCDF <= 0.025)
  
  # As defined, the Exact method's upper bound is the minimum value among the 
  # lambdas whose indices are in upperBoundCDF.requiredIndices
  exact.confidenceIntervals$UB[i] <- min(
    listOfLambdas[upperBoundCDF.requiredIndices]
  )
  # and the lower bound is the maximum value among the lambdas whose indices are
  # in lowerBoundCDF.requiredIndices
  exact.confidenceIntervals$LB[i] <- max(
    listOfLambdas[lowerBoundCDF.requiredIndices]
  )
  
  # Calculate the length of the obtained confidence interval
  exact.confidenceIntervals$length <- 
    exact.confidenceIntervals$UB - exact.confidenceIntervals$LB
  
  # Check whether the CI contains the mean or not
  exact.confidenceIntervals$containsMean[i] <- ifelse(
    numvisit.MLE.fit$estimate >= exact.confidenceIntervals$LB[i] &
      numvisit.MLE.fit$estimate <= exact.confidenceIntervals$UB[i],
    1, 0
  )
}

# 7.4. Compute the Empirical Coverage Probability defined as the empirical mean
# of the number of times the true mean was found in the Exact CI.
# ECP -> Empirical Coverage Probability
exact.ECP <- mean(exact.confidenceIntervals$containsMean)

# 7.5. Compute the Empirical Expected Length of the Exact CI which is the
# empirical mean of the lengths of all obtained CIs.
# EEL -> Empirical Expected Length
exact.EEL <- mean(exact.confidenceIntervals$length)

# 7.6. Compute the numerical approximation of the population coverage 
# probability (APCP) as defined in (Barker, 2002, p. 85)
exact.APCP <- 0
# Create list that will hold the upper and lower bounds computed using the Exact
# method
upperBounds <- list()
lowerBounds <- list()
# We choose the value of 6000 in order to keep the computations tractable...
# We choose the value of 6000 in order to keep the computations tractable...
# Create list that will hold the upper and lower bounds calculated based on the
# Exact Method
upperBounds <- list()
lowerBounds <- list()
for (i in 1:6000) {
  
  # Calculate the upper and lower bounds following the Exact Method
  for (j in seq(from = 2, to = 3, by = 0.001)) {
    # Compute the lower bound CDF
    lowerBoundCDF <- ppois(i, n * j, lower = FALSE)
    # Save computed lower bound
    lowerBounds[i] <- n * (j - 0.001)
    # We stop when the lowerBoundCDF is smaller or equal to alpha/2
    if (lowerBoundCDF > 0.025) break
  }
  
  for (j in seq(from = 2, to = 3, by = 0.001)) {
    # Compute the upper bound CDF
    upperBoundCDF <- ppois(i, n * j)
    # Save computed upper bound
    upperBounds[i] <- n * j
    # We stop when the upperBoundCDF is greater than alpha/2
    if (upperBoundCDF <= 0.025) break
  }

  # Create create an indicator of the event n*\theta between upper & lower bound
  indic <- ifelse(n * numvisit.MLE.fit$estimate >= lowerBounds[i] &
                    n * numvisit.MLE.fit$estimate <= upperBounds[i], 1, 0)
  # Update value
    exact.APCP <- exact.APCP + indic * dpois(i, n * numvisit.MLE.fit$estimate)
}

# Let's take a look at the obtained approximation of the population coverage
# probability
exact.APCP # Nice!



# ##### ---------------------------------------------------------------------- #
# 8. 
# BOOTSTRAP CONFIDENCE INTERVAL
# ---------------------------------------------------------------------------- #

# For each of the 500 previously generated samples from Poisson distribution, 
# generate B bootstrap samples.
# Review last slide in the slides set "Confidence Intervals" by Prof. Thomas

# 8.1. This will be a "list of lists" containing, for each of the 500 samples we
# have, B bootstrap samples
bootstrapSamples <- list()

# 8.2. This will be a "list of vectors" containing, for each of the 500 samples we 
# have, B means for each of the B bootstrap samples generated.
bootstrapMeans <- list()

# 8.3. We create a data frame in which we save the lower and upper bounds, and 
# the length of each confidence interval, calculated from each of the 500 
# samples. We also save a binary variable "containsMean" which takes the value
# 1 if the true mean is in the interval and 0 if not.
boot.confidenceIntervals <- data.frame(
  "LB" = rep(0, numberOfSamples),
  "UB" = rep(0, numberOfSamples),
  "length" = rep(0, numberOfSamples),
  "containsMean" = rep(0, numberOfSamples)
)

# 8.4. Set number of bootstrap samples to generate
B <- 100

# 8.5. For each of the 500 generated samples, generate B bootstrap samples and 
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
  boot.confidenceIntervals$LB[i] <- lowerBound
  boot.confidenceIntervals$UB[i] <- upperBound
  
  # Save the length of the calculated interval
  boot.confidenceIntervals$length[i] <- upperBound - lowerBound
  
  # Check whether the CI contains the mean or not
  boot.confidenceIntervals$containsMean[i] <- ifelse(
    numvisit.MLE.fit$estimate >= lowerBound &
      numvisit.MLE.fit$estimate <= upperBound,
    1, 0
  )
  
}

# Let's take a look at the obtained results
View(boot.confidenceIntervals)
# Nice!

# 8.6. Compute the Empirical Coverage Probability defined as the empirical mean
# of the number of times the true mean was found in the Bootstrap CI.
# ECP -> Empirical Coverage Probability
boot.ECP <- mean(boot.confidenceIntervals$containsMean)

# 8.7. Compute the Empirical Expected Length of the variance-stablizing CI which
# is the empirical mean of the lengths of all obtained CIs.
# EEL -> Empirical Expected Length
boot.EEL <- mean(boot.confidenceIntervals$length)



# ##### ---------------------------------------------------------------------- #
# 9.
# COMPARISON OF BETWEEN THE DIFFERENT CONFIDENCE INTERVALS 
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# 9.1. Comparison of the CIs lengths
# ---------------------------------------------------------------------------- #

# We create a boxplot showing the distribution of the lengths of the different 
# CIs implemented above. We implement the plots using "ggplot2"

# 9.1.1. We prepare the data in the correct format for plotting
dataToPlot.lengths <- data.frame(
  "confidenceInterval" = rep(NA, 5 * numberOfSamples),
  "length" = rep(0, 5 * numberOfSamples)
)
dataToPlot.lengths$confidenceInterval <- c(
  rep("Wald", numberOfSamples), 
  rep("Score", numberOfSamples),
  rep("VS", numberOfSamples),
  rep("Exact", numberOfSamples),
  rep("Bootstrap", numberOfSamples)
)
dataToPlot.lengths$length <- c(
  wald.confidenceIntervals$length,
  score.confidenceIntervals$length,
  vs.confidenceIntervals$length,
  exact.confidenceIntervals$length,
  boot.confidenceIntervals$length
)

# 9.1.2. Plot!
ggplot(data = dataToPlot.lengths, 
       aes(x = confidenceInterval, y = length,
           fill = confidenceInterval)) + 
  geom_boxplot(notch = TRUE) + 
  scale_fill_brewer("Confidence Interval", palette = "GnBu") + 
  xlab("") + ylab("Length of CI") + 
  ggtitle("Distribution of the CIs lengths") +
  theme(plot.title = element_text(face = "bold"), 
        plot.subtitle = element_text(face = "italic"), 
        legend.text = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  ggsave("./Plots/Boxplots_CIs_All.png", type = "cairo", scale = 1.1)

# We notice that there is a huge difference between the lengths of the bootstrap
# CIs and the other CIs in such a way that it is impossible to see the 
# the difference between the other different CIs. For this reason, repeat the 
# plot again, ignoring the bootstrap CIs.

# 9.1.3. Repeat the same as above, ignoring the bootstrap CIs.
dataToPlot.lengths <- data.frame(
  "confidenceInterval" = rep(NA, 4 * numberOfSamples),
  "length" = rep(0, 4 * numberOfSamples)
)
dataToPlot.lengths$confidenceInterval <- c(
  rep("Wald", numberOfSamples),
  rep("Score", numberOfSamples), 
  rep("VS", numberOfSamples),
  rep("Exact", numberOfSamples)
)
dataToPlot.lengths$length <- c(
  wald.confidenceIntervals$length,
  score.confidenceIntervals$length,
  vs.confidenceIntervals$length,
  exact.confidenceIntervals$length
)

# 9.1.4. Plot!
ggplot(data = dataToPlot.lengths, aes(x = confidenceInterval, y = length, 
                                      fill = confidenceInterval)) + 
  geom_boxplot(notch = TRUE) + 
  scale_fill_brewer("Confidence Interval", palette = "GnBu") + 
  xlab("") + ylab("Length of CI") + 
  ggtitle("Distribution of the CIs lengths", 
          subtitle = "Values are extremely similar, the difference is barely apparent") +
  theme(plot.title = element_text(face = "bold"), 
        plot.subtitle = element_text(face = "italic"), 
        legend.text = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  ggsave("./Plots/Boxplots_CIs_noBootstrap.png", type = "cairo", scale = 1.1)
