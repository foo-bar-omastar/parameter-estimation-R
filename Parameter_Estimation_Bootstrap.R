#We load the dataset using the 'import dataset' option in RStudio and then perform 
#all our calculations and solutions

gamma.arrivals <- read.table("S:/dataset/ASCII Comma/Chapter 8/gamma-arrivals.txt", quote="\"", comment.char="")

#---------------------------------------------------------------------------------

# PART A 
#-------------
# Make a histogram of the interarrival times. Does it appear that a gamma
# distribution would be a plausible model?

#To plot the estimates, we use the hist() function as below
#Frequency Histogram of the estimates
hist(gamma.arrivals$V1,
     main="Frequency Histogram: Interarrival Times",
     xlab="Interarrival Times",
     col="skyblue3",
     labels = TRUE,
     breaks = seq(min(gamma.arrivals$V1), max(gamma.arrivals$V1), length.out = 11))

#We can use the below code snippet to highlight the intervals of the interarrival
#times that were used to classify the 10 intervals of the histogram.
hist(gamma.arrivals$V1,
     breaks = seq(min(gamma.arrivals$V1), max(gamma.arrivals$V1), length.out = 11),
     plot = FALSE)



#---------------------------------------------------------------------------------

# PART B
#-------------
# Let's assume that the interarrival times are gamma distributed.


#-----------------------------------------------------------
#Part (i) - Fit the parameters of the gamma distribution by the method of moments
#-----------------------------------------------------------
# In order to calculate the parameters alpha_hat and lambda_hat, we need to
# calculate the mean and variance of all the datapoints available to us. Then
# we substitute these calculated values into formulae to get our answer

#Calculating the mean of all the data available to us
x_bar <- mean(gamma.arrivals$V1)

#Calculating the variance of all the data available to us
sum_xbar<-0
for(i in 1:nrow(gamma.arrivals)) {
        sum_xbar = sum_xbar + (gamma.arrivals$V1[i] - x_bar)^2
}

var_x <- (1/nrow(gamma.arrivals))*sum_xbar


# Now we put the above values in the following formulae to calculate alpha_hat
# and lambda_hat
alpha_hat <- (x_bar^2)/var_x
lambda_hat <- (x_bar)/var_x





#-----------------------------------------------------------
#Part (ii) - Fit the parameters of the gamma distribution by maximum
#likelihood. How do these values compare to those found before?
#-----------------------------------------------------------
# In order to calculate the parameters alpha_tilde and lambda_tilde,
# for the MLE Method, we need to leverage the value of alpha_hat 
# obtained previously in the method of moments method. We use that
# in conjunction with a formula listed below and calculate the roots
# of the equation to get the value of alpha_tilde. We will then use
# that value to get lambda_tilde as well.


#Calculating the (1/n)*sum log X_i term of the non linear equation
# that calculates the MLE of alpha

sum_log_xi <- 0
for(i in 1:nrow(gamma.arrivals)) {
        sum_log_xi = sum_log_xi + log(gamma.arrivals$V1[i])
}

sum_log_xi = sum_log_xi/nrow(gamma.arrivals)

#Calculating the log x_bar term of the non linear equation
# that calculates the MLE of alpha

log_x_bar <- log(x_bar)

#Now we could use these two values in the non linear equation and 
# solve the equation in order to get the value of alpha.
# OR We can also calculate the value of alpha using the fitdistr() function
# and use that instead. 

#Formula to fit the gamma distribuition based on MLE and yield alpha_tilde
fitdistr(gamma.arrivals$V1, "gamma")$estimate[1]

#We obtain the value of alpha_tilde, we load it into a variable,
# and then we calculate lambda_tilde as below

alpha_tilde <- fitdistr(gamma.arrivals$V1, "gamma")$estimate[1]
lambda_tilde <- alpha_tilde/x_bar




#-----------------------------------------------------------
#Part (iii) - Plot the two fitted gamma densities on top of the histogram. 
# Do the fits look reasonable?
#-----------------------------------------------------------
# We plot the fitted gamma densities on two separate histograms using the 
# values we have calculated in the previous exercises. 

#Load data-set into the variable X
X<-gamma.arrivals$V1

#Using the values from Method of Moments Estimate Exercise
den.x <- seq(0,max(X))
den.y <- dgamma(den.x, 1.012352, 0.01266466)

hist(X, col="skyblue3", probability=TRUE, ylim = c(0,1.1*max(den.y)),
     breaks = seq(min(gamma.arrivals$V1), max(gamma.arrivals$V1), length.out = 11),
     main="Method of Moments: Gamma Density vs Histogram",
     xlab="Interarrival Times")
lines(den.x, den.y, col="green", lwd=2)


#Using the values from Maximum Likelihood Estimate Exercise
den.x <- seq(0,max(X))
den.y <- dgamma(den.x, 1.02633, 0.01283952)

hist(X, col="skyblue3", probability=TRUE, ylim = c(0,1.1*max(den.y)),
     breaks = seq(min(gamma.arrivals$V1), max(gamma.arrivals$V1), length.out = 11),
     main="Maximum Likelihood Estimate: Gamma Density vs Histogram",
     xlab="Interarrival Times")
lines(den.x, den.y, col="tomato2", lwd=2)







#---------------------------------------------------------------------------------

# PART C
#-------------
#In this exercise we will estimate the standard errors of the parameter estimates


#-----------------------------------------------------------
#Part (i) - Estimate the standard errors of the parameters fit by method of 
#moments by using bootstrap.
#-----------------------------------------------------------
# RStudio makes it easy to perform parametric bootstrap sampling by using the
# boot() function. First, we create a function that will be used to calculate
# the required statistic, which is the Method of Moment Estimate in this case.


# The below function takes the dataset and creates a gamma sample that uses the
# parameter estimates that were calculated in the previous exercise. The fitdist
# function then fits a gamma distribution using the method of moments estimate
# and then returns the alpha and lambda parameter estimates for a specific sample

calcMME = function(data,sample){
        data <- rgamma(data, alpha_hat, lambda_hat)
        temp <-fitdist(data, "gamma", method = "mme")
        return(temp$estimate)
}

# Creating a boot class variable to store the values of 1000 bootstrap simulations
MME_results<-0
# Boot function to perform parametric bootstrap for 1000 samples
MME_results = boot(data = gamma.arrivals$V1, statistic = calcMME, R = 1000, sim = "parametric")
# Print the Results of the Parametric Bootstrap calculation
MME_results



#-----------------------------------------------------------
#Part (ii) - Estimate the standard errors of the parameters fit by maximum
# likelihood by using bootstrap. How do they compare to the results found previously?
#-----------------------------------------------------------


# The below function takes the dataset and creates a gamma sample that uses the
# parameter estimates that were calculated in the previous exercise. The fitdist
# function then fits a gamma distribution using the maximum likelihood estimate
# and then returns the alpha and lambda parameter estimates for a specific sample

calcMLE = function(data,sample){
        data <- rgamma(data, alpha_tilde, lambda_tilde)
        temp <-fitdist(data, "gamma", method = "mle")
        return(temp$estimate)
}


# Creating a boot class variable to store the values of 1000 bootstrap simulations
MLE_results<-0
# Boot function to perform parametric bootstrap for 1000 samples
MLE_results = boot(data = gamma.arrivals$V1, statistic = calcMLE, R = 1000, sim = "parametric")
# Print the Results of the Parametric Bootstrap calculation
MLE_results






#---------------------------------------------------------------------------------

# PART D
#-------------
#In this exercise we will create approximate confidence intervals for the 
# parameter estimates.


#-----------------------------------------------------------
#Part (i) - Use bootstrap to form 95% approximate confidence intervals for the
#parameter estimates obtained by the method of moments.
#-----------------------------------------------------------
# RStudio makes it easy to extract the confidence intervals of a parametric
# bootstrap sampling by using the boot.ci function. Here we will leverage the 
# bootstrap calculation made already in the previous exercise and use those to 
# calculate the confidence intervals. 


# Boot.ci function to calculate the confidence intervals for the bootstraps 
# created in the previous exercise. We specify the 95% confidence interval,
# set index = 1 to select the first output of the boot function which is alpha,
# then we set index=2 for the confidence interval for lambda

MME_CI_alpha = boot.ci(MME_results, conf=0.95, type ="basic", index = 1)
MME_CI_alpha

MME_CI_lambda = boot.ci(MME_results, conf=0.95, type ="basic", index = 2)
MME_CI_lambda



#-----------------------------------------------------------
#Part (ii) - Use bootstrap to form 95% approximate confidence intervals for 
# the parameter estimates obtained by maximum likelihood. How do the 
# confidence intervals for the two methods compare?
#-----------------------------------------------------------

# Boot.ci function to calculate the confidence intervals for the bootstraps 
# created in the previous exercise. We specify the 95% confidence interval,
# set index = 1 to select the first output of the boot function which is alpha,
# then we set index=2 for the confidence interval for lambda

MLE_CI_alpha = boot.ci(MLE_results, conf=0.95, type ="basic", index = 1)
MLE_CI_alpha

MLE_CI_lambda = boot.ci(MLE_results, conf=0.95, type ="basic", index = 2)
MLE_CI_lambda
