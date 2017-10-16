rm(list=ls()) # clears all

library(ggplot2)
library(reshape2)

#### GBM simulation ####

## initialize

x0    <- 100
mu    <- 0.01
sigma <- 0.1


## simulate short rate paths

n  <- 10^3    # monte carlo simulation paths 
endTime  <- 10
dt <- 1 / 12
time_vector <- seq(0, endTime, by = dt)
m  <- length(time_vector)      # subintervals



x <- matrix(NA, m, n)
x[1,] <- x0


for (k in 1:n) {
  for (i in 2:m) {
    dx <- mu * x[i - 1, k] * dt + x[i - 1, k] * sigma * sqrt(dt) * rnorm(1,0,1)
    x[i, k] <- x[i - 1, k] + dx
    # x[i, k] <- x0 * exp((mu - 0.5 * sigma ^ 2) * time_vector[i] + sigma * sqrt(time_vector[i]) * rnorm(1, 0, 1))
  }
}



# ## plot
# 
# data <- data.frame('time' = seq(0, m), x)
# data <- melt(data,  id = c('time'))
# 
# ggplot(data, aes(time, value)) +
#   geom_line(aes(colour = variable))

## present value of future expected cash flow

r <- 0.01
price_MC <- rowMeans(1/(1 + r * dt) ^ (0:(m - 1))  * x)
plot(time_vector, price_MC)

## present value of expected future profits

asset_expected <- rowMeans(x)
liability_expected <- 100
profit_expected <- (asset_expected - liability_expected)
pv_profit <- sum(1/(1 + r * dt) ^ (0:(m - 1)) * profit_expected)

## thoeretical expectation and variance of GBM

expectation_theoretical <- x0 * exp(mu * time_vector)
variance_theoretical <- x0^2 * exp(2 * mu * time_vector) * (exp(sigma ^ 2 * time_vector) - 1)
std_theoretical <- sqrt(variance_theoretical)

## compare empirical and theoretical moments

plot(time_vector, asset_expected) 
points(time_vector, expectation_theoretical, type = 'l')

plot(time_vector, apply(x, 1, var)) 
points(time_vector, variance_theoretical, type = 'l')

## compare histogram and thoeretical distribution of GBM

density_theoretical_GBM <- function(s, t, mu, sigma, s0) {
  
  1 / sqrt(2 * pi * s ^ 2 * sigma ^ 2 * t) * 
    exp(
      -((log(s) - log(s0) - (mu - 0.5 * sigma ^ 2) * t) ^ 2 / (2 * sigma ^ 2 * t))
    )
  
}

density <- sapply(sort(x[2,]), function(i) density_theoretical_GBM(s = i, 
                                                                   t = time_vector[2], 
                                                                   mu = mu, 
                                                                   sigma = sigma,
                                                                   s0 = x0))
hist(x[2,], freq = F)
lines(sort(x[2,]), density, lwd = 2)



## confidence interval - forklar hvorfor det vokser med tiden

# sapply(1:m, function(i) quantile(x[i,], c(0.025,0.975)) ) 


## does the simulation method matter? YES! x_closed_form is wrong because W(t_i) and W(t_{i-1}) are not independent. But in the simulation scheme below we assume independence through the independent normal rv's, and hence the resulting path for x_closed_form is incorrect 

n  <- 1    # monte carlo simulation paths 
endTime  <- 10
dt <- 1 / 12
time_vector <- seq(0, endTime, by = dt)
m  <- length(time_vector)     

x_closed_form <- x_increment <- x_closed_form_increment <- matrix(NA, m, n)
x_closed_form[1, ] <- x_increment[1, ] <- x_closed_form_increment[1, ] <- x0

for (k in 1:n) {
  for (i in 2:m) {
    Z <- rnorm(1,0,1)
    dx <- mu * x_increment[i - 1, k] * dt + x_increment[i - 1, k] * sigma * sqrt(dt) * Z
    x_increment[i, k] <- x_increment[i - 1, k] + dx
    x_closed_form[i, k] <- x0 * exp((mu - 0.5 * sigma ^ 2) * time_vector[i] + sigma * sqrt(time_vector[i]) * rnorm(1, 0, 1))
    x_closed_form_increment[i, k] <- x_closed_form_increment[i - 1, k] * exp((mu - 0.5 * sigma ^ 2) * dt + sigma * sqrt(dt) * Z)
  }
}

plot(time_vector, x_closed_form[, 1], type = 'l')
points(time_vector, x_increment[, 1], type = 'l', col = 'red')
points(time_vector, x_closed_form_increment[, 1], type = 'l', col = 'blue')


