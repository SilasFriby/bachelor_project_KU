
rm(list = ls())

source('functions.R')
source('present_value_functions.R')
library(reshape2)
library(ggplot2)

## initialize

# simulation
n  <- 10^4 # monte carlo simulation paths 
end_time  <- 10
dt <- 1/12
time_vector <- seq(0, end_time, by = dt)

# stock parameters
s0    <- 100
sigma_s <- 0.1

# short rate parameters
r0 <- 0
kappa <- 0.01
theta_P <- 0.04
theta_Q <- 0.06
sigma_r <- 0.01


## simulate short rate under Q

short_rate <- simulate_vasicek(x0 = r0,
                               n = n,
                               dt = dt,
                               end_time = end_time,
                               kappa = kappa,
                               theta = theta_P,
                               sigma = sigma_r)

## plot short rate paths

data <- data.frame('time' = time_vector, short_rate[, 1:10])
data <- melt(data,  id = c('time'))

ggplot(data, aes(time, value)) +
  geom_line(aes(colour = variable))



## simulate stock prices under U

stock_prices <- simulate_gmb(x0 = s0,
                             n = n,
                             dt = dt,
                             end_time = end_time,
                             mu = r0,#short_rate,
                             sigma = sigma_s,
                             drift = 'constant')

## plot stock price paths

data <- data.frame('time' = time_vector, stock_prices[, 1:10])
data <- melt(data,  id = c('time'))

ggplot(data, aes(time, value)) +
  geom_line(aes(colour = variable))



## firm asset

number_of_call_options <- 1100
asset_value_scenarios <- sapply(time_vector, function(i) {
  
  pv_european_option_MC(t = i,
                        time_vector = time_vector,
                        dt = dt,
                        underlying_paths = stock_prices,
                        strike = 120,
                        short_rate_paths = short_rate,
                        type = 'call')
  
})

## firm liability

number_of_non_callable_loans <- 1
liability_value_expected <- sapply(time_vector, function(i){ # VERY SLOW!
  
  pv_non_callable_bond_MC(t = i,
                          time_vector = time_vector,
                          dt = dt,
                          payment_dt = 100,
                          short_rate_paths = short_rate)
  
}) * number_of_non_callable_loans


## expected future profit

expected_future_profit <- asset_value_expected - liability_value_expected 
plot(time_vector, expected_future_profit)

## expected future profit

# yield curve today
yield_curve <- VASICEKyield(r0, time_vector[-1], c(kappa, theta_P, sigma_r), theta_Q - theta_P)
plot(time_vector[-1], yield_curve)

pv_expected_future_profit <- sum( 1/(1 + c(0, yield_curve)) * expected_future_profit )



# 
# 
# 
# 
# 
# 
# annuity <- 10^2 # payment each year
# 
# number_of_stocks <- 1
# 
# # expected scenario
# asset_value_expected <- rowMeans(stock_prices[which(time_vector %in% times_of_interest), ]) * number_of_stocks
# liability_value_expected <- loan_annuity
# profit_expected <- asset_value_expected - liability_value_expected  
# pv_profit_expected <- sum(1/(1 + yield_curve) * profit_expected)
# 
# 
# # very good scenario
# asset_value_expected <- apply(stock_prices[which(time_vector %in% times_of_interest), ], 1, FUN = function(x) quantile(x, 0.99)) * number_of_stocks
# liability_value_expected <- loan_annuity
# profit_expected <- asset_value_expected - liability_value_expected  
# pv_profit_good <- sum(1/(1 + yield_curve) * profit_expected)
# 
# # very bad scenario
# asset_value_expected <- apply(stock_prices[which(time_vector %in% times_of_interest), ], 1, FUN = function(x) quantile(x, 0.01)) * number_of_stocks
# liability_value_expected <- loan_annuity
# profit_expected <- asset_value_expected - liability_value_expected  
# pv_profit_bad <- sum(1/(1 + yield_curve) * profit_expected)
# 
# ## proit loss distribution
# 
# profit_all_scenarios <- stock_prices[which(time_vector %in% times_of_interest), ] * number_of_stocks - loan_annuity
# pv_profit_all_scenarios <- colSums(1/(1 + yield_curve) * profit_all_scenarios)
# plot(pv_profit_all_scenarios)
# hist(pv_profit_all_scenarios)
# 
# 
# ## assets in call options
# 
# strike <- 120
# number_of_call_options <- 10
# 
# # expected scenario
# asset_value_expected <- mean(pmax(stock_prices[which(time_vector == end_time),] - strike, 0)) * number_of_call_options
# liability_value_expected <- loan_annuity
# 
# profit_expected <- asset_value_expected - liability_value_expected  
# pv_profit_expected <- sum(1/(1 + yield_curve) * profit_expected)


