
rm(list = ls())

source('black_scholes_formula.R')
library(quantmod)
library(stats4) # mle function 


## stock prices from yahoo using the package quantmod

# get SPY stock prices
setSymbolLookup(SPY='yahoo')
getSymbols('SPY')

# convert data from xts table to data frame
SPY <- as.data.frame(SPY) 
SPY <- cbind.data.frame(as.Date(row.names(SPY)), SPY$SPY.Close)
rownames(SPY) <- NULL
colnames(SPY) <- c('date', 'price')
head(SPY)
str(SPY)

# subset data

SPY_sub <- subset(SPY, date >= as.Date('2012-10-17'))
plot(1:dim(SPY_sub)[1], SPY_sub$price, type = 'l')


## log returns

log_return <- log(SPY_sub$price[-1] / SPY_sub$price[-dim(SPY_sub)[1]])
hist(log_return)


## test independence of log returns - seems OK

acf(log_return)
Box.test(log_return, lag = 1, type = "Ljung-Box")

## test normal distribution - not OK in general! Heavier tales than a normal distribution

qqnorm(log_return);qqline(log_return, col = 2)
shapiro.test(log_return)

## estimate Q-parameters in GBM - assuming equidistant data

sigma_mle <-0.1
r_mle <- 0
dt <- 1/200

# dt <- 1/(dim(SPY_sub)[1]/ 5) # one over the number of observations per year
# sigma_mle <- sqrt(var(log_return) / dt)
# r_mle <- mean(log_return) / dt + sigma_mle / 2

hist(log_return, freq = FALSE)
lines(sort(log_return), dnorm(sort(log_return), mean = (r_mle - sigma_mle / 2) * dt, sd = sigma_mle * sqrt(dt)), lwd = 2)

# # approach with likelihood function - based on normality assumption
# 
# LL <- function(obs, mu, sigma) {
#    -sum(log(dnorm(obs, mu, sigma)))
# }
# 
# mle_estimates <- mle(function(mu, sigma) LL(obs = log_return, mu, sigma),
#                      start = list(mu = 1, sigma=1))
# 
# r_mle <- 0.0004428211
# sigme_mle <- 0.0076845548
# 
# sigme_mle <- sqrt(0.0076845548^2 / dt)
# r_mle <- 0.0004428211 / dt + sigma_mle / 2
# 
# LL <- function(obs, theta) {
#   -sum(log(dnorm(obs, theta[1], theta[2])))
# } 
# 
# nlm(function(theta) LL(obs = log_return, theta), c(1, 1), hessian = TRUE)
# 
# optim(c(0.1, 0.1), 
#       function(theta) LL(obs = log_return, theta),
#       method = "L-BFGS-B",
#       lower = c(0.001, 0.001),
#       upper = c(1, 1))
# 
# LL(obs = log_return, theta = c(0.001, 0.001))

## option data from Yahoo using the package quantmod

spy_options <- getOptionChain('SPY') # SPY is one of the most liquid option markets - each options contract gives the owner the right to 100 shares of the underlying ETF (exchange-traded fund)
head(spy_options$calls)
str(spy_options)

# call option data - set price to average between bid and ask
spy_options_call <- cbind.data.frame(rownames(spy_options$calls), 
                                     spy_options$calls$Strike, 
                                     rowMeans(cbind(spy_options$calls$Bid, spy_options$calls$Ask)))
colnames(spy_options_call) <- c('symbol', 'strike', 'price')

# add column maturity_date maturity (in years)
spy_options_call$maturity_date <- 
  as.Date(substring(spy_options_call$symbol, 4, 9), format = '%d%m%y')

# add column maturity (in years)
spy_options_call$maturity <- as.numeric((spy_options_call$maturity_date - max(SPY_sub$date))/365.25)

price_bs <- rep(NA, dim(spy_options_call)[1])
for (k in 1:dim(spy_options_call)[1]) {

  price_underlying <- SPY_sub$price[SPY_sub$date == max(SPY_sub$date)]
  maturity <- spy_options_call$maturity[k]
  strike <- spy_options_call$strike[k]
  
  # black-scholes call price
  price_bs[k] <- black_scholes_formula(t = 0, 
                                    maturity = maturity,
                                    s = price_underlying, 
                                    K = strike, 
                                    r = r_mle, 
                                    sigma = sigma_mle)
  
}


plot(price_bs, type = 'l')
points(spy_options_call$price)
# plot((price_bs - spy_options_call$price) / price_bs)
