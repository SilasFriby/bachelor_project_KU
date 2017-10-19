
rm(list = ls())

source('black_scholes_formula.R')
library(quantmod)
library(RND)

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

SPY_sub <- subset(SPY, date >= as.Date('2016-10-17'))
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

## estimate Q-parameters in GBM - assuming equidistant data - recall drift of stock is not neccessary for pricing!

# mle
r <- 0 # the variance on the drift estimate in GBM is very high - hence I have just chosen a value based on the market interest level
dt <- 1/250 # approximate number of observations per year
sigma_mle <- sqrt(var(log_return) / dt)

# mle variance
var_r <- sigma_mle ^ 2 * (2 + sigma_mle ^ 2 * dt) / (2 * dt) 
var_sigma_mle <- sigma_mle ^ 2 /2


# dt <- 1#/(dim(SPY_sub)[1]/ 5) # one over the number of observations per year
# sigma_mle <- sqrt(var(log_return) / dt)
# r <- mean(log_return) / dt + sigma_mle ^ 2 / 2

hist(log_return, freq = FALSE)
lines(sort(log_return), dnorm(sort(log_return), mean = (r - sigma_mle ^ 2 / 2) * dt, sd = sigma_mle * sqrt(dt)), lwd = 2)


## option data from Yahoo using the package quantmod

spy_options <- getOptionChain('SPY') # SPY is one of the most liquid option markets - each options contract gives the owner the right to 100 shares of the underlying ETF (exchange-traded fund)
head(spy_options$calls)
str(spy_options)

# call option data - set price to average between bid and ask
spy_options_call <- cbind.data.frame(rownames(spy_options$calls), 
                                     spy_options$calls$Strike, 
                                     spy_options$calls$Last) #rowMeans(cbind(spy_options$calls$Bid, spy_options$calls$Ask)))
colnames(spy_options_call) <- c('symbol', 'strike', 'price')

# add column maturity_date maturity (in years)
spy_options_call$maturity_date <- 
  as.Date(substring(spy_options_call$symbol, 4, 9), format = '%d%m%y')

# add column maturity (in years)
spy_options_call$maturity <- as.numeric((spy_options_call$maturity_date - max(SPY_sub$date))/365.25)

price_bs <- rep(NA, dim(spy_options_call)[1])
price_underlying <- SPY_sub$price[SPY_sub$date == max(SPY_sub$date)]
maturity <- spy_options_call$maturity[1] # same maturity for all
for (k in 1:dim(spy_options_call)[1]) {

  strike <- spy_options_call$strike[k]
  
  # black-scholes call price
  price_bs[k] <- black_scholes_formula(ttm = maturity,
                                       s = price_underlying, 
                                       K = strike, 
                                       r = r, 
                                       sigma = sigma_mle,
                                       dividend = 0)
  
}


plot(price_bs, type = 'l')
points(spy_options_call$price)
# plot((price_bs - spy_options_call$price) / price_bs)



## implied volatility


vol_implied <- rep(NA, dim(spy_options_call)[1])
price_underlying <- SPY_sub$price[SPY_sub$date == max(SPY_sub$date)]
maturity <- spy_options_call$maturity[1] # same maturity for all
for (k in 1:dim(spy_options_call)[1]) {

  price_call_obs <- spy_options_call$price[k]
  
  strike <- spy_options_call$strike[k]

  f <- function(sigma) black_scholes_formula(ttm = maturity,
                                             s = price_underlying,
                                             K = strike,
                                             r = r,
                                             sigma,
                                             dividend = 0.02)

  vol_implied[k] <- uniroot(function(sigma) f(sigma) - price_call_obs, c(-1, 1))$r


}


plot(spy_options_call$strike / price_underlying, vol_implied, type = 'l')


## perfect match when using implied volatility - OF COURSE!

price_bs_implied <- rep(NA, dim(spy_options_call)[1])
for (k in 1:dim(spy_options_call)[1]) {
  
  strike <- spy_options_call$strike[k]
  
  # black-scholes call price
  price_bs_implied[k] <- black_scholes_formula(ttm = maturity,
                                       s = price_underlying, 
                                       K = strike, 
                                       r = r, 
                                       sigma = vol_implied[k],
                                       dividend = 0.02)
  
}

plot(price_bs_implied, type = 'l')
points(spy_options_call$price)

# ## VIX from yahoo using the package quantmod
#
# # get VIX
# setSymbolLookup('^VIX'='yahoo')
# getSymbols('^VIX')
# str(VIX)
#
# # convert data from xts table to data frame
# VIX <- as.data.frame(VIX)
# VIX <- cbind.data.frame(as.Date(row.names(VIX)), VIX$VIX.Close)
# rownames(VIX) <- NULL
# colnames(VIX) <- c('date', 'vol')
# head(VIX)
# str(VIX)
# plot(VIX$vol, type = 'l')



