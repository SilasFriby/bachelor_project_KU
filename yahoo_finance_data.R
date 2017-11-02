
rm(list = ls())

source('black_scholes_formula.R')
library(quantmod)


## stock prices from yahoo using the package quantmod

#' Title
#'
#' @param symbol - data symbol, e.g '^OEX' for S&P 100 (charater) 
#' @param option_maturity_year - maturity years for options, e.g c('2017', '2018') (character vector)
#' @param data_source - where to extract data from. Works only for 'yahoo'
#'
#' @return list of symbol prices and list of option prices with the symbol as underlying

yahoo_data_fct <- function(symbol, option_maturity_year, data_source = 'yahoo') {
  
  # stock prices
  setSymbolLookup(symbol = data_source)
  data_env <- new.env() # create environment
  getSymbols(symbol, env = data_env)
  data_name <- get(ls(data_env), envir = data_env)
  
  
  
  # convert data from xts table to data frame
  data_stock <- as.data.frame(data_name) 
  data_stock <- cbind.data.frame('date' = as.Date(row.names(data_stock)), data_stock)
  rownames(data_stock) <- NULL
  
  # option prices
  data_option <- getOptionChain(symbol, Exp = option_maturity_year , src = data_source)
  
  # add the two columns maturity date and time to maturity (in years)
  for (i in 1:length(data_option)) {
    
    for (option_type in c('calls', 'puts')) {
      
      data_option[[names(data_option)[i]]][[option_type]] <- cbind.data.frame('option_symbol' = row.names(data_option[[names(data_option)[i]]][[option_type]]),
                                                                              data_option[[names(data_option)[i]]][[option_type]])
      row.names(data_option[[names(data_option)[i]]][[option_type]]) <- NULL
      
      today <- max(data_stock$date)
      data_option[[names(data_option)[i]]][[option_type]]$maturity_date <- as.Date(substring(data_option[[names(data_option)[i]]][[option_type]]$option_symbol, nchar(ls(data_env)) + 1, nchar(ls(data_env)) + 6), format = '%y%m%d')
      data_option[[names(data_option)[i]]][[option_type]]$ttm <- as.numeric((data_option[[names(data_option)[i]]][[option_type]]$maturity_date - today) / 365.25)
      
    }
    
  }
  
  
  return(
    list(
      data_stock = data_stock,
      data_option = data_option
    )
  )
  
} 



#### get S&P 100 data from yahoo


symbol <- '^OEX'
option_maturity_year <- c('2017', '2018', '2019')

data_yahoo <- yahoo_data_fct(symbol = symbol,
                             option_maturity_year = option_maturity_year)



#### test black scholes assumptions on stock prices


## subset data

data_stock <- subset(data_yahoo$data_stock, date >= as.Date('2010-01-01'))
plot(data_stock$date, data_stock$OEX.Close, type = 'l')


## log returns

log_return <- log(data_stock$OEX.Close[-1] / data_stock$OEX.Close[-dim(data_stock)[1]])
hist(log_return)


## test independence of log returns - seems OK

acf(log_return)
Box.test(log_return, lag = 1, type = "Ljung-Box")

## test normal distribution - not OK in general! Heavier tales than a normal distribution

qqnorm(log_return);qqline(log_return, col = 2)
shapiro.test(log_return)




#### estimate Q-parameters in GBM - assuming equidistant data - recall drift of stock is not neccessary for pricing!


## mle

dt <- 1/250 # approximate number of observations per year
sigma_mle <- sqrt(var(log_return) / dt)
mu_mle <- mean(log_return) / dt + sigma_mle ^ 2 / 2  


## mle variance

var_sigma_mle <- sigma_mle ^ 2 /2
var_mu_mle <- sigma_mle ^ 2 * (2 + sigma_mle ^ 2 * dt) / (2 * dt) # the variance on the drift estimate in GBM is very high

## hist

hist(log_return, freq = FALSE)
lines(sort(log_return), dnorm(sort(log_return), mean = (mu_mle - sigma_mle ^ 2 / 2) * dt, sd = sigma_mle * sqrt(dt)), lwd = 2)





#### compute option prices using black scholes and compare to observed prices

data_option <- data_yahoo$data_option

r <- 0 # constant market rate
price_underlying <- tail(data_stock$OEX.Close, 1)
price_bs_list <- list()

for (option_type in c('calls')) {
  
  for (list_date in names(data_option)) { # number of different maturity dates
  
    data_temp <- data_option[[list_date]][[option_type]]
    price_bs_temp <- rep(NA, dim(data_temp)[1])
    
    for (k in 1:dim(data_temp)[1]) {

      # black-scholes option price
      price_bs_temp[k] <- black_scholes_formula(ttm = data_temp$ttm[k],
                                                s = price_underlying, 
                                                K = data_temp$Strike[k], 
                                                r = r, 
                                                sigma = sigma_mle,
                                                dividend = 0)
      
    }
   
    price_bs_list[[list_date]][[option_type]] <- data.frame(price = price_bs_temp, strike = data_temp$Strike)
     
  }
  
}


## choose the n maturity dates with most observed option prices

n <- 9
number_of_options <- sapply(names(price_bs_list), function(i) {
  dim(price_bs_list[[i]][['calls']])[1]
})
number_of_options <- rev(sort(number_of_options))[1:n]
index <- which(names(price_bs_list) %in% names(number_of_options))


## plot observed prices and bs prices

par( mfrow = c( 3, 3 ) )
for (i in index) {
  
  x_obs <- data_option[[names(data_option)[i]]][['calls']]$Strike
  y_obs <- data_option[[names(data_option)[i]]][['calls']]$Ask
  x_bs <- price_bs_list[[names(data_option)[i]]][['calls']]$strike
  y_bs <- price_bs_list[[names(data_option)[i]]][['calls']]$price
  
  plot(x_obs, y_obs, main = names(data_option)[i], xlab = 'strike', ylab = 'price')
  points(x_bs, y_bs, type = 'l')
  
}




# # plot((price_bs - spy_options_call$price) / price_bs)
# 
# 
# 
# ## implied volatility
# 
# 
# vol_implied <- rep(NA, dim(spy_options_call)[1])
# price_underlying <- data_stock$price[data_stock$date == max(data_stock$date)]
# maturity <- spy_options_call$maturity[1] # same maturity for all
# for (k in 1:dim(spy_options_call)[1]) {
# 
#   price_call_obs <- spy_options_call$price[k]
#   
#   strike <- spy_options_call$strike[k]
# 
#   f <- function(sigma) black_scholes_formula(ttm = maturity,
#                                              s = price_underlying,
#                                              K = strike,
#                                              r = r,
#                                              sigma,
#                                              dividend = 0.02)
# 
#   vol_implied[k] <- uniroot(function(sigma) f(sigma) - price_call_obs, c(-1, 1))$r
# 
# 
# }
# 
# 
# plot(spy_options_call$strike / price_underlying, vol_implied, type = 'l')
# abline(h = sigma_mle)
# 
# ## perfect match when using implied volatility - OF COURSE!
# 
# price_bs_implied <- rep(NA, dim(spy_options_call)[1])
# for (k in 1:dim(spy_options_call)[1]) {
#   
#   strike <- spy_options_call$strike[k]
#   
#   # black-scholes call price
#   price_bs_implied[k] <- black_scholes_formula(ttm = maturity,
#                                        s = price_underlying, 
#                                        K = strike, 
#                                        r = r, 
#                                        sigma = vol_implied[k],
#                                        dividend = 0.02)
#   
# }
# 
# plot(price_bs_implied, type = 'l')
# points(spy_options_call$price)





