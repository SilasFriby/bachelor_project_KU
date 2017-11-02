
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


## choose the n maturity dates with most observed option prices and choose options which has a price above x

n <- 4
x <- 100

number_of_options <- sapply(names(price_bs_list), function(i) {
  dim(price_bs_list[[i]][['calls']])[1]
})
number_of_options <- rev(sort(number_of_options))[1:n]
number_of_options_index <- which(names(price_bs_list) %in% names(number_of_options))

max_prices <- sapply(names(price_bs_list), function(i) {
  max(price_bs_list[[i]][['calls']]$price) > x
})
max_prices_index <- which(max_prices)

index_options <- Reduce(intersect, list(max_prices_index, number_of_options_index))# common elements



## plot observed prices and bs prices

par( mfrow = c( 2, 2 ) )
for (i in index_options) {
  
  x_obs <- data_option[[names(data_option)[i]]][['calls']]$Strike
  y_obs <- (data_option[[names(data_option)[i]]][['calls']]$Ask + data_option[[names(data_option)[i]]][['calls']]$Bid) / 2
  y_bs <- price_bs_list[[names(data_option)[i]]][['calls']]$price
  mean_squared_error <- round(sum((y_obs - y_bs) ^ 2) / length(y_obs), 2)
  bid_ask <- round(sum(data_option[[names(data_option)[i]]][['calls']]$Ask - data_option[[names(data_option)[i]]][['calls']]$Bid) / length(y_obs), 2)
  
  plot(x_obs, y_obs, 
       main = paste('maturity =', names(data_option)[i], '\nmean squared error =', mean_squared_error, '\naverage bid-ask =', bid_ask, sep = ' '), 
       xlab = 'strike', 
       ylab = 'price')
  points(x_obs, y_bs, type = 'l')
  
}


## remove options that have not been traded (i.e vol = 0)

par( mfrow = c( 2, 2 ) )
for (i in index_options) {
  
  volume_index <- which(data_option[[names(data_option)[i]]][['calls']]$Vol > 0)
  x_obs <- data_option[[names(data_option)[i]]][['calls']]$Strike[volume_index]
  y_obs <- (data_option[[names(data_option)[i]]][['calls']]$Ask + data_option[[names(data_option)[i]]][['calls']]$Bid)[volume_index] / 2
  y_bs <- price_bs_list[[names(data_option)[i]]][['calls']]$price[volume_index]
  mean_squared_error <- round(sum((y_obs - y_bs) ^ 2) / length(y_obs), 2)
  bid_ask <- round(sum(data_option[[names(data_option)[i]]][['calls']]$Ask - data_option[[names(data_option)[i]]][['calls']]$Bid) / length(y_obs), 2)
  
  plot(x_obs, y_obs, 
       main = paste('maturity =', names(data_option)[i], '\nmean squared error =', mean_squared_error, '\naverage bid-ask =', bid_ask, sep = ' '), 
       xlab = 'strike', 
       ylab = 'price')
  points(x_obs, y_bs, type = 'l')
  
}




#### implied volatility


vol_implied_list <- list() 
for (option_type in c('calls')) {
  
  for (list_date in names(data_option)[index_options]) { # number of different maturity dates
    
    volume_index <- which(data_option[[list_date]][['calls']]$Vol > 0)
    data_temp <- data_option[[list_date]][[option_type]][volume_index, ]
    vol_implied <- rep(NA, dim(data_temp)[1])
    
    for (k in 1:dim(data_temp)[1]) {
      
      f <- function(sigma) black_scholes_formula(ttm = data_temp$ttm[k],
                                                 s = price_underlying, 
                                                 K = data_temp$Strike[k], 
                                                 r = r, 
                                                 sigma,
                                                 dividend = 0)
      
      price_obs <- mean(data_temp$Ask[k], data_temp$Bid[k])
      vol_implied[k] <- uniroot(function(sigma) f(sigma) - price_obs, c(-10, 10))$r
      
    }
    
    vol_implied_list[[list_date]][[option_type]] <- data.frame(vol_implied = vol_implied, strike = data_temp$Strike)
    
  }
  
}


## plot implied volatility


par( mfrow = c( 2, 2 ) )
for (i in index_options) {
  
  volume_index <- which(data_option[[names(data_option)[i]]][['calls']]$Vol > 0)
  x <- data_option[[names(data_option)[i]]][['calls']]$Strike[volume_index] / price_underlying
  y <- vol_implied_list[[names(data_option)[i]]][['calls']]$vol_implied
  
  plot(x, y, 
       main = paste('maturity =', names(data_option)[i], sep = ' '), 
       xlab = 'strike', 
       ylab = 'implied volatility',
       type = 'l')
  abline(h = sigma_mle)
  
}


plot(spy_options_call$strike / price_underlying, vol_implied, type = 'l')
abline(h = sigma_mle)



## perfect match when using implied volatility - OF COURSE! (sanity check)

price_bs_implied_list <- list() 
for (option_type in c('calls')) {
  
  for (list_date in names(data_option)[index_options]) { # number of different maturity dates
    
    volume_index <- which(data_option[[list_date]][['calls']]$Vol > 0)
    data_temp <- data_option[[list_date]][[option_type]][volume_index, ]
    price_bs_implied <- rep(NA, dim(data_temp)[1])
    
    for (k in 1:dim(data_temp)[1]) {

      price_bs_implied[k] <- black_scholes_formula(ttm = data_temp$ttm[k],
                                                   s = price_underlying,
                                                   K = data_temp$Strike[k],
                                                   r = r,
                                                   sigma = vol_implied_list[[list_date]][[option_type]]$vol_implied[k],
                                                   dividend = 0)
      
    }
    
    price_bs_implied_list[[list_date]][[option_type]] <- data.frame(price = price_bs_implied, strike = data_temp$Strike)
    
  }
  
}


## plot implied prices

par( mfrow = c( 2, 2 ) )
for (i in index_options) {
  
  volume_index <- which(data_option[[names(data_option)[i]]][['calls']]$Vol > 0)
  x <- data_option[[names(data_option)[i]]][['calls']]$Strike[volume_index]
  y_obs <- (data_option[[names(data_option)[i]]][['calls']]$Ask + data_option[[names(data_option)[i]]][['calls']]$Bid)[volume_index] / 2
  y_implied <- price_bs_implied_list[[names(data_option)[i]]][['calls']]$price
  
  plot(x, y_obs, 
       main = paste('maturity =', names(data_option)[i], sep = ' '), 
       xlab = 'strike', 
       ylab = 'price')
  points(x, y_implied, type = 'l')
  
}




