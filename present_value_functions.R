

pv_european_option_MC <- function(t, time_vector, dt, underlying_paths, strike, short_rate_paths, type) {
  
  index <- which(time_vector == t)
  indexes <- index:dim(short_rate_paths)[1]
  
  # integral short rate over time
  if (length(indexes) > 1) {
    int <- colSums(short_rate_paths[indexes, ]) * dt
  } else {
    int <- 0 # on the last time point
  }
    
  
  
  # discount factor
  discount_factor <- exp(-int)
  expected_discount_factor <- mean(discount_factor)
  
  # payoff from option
  if (type == 'call') payoff <- pmax(underlying_paths[tail(indexes, 1), ] - strike, 0)
  if (type == 'put') payoff <- pmax(strike - underlying_paths[tail(indexes, 1), ], 0)
  expected_payoff <- mean(payoff)
  
  # expected prsent value
  pv <- expected_discount_factor * expected_payoff # only correct if the short rate and stock are assumed to be independent 
  
  return(list(
    pv = pv,
    quantile_1_pct = quantile(discount_factor, 0.01) * quantile(payoff, 0.01), # giver det mening?
    quantile_99_pct = quantile(discount_factor, 0.99) * quantile(payoff, 0.99) # giver det mening?
  ))
  
}


# test <- sapply(time_vector, function(i) {
#   pv_european_option_MC(t = i, 
#                         time_vector = time_vector, 
#                         dt = dt, 
#                         underlying_paths = stock_prices, 
#                         strike = 120, 
#                         short_rate_paths = short_rate, 
#                         type = 'call')
# })


pv_non_callable_bond_MC <- function(t, time_vector, dt, payment_dt, short_rate_paths) {
  
  index <- which(time_vector == t)
  indexes <- index:dim(short_rate_paths)[1]
  
    
  # expected values for ZCB's
  expected_zcb_values <- sapply(indexes, function(i){
    
    temp_indexes <- index:i 
    
    # integral short rate over time for ZCB (t, T_i)
    if(length(temp_indexes) > 1) {
      int <- colSums(short_rate_paths[temp_indexes, ]) * dt
    } else {
      int <- 0 # the ZCB 'today'
    }
    
    
    # discount factor / ZCB value
    discount_factor <- exp(-int)
    expected_discount_factor <- mean(discount_factor)
    
  })
  
  # expected prsent value
  pv <- payment_dt * sum(expected_zcb_values)
    
  
  return(pv)
  
}



test <- rep(NA, length(time_vector))
for( i in time_vector) {
  test[which(time_vector == i)] <- pv_non_callable_bond_MC(t = i,
                                     time_vector = time_vector,
                                     dt = dt,
                                     payment_dt = 100,
                                     short_rate_paths = short_rate)
}

plot(time_vector, test)
